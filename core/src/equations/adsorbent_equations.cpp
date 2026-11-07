/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "adsim/equations/adsorbent_equations.h"
#include "adsim/objects/user_data.h"
#include "adsim/boundary_conditions/step.h"
#include "adsim/math/constants.h"  // kMinDenomPressure/Temperature, detail::DenomFloor

// STL includes
#include <vector>
#include <functional>
#include <algorithm>

// SUNDIALS includes
#include <ida/ida.h>


namespace adsim {

  // Cell-centred independent-equation coefficients (hi, ho, Lam, Dax, hp, Ki)
  // have no physical meaning in void cells (no adsorbent, no packed-bed
  // correlation). Zero them at void cells, then mirror the ghost cells from
  // the first/last *bed* cell so BC residuals reference representative bed
  // values. Bit-identical to the legacy `var[0]=var[1]; var[N-1]=var[N-2]`
  // when no voids are enabled.
  static inline void ZeroVoidCellsAndMirrorGhosts(std::vector<double>& var, const Bed& bed) {
    const int N = static_cast<int>(var.size());
    if (N < 2) return;
    if (bed.bottom_void_.enabled_) {
      const int iv = bed.GetBottomVoidCellIndex();
      if (iv >= 0 && iv < N) var[iv] = 0.0;
    }
    if (bed.top_void_.enabled_) {
      const int iv = bed.GetTopVoidCellIndex();
      if (iv >= 0 && iv < N) var[iv] = 0.0;
    }
    const int first_bed = bed.GetFirstBedCellIndex();
    const int last_bed = bed.GetLastBedCellIndex();
    if (first_bed >= 0 && first_bed < N) var[0] = var[first_bed];
    if (last_bed >= 0 && last_bed < N) var[N - 1] = var[last_bed];
  }

  // The running step's heating medium on the wall's outer face, or nullptr
  // when the step has none (or no step is set): the wall then exchanges heat
  // with Bed::ambient_temperature_ through each layer's external coefficient.
  static inline const HeatingMedium* ActiveHeatingMedium(const Equations& state) {
    const Step* step = state.CurrentStep();
    return step ? step->ActiveHeatingMedium() : nullptr;
  }

  // The running step's internal heat exchanger: its volumetric conductance
  // a_x h_x [W/(m3 K)] and the medium temperature at the residual time. ua is
  // 0 (and the bed balances take no exchanger term at all) without a medium
  // or with a_x h_x = 0.
  struct InternalExchanger {
    double ua = 0.0;
    double T_m = 0.0;
  };
  static inline InternalExchanger ActiveInternalExchanger(const UserData* udata) {
    InternalExchanger x;
    const HeatingMedium* medium = ActiveHeatingMedium(udata->state);
    if (!medium) return x;
    const double ua = medium->ExchangerConductance();
    if (!(ua > 0.0)) return x;
    x.ua = ua;
    x.T_m = medium->Temperature(udata->residual_time_ - udata->step_start_time_);
    return x;
  }



  namespace {
    // Smoothed upwind face value: the left value for u >> 0, the right value
    // for u << 0, blended through stagnant flow (see kUpwindSmoothingVelocity).
    // Evaluates only the candidate(s) that carry weight.
    template <class Left, class Right>
    double UpwindFace(double u, Left left, Right right) {
      const double w = detail::UpwindWeight(u);
      if (w == 1.0) return left();
      if (w == 0.0) return right();
      return w * left() + (1.0 - w) * right();
    }

    // Danckwerts boundary faces for the fluid energy balance. The boundary
    // conditions set the fluid-temperature ghost cell so that the FACE
    // temperature is the boundary value (T_ghost + T_first = 2 T_in at an
    // inlet; T_ghost = T_first at an outlet or closed end). Across a boundary
    // face of the bed the enthalpy flux is the convective rho Cp u T_face
    // alone, with no conduction (the residuals drop Lam dT/dz there), so the
    // upstream value of a boundary face is the face temperature, not the ghost
    // value: upwinding to the ghost, 2 T_in - T_1, would add
    // rho Cp u (T_in - T_1) to the feed's enthalpy. For a face of an interior
    // cell that is not on the boundary this is the neighbour's temperature.
    // side = -1: the west face of cell i, +1: the east face.
    inline double BoundaryFaceUpstreamT(const std::vector<double>& T, int i, int side) {
      const int last = static_cast<int>(T.size()) - 1;
      const int nb = i + side;
      if (nb == 0 || nb == last) return 0.5 * (T[nb] + T[i]);
      return T[nb];
    }

    // Convective face temperatures of the energy balances, one per face
    // (face n between cells n and n + 1), reconstructed with the advection
    // scheme of the mole fractions. Defined after the reconstructions below.
    void ConvectiveFaceTemperatures(const UserData* udata, const std::vector<double>& T,
                                    std::vector<double>& T_face);

    // Heat released per mole adsorbed, as it enters the energy balances:
    // dH for the ISOSTERIC convention (dH is the isosteric heat q_st);
    // dH - Cpa*T for the INTERNAL_ENERGY convention (Haghpanah/Ward), which
    // adds the adsorbed-phase sensible-enthalpy transfer term; dH + R*T for
    // the ISOSTERIC_RT convention (dH is -dU, q_st = -dU + R*T). T is the
    // temperature of the adsorbed phase (fluid T in LTE, solid Ts in LTNE).
    // See HeatOfAdsorptionConvention in fluid.h. Shared by the energy
    // residuals and the logged component energy source so the two agree.
    // Axial position of the centre of cell n [m], from the bottom bed-end face
    // (the coordinate of DomainVectors::positions): where the manufactured
    // sources of the cell balances are evaluated.
    inline double CellCentre(const DomainVectors& dv, int n) {
      return dv.positions[n - 1] + 0.5 * dv.cell_width[n];
    }

    // The manufactured source of `component` in `map` (verification hook,
    // see manufactured_sources.h), or null when there is none.
    inline const ManufacturedSources::Field* SourceOf(
        const std::unordered_map<std::string, ManufacturedSources::Field>& map, const std::string& component) {
      if (map.empty()) return nullptr;
      auto it = map.find(component);
      return (it != map.end() && it->second) ? &it->second : nullptr;
    }

    inline double AdsorptionHeat(double dH, double Cpa, double T, HeatOfAdsorptionConvention convention) {
      switch (convention) {
        case HeatOfAdsorptionConvention::INTERNAL_ENERGY: return dH - Cpa * T;
        case HeatOfAdsorptionConvention::ISOSTERIC_RT: return dH + constants::kIdealGasConstant * T;
        case HeatOfAdsorptionConvention::ISOSTERIC: break;
      }
      return dH;
    }
  }

  void UpdateConcentration(UserData* udata) {
    auto& bed = udata->bed;
    auto& state = udata->state;

    auto& concentration = state.concentration;
    auto& C = concentration.var;
    auto& res_C = concentration.res_var;
    const auto& pressure = state.pressure;
    const auto& P = pressure.var;
    const auto& temperature = state.fluid_temperature;
    const auto& T = temperature.var;

    const double R = constants::kIdealGasConstant;
    int N = C.size();
    for (int i = 0; i < N; ++i) {
      // Floor T in the ideal-gas denominator so a transient non-positive
      // temperature cannot produce Inf/NaN (no-op for physical T).
      const double Tc = detail::DenomFloor(T[i], constants::kMinDenomTemperature);
      const double conc = P[i] / (R * Tc);
      res_C[i] = C[i] - conc;
    }
  }

  bool ErgunFaceCoefficients(const Equations& state, const Bed* bed, int face, double& A, double& B) {
    const auto& dv = state.domain_vectors;
    const auto& dx = dv.cell_width;
    const auto& K = dv.permeability;
    const auto& beta = dv.inertial_coefficient;
    const auto& is_void = dv.is_void_cell;
    const auto& mu = state.viscosity.var;
    const auto& rho = state.fluid_density.var;
    const int i = face;
    // A state without the momentum fields (or a face outside them) has no
    // Ergun face to describe.
    const int n = static_cast<int>(beta.size());
    if (i < 0 || i + 1 >= n || (int)dx.size() < n || (int)K.size() < n
        || (int)mu.size() < n || (int)rho.size() < n) return false;
    const bool void_left  = (i     < (int)is_void.size()) && is_void[i];
    const bool void_right = (i + 1 < (int)is_void.size()) && is_void[i + 1];

    // Transport properties (viscosity, density) of each side. A ghost cell's
    // gas state is the boundary condition's mirror construct (P = 2 Pf - P1,
    // T = 2 Tin - T1), not a gas that exists: when a step switch moves the
    // boundary pressure, the ghost pressure and density jump by the full swing
    // (to ~2 bar for a 1 bar feed into an evacuated bed), and the inertial
    // coefficient's 1 / sqrt(rho) sensitivity to that jump made IDA's Newton
    // corrector diverge at the start of every pressurisation step (the
    // Jacobian is evaluated at the predicted, pre-switch state). The two
    // boundary faces therefore take their transport properties from the
    // interior neighbour, whose state is continuous through the switch. At
    // steady flow the difference from the ghost-averaged value is the
    // density change over half a cell (well below 0.1 %).
    const int n_cells = static_cast<int>(is_void.size());
    const int tl = (i == 0) ? 1 : i;
    const int tr = (i + 1 == n_cells - 1) ? n_cells - 2 : i + 1;

    if (!void_left && !void_right) {
      if (beta[i] == 0.0 && beta[i + 1] == 0.0) return false;
      // The two half-cells between the cell centres act in series, each with
      // its own viscous (mu dx / 2K) and inertial (beta rho dx / 2)
      // resistance. With beta = 0 this is exactly the harmonic-K Darcy face,
      // so a Darcy layer next to an Ergun layer is treated consistently.
      A = 0.5 * (mu[tl] * dx[i] / K[i] + mu[tr] * dx[i + 1] / K[i + 1]);
      B = 0.5 * (beta[i] * rho[tl] * dx[i] + beta[i + 1] * rho[tr] * dx[i + 1]);
      return true;
    }

    // Exactly one side is a void cell: the bed-side cell's resistance over its
    // full width (Darcy's law does not apply inside a lumped tank).
    const int bed_idx = void_left ? (i + 1) : i;
    if (beta[bed_idx] == 0.0) return false;
    const int t_idx = void_left ? tr : tl;  // the bed side, or its interior neighbour at a boundary face
    double R_face = mu[t_idx] * dx[bed_idx] / K[bed_idx];
    // Optional user-specified void resistance: applied ONLY at the void-bed
    // interface (where a screen or distributor would sit), not at the
    // ghost-void face (the BC face, where the boundary condition controls the
    // flow). It is a linear Pa.s/m term, so it adds no inertial resistance.
    const int void_idx = void_left ? i : (i + 1);
    const int N_cells = (int)is_void.size();
    const bool bed_is_ghost = (bed_idx == 0) || (bed_idx == N_cells - 1);
    if (!bed_is_ghost && bed) {
      if (const VoidRegion* vr = bed->VoidRegionAt(void_idx)) {
        const double R_user = vr->resistance_coefficient_.Get();
        if (R_user > 0.0) R_face += R_user;
      }
    }
    A = R_face;
    B = beta[bed_idx] * rho[t_idx] * dx[bed_idx];
    return true;
  }

  void UpdateVelocity(UserData* udata) {
    auto& bed = udata->bed;
    auto& state = udata->state;
    auto& domain_vectors = state.domain_vectors;

    auto& velocity = state.velocity;
    auto& u = velocity.var;
    auto& res_u = velocity.res_var;
    auto& pressure = state.pressure;
    auto& P = pressure.var;

    const auto& dx = domain_vectors.cell_width;
    auto& K = domain_vectors.permeability;
    const auto& mu = state.viscosity.var;
    const auto& is_void = domain_vectors.is_void_cell;

    // Manufactured momentum source [Pa/m] (verification hook; empty in a
    // run): -dP/dz = (mu/K) u + beta rho |u| u - S, i.e. S times the distance
    // between the two cell centres is added to the face's pressure drop.
    const auto& mms_momentum = udata->manufactured_sources_.momentum;
    const double t_mms = udata->residual_time_;
    auto mms_face_drop = [&](int i) {
      return mms_momentum(domain_vectors.positions[i], t_mms) * 0.5 * (dx[i] + dx[i + 1]);
    };

    int N = u.size();
    for (int i = 0; i < N; ++i) {
      // Ergun faces (an inertial resistance on either side) carry the implicit
      // residual u - (g - B|u|u)/A (see ErgunMomentumResidual). A pressure or
      // closed boundary condition replaces it on its own boundary face (see
      // PressureBoundaryCondition / ClosedBoundaryCondition::Apply).
      double A = 0.0;
      double B = 0.0;
      if (ErgunFaceCoefficients(state, &bed, i, A, B)) {
        double g = P[i] - P[i + 1];
        if (mms_momentum) g += mms_face_drop(i);
        res_u[i] = ErgunMomentumResidual(u[i], g, A, B);
        continue;
      }

      // Darcy faces (Kozeny-Carman or Blake-Kozeny, no inertial term): the
      // legacy closed forms, verbatim, so Kozeny-Carman decks reproduce their
      // previous results bit for bit.
      const bool void_left  = (i     < (int)is_void.size()) && is_void[i];
      const bool void_right = (i + 1 < (int)is_void.size()) && is_void[i + 1];

      if (!void_left && !void_right) {
        // Standard porous-media face: harmonic-K averaging.
        double dP_face = P[i + 1] - P[i];
        if (mms_momentum) dP_face -= mms_face_drop(i);
        const double vel = -2 * K[i] * K[i + 1] * dP_face
              / (K[i] * mu[i + 1] * dx[i + 1] + K[i + 1] * mu[i] * dx[i]);
        res_u[i] = u[i] - vel;
      } else {
        // Exactly one side is a void cell. Darcy's law does not apply inside
        // a lumped tank — use the bed-side cell's one-sided Darcy as the
        // resistance at the face. The void contributes no resistance by
        // default.
        const int bed_idx = void_left ? (i + 1) : i;
        const int void_idx = void_left ? i : (i + 1);
        double R_face = mu[bed_idx] * dx[bed_idx] / K[bed_idx];

        // Optional user-specified void resistance: applied ONLY at the
        // void↔bed interface (the manifold-to-bed face, where a screen or
        // distributor would physically sit), not at the ghost↔void interface
        // (the BC face, where the FlowBC already controls the inlet flow).
        const int N_cells = (int)is_void.size();
        const bool bed_is_ghost = (bed_idx == 0) || (bed_idx == N_cells - 1);
        if (!bed_is_ghost) {
          if (const VoidRegion* vr = bed.VoidRegionAt(void_idx)) {
            const double R_user = vr->resistance_coefficient_.Get();
            if (R_user > 0.0) R_face += R_user;
          }
        }

        const double vel = -(P[i + 1] - P[i]) / R_face;
        res_u[i] = u[i] - vel;
      }
    }
  }

  void UpdatePressure(UserData* udata) {
    auto& fluid = udata->fluid;
    auto& bed = udata->bed;
    auto& state = udata->state;
    auto& domain_vectors = state.domain_vectors;

    auto& pressure = state.pressure;
    auto& P = pressure.var;
    auto& ddt_P = pressure.ddt_var;
    auto& res_P = pressure.res_var;
    const auto& concentration = state.concentration;
    const auto& C = concentration.var;
    const auto& ddt_C = concentration.ddt_var;
    const auto& res_C = concentration.res_var;
    const auto& velocity = state.velocity;
    const auto& u = velocity.var;
    const auto& temperature = state.fluid_temperature;
    const auto& T = temperature.var;
    const auto& ddt_T = temperature.ddt_var;
    const auto& vis = state.viscosity.var;
    
    const auto& dx = domain_vectors.cell_width;
    const auto& A = domain_vectors.area;
    const auto& dV = domain_vectors.volume;
    const auto& K = domain_vectors.permeability;
    const auto& rho_b = domain_vectors.bulk_density;
    const auto& et = domain_vectors.et;

    // Manufactured overall-balance source [mol/(m3 s)] (verification hook).
    const auto& mms_total = udata->manufactured_sources_.total;

    int N = P.size(); // Includes ghost cells
    for (int n = 1; n < N - 1; ++n) { // Interior
      // Floor every pressure/temperature that appears as a denominator so a
      // transient non-positive state the BDF solver iterates through cannot
      // spawn Inf/NaN. No-op for physical states (floor << any real P,T).
      const double Pc = detail::DenomFloor(P[n], constants::kMinDenomPressure);
      const double Tc = detail::DenomFloor(T[n], constants::kMinDenomTemperature);

      double Te = 0.5 * (T[n] + T[n + 1]);
      double Tw = 0.5 * (T[n] + T[n - 1]);
      Te = detail::DenomFloor(Te, constants::kMinDenomTemperature);
      Tw = detail::DenomFloor(Tw, constants::kMinDenomTemperature);
      double dT = dx[n] * et[n] * ddt_T[n] / Tc;

      double Pe = 0.5 * (P[n] + P[n + 1]);
      double Pw = 0.5 * (P[n] + P[n - 1]);
      double dP = dx[n] * et[n] * ddt_P[n] / Pc;

      double u_e = u[n];
      double u_w = u[n - 1];

      double F_e = (Tc / Pc) * (Pe * u_e / Te);
      double F_w = (Tc / Pc) * (Pw * u_w / Tw);

      // Total adsorption source: sum over every adsorbent present and every
      // component, each loading rate weighted by that adsorbent's bulk density.
      double rho_b_ddt_q = 0.0;
      for (const auto& ca : domain_vectors.cell_adsorbents[n]) {
        for (const auto* ddt_q : state.GetLoadingRateTerms(ca.id)) {
          rho_b_ddt_q += ca.bulk_density * (*ddt_q)[n];
        }
      }

      double S = R * T[n] * dx[n] * rho_b_ddt_q / Pc;

      res_P[n] = 
        dP - dT 
        + (F_e - F_w) 
        + S;

      // The residual is the balance times R T dx / P, as S above.
      if (mms_total)
        res_P[n] -= R * T[n] * dx[n] * mms_total(CellCentre(domain_vectors, n), udata->residual_time_) / Pc;
    }
  }

  void UpdateTemperature(UserData* udata) {
    if (udata->bed.equations_.thermal_equilibrium_mode_ == ThermalEquilibriumMode::NON_LOCAL) {
      UpdateFluidTemperature(udata);
      return;
    }

    auto& state = udata->state;
    const auto& bed = udata->bed;
    const auto& fluid = udata->fluid;

    auto& temperature = state.fluid_temperature;
    auto& res_T = temperature.res_var;
    const auto& T = temperature.var;
    const auto& ddt_T = temperature.ddt_var;

    const auto& wall_temperature = state.wall_temperature;
    const auto& Tw = wall_temperature.var;

    const auto& pressure = state.pressure;
    const auto& P = pressure.var;
    const auto& ddt_P = pressure.ddt_var;

    const auto& velocity = state.velocity;
    const auto& u = velocity.var;

    const auto& density = state.fluid_density;
    const auto& rho = density.var;
    
    const auto& heat_capacity = state.heat_capacity;
    const auto& Cpg = heat_capacity.var;

    const auto& axial_thermal_conductivity = state.axial_thermal_conductivity;
    const auto& Lam = axial_thermal_conductivity.var;

    const auto& internal_htc = state.internal_heat_transfer_coefficient;
    const auto& hi = internal_htc.var;

    const auto& Cps = state.domain_vectors.Cps;
    const auto& et = state.domain_vectors.et;
    const auto& dV = state.domain_vectors.volume;
    const auto& dx = state.domain_vectors.cell_width;
    const auto& rho_b = state.domain_vectors.bulk_density;
    const auto& A = state.domain_vectors.area;

    const auto& cell_ads = state.domain_vectors.cell_adsorbents;
    const auto& di = state.domain_vectors.di;
    const auto& Qe = state.domain_vectors.energy_source_power;

    int N = T.size();

    // Convective face temperatures (the advection scheme's reconstruction).
    thread_local std::vector<double> T_face;
    ConvectiveFaceTemperatures(udata, T, T_face);

    // heat = dH (ISOSTERIC) or dH - Cpa*T (INTERNAL_ENERGY); see AdsorptionHeat.
    const HeatOfAdsorptionConvention hoa_convention = fluid.heat_of_adsorption_convention_;

    // Manufactured energy source [W/m3] (verification hook).
    const auto& mms_energy = udata->manufactured_sources_.energy;

    // The step's internal heat exchanger (medium -> bed), in bed cells only.
    const InternalExchanger exchanger = ActiveInternalExchanger(udata);
    const auto& is_void = state.domain_vectors.is_void_cell;

    for (int i = 1; i < N - 1; ++i) {
      // Adsorbed-phase heat capacity contribution, per adsorbent:
      // sum_k rho_b,k * sum_comp Cpa_comp * q_{k,comp}.
      double ads_heat_cap = 0.0;
      for (const auto& ca : cell_ads[i]) {
        for (const auto& t : state.GetEnergyTerms(ca.id)) {
          ads_heat_cap += ca.bulk_density * t.heat_capacity * (*t.loading)[i];
        }
      }
      double Cbed = et[i] * rho[i] * Cpg[i] + rho_b[i] * Cps[i] + ads_heat_cap;

      double u_e = u[i];
      double u_w = u[i - 1];

      double T_e = T_face[i];
      double T_w = T_face[i - 1];

      double u_avr = 0.5 * (u_e + u_w);

      double F_e = rho[i] * Cpg[i] * u_avr * T_e;
      double F_w = rho[i] * Cpg[i] * u_avr * T_w;

      double work = dV[i] * ddt_P[i] * et[i];

      // Adsorption energy source: sum_k sum_comp heat_{k,comp} * rho_b,k * dq/dt.
      double Se = 0.0;
      for (const auto& ca : cell_ads[i]) {
        for (const auto& t : state.GetEnergyTerms(ca.id)) {
          const double heat = AdsorptionHeat((*t.heat_of_adsorption)[i], t.heat_capacity, T[i], hoa_convention);
          Se += heat * ca.bulk_density * (*t.loading_rate)[i];
        }
      }
      Se *= dV[i];
      double Acc = dV[i] * Cbed * ddt_T[i];

      double Hi = 4 * hi[i] * dV[i] / di[i];

      double Lam_e = 0.5 * (Lam[i] + Lam[i + 1]);
      double Lam_w = 0.5 * (Lam[i] + Lam[i - 1]);
      
      double fe = dx[i + 1] / (dx[i] + dx[i + 1]);
      double fw = dx[i] / (dx[i - 1] + dx[i]);

      // double Lam_e = 1 / ((1 - fe) / Lam[i] + fe / Lam[i + 1]);
      // double Lam_w = 1 / ((1 - fw) / Lam[i - 1] + fw / Lam[i]);

      // No conduction through the bed's boundary faces (Danckwerts; see
      // BoundaryFaceUpstreamT).
      double D_e = (i == N - 2) ? 0.0 : Lam_e * (T[i + 1] - T[i]) / dx[i];
      double D_w = (i == 1) ? 0.0 : Lam_w * (T[i] - T[i - 1]) / dx[i];

      res_T[i] =
        Acc
        + (F_e - F_w) * A[i]
        - (D_e - D_w) * A[i]
        + Hi * (T[i] - Tw[i])
        - work
        - Se
        - Qe[i] * dV[i];

      // A volumetric source, as the heater's power.
      if (mms_energy)
        res_T[i] -= dV[i] * mms_energy(CellCentre(state.domain_vectors, i), udata->residual_time_);

      // Internal exchanger: a_x h_x (T_m - T) per bed volume.
      if (exchanger.ua > 0.0 && !(i < (int)is_void.size() && is_void[i]))
        res_T[i] -= exchanger.ua * dV[i] * (exchanger.T_m - T[i]);
    }
  }

  void UpdateFluidTemperature(UserData* udata) {
    auto& state = udata->state;
    const auto& fluid = udata->fluid;

    auto& temperature = state.fluid_temperature;
    auto& res_T = temperature.res_var;
    const auto& T = temperature.var;
    const auto& ddt_T = temperature.ddt_var;

    const auto& Ts = state.solid_temperature.var;

    const auto& wall_temperature = state.wall_temperature;
    const auto& Tw = wall_temperature.var;

    const auto& pressure = state.pressure;
    const auto& ddt_P = pressure.ddt_var;

    const auto& velocity = state.velocity;
    const auto& u = velocity.var;

    const auto& density = state.fluid_density;
    const auto& rho = density.var;

    const auto& heat_capacity = state.heat_capacity;
    const auto& Cpg = heat_capacity.var;

    const auto& axial_thermal_conductivity = state.axial_thermal_conductivity;
    const auto& Lam = axial_thermal_conductivity.var;

    const auto& internal_htc = state.internal_heat_transfer_coefficient;
    const auto& hi = internal_htc.var;

    const auto& particle_htc = state.particle_heat_transfer_coefficient;
    const auto& hp = particle_htc.var;

    const auto& et = state.domain_vectors.et;
    const auto& eb = state.domain_vectors.eb;
    const auto& wall_gas_share = state.domain_vectors.wall_gas_share;
    const auto& dV = state.domain_vectors.volume;
    const auto& dx = state.domain_vectors.cell_width;
    const auto& A = state.domain_vectors.area;
    const auto& di = state.domain_vectors.di;
    const auto& Qe = state.domain_vectors.energy_source_power;

    int N = T.size();

    // Convective face temperatures (the advection scheme's reconstruction).
    thread_local std::vector<double> T_face;
    ConvectiveFaceTemperatures(udata, T, T_face);

    // The step's internal heat exchanger: the gas takes the eb share of its
    // heat, the solid the (1 - eb) share (UpdateSolidTemperature), as for the
    // wall. Bed cells only.
    const InternalExchanger exchanger = ActiveInternalExchanger(udata);
    const auto& is_void = state.domain_vectors.is_void_cell;

    for (int i = 1; i < N - 1; ++i) {
      // Fluid-only heat capacity (no solid contribution)
      double Cfluid = et[i] * rho[i] * Cpg[i];

      double u_e = u[i];
      double u_w = u[i - 1];
      double T_e = T_face[i];
      double T_w = T_face[i - 1];
      double u_avr = 0.5 * (u_e + u_w);

      double F_e = rho[i] * Cpg[i] * u_avr * T_e;
      double F_w = rho[i] * Cpg[i] * u_avr * T_w;

      double work = dV[i] * ddt_P[i] * et[i];
      double Acc = dV[i] * Cfluid * ddt_T[i];

      // Bed-to-wall exchange over the internal wall area pi*di*dx (= 4*dV/di).
      // In LTNE the wall contact is split: the gas carries the share
      // wall_gas_share (eb in a packed bed, 0 in a structured layer, whose gas
      // reaches the column wall only through the channel walls), the solid the
      // rest (see UpdateSolidTemperature and the NON_LOCAL branch of
      // UpdateWallTemperature, which absorbs the sum).
      double Hi = 4 * hi[i] * dV[i] / di[i];
      double Hw = wall_gas_share[i] * Hi * (T[i] - Tw[i]);

      // Gas-to-solid heat transfer (hp is volumetric: W/(m3.K))
      double Hp = hp[i] * dV[i] * (T[i] - Ts[i]);

      double Lam_e = 0.5 * (Lam[i] + Lam[i + 1]);
      double Lam_w = 0.5 * (Lam[i] + Lam[i - 1]);

      // No conduction through the bed's boundary faces (Danckwerts; see
      // BoundaryFaceUpstreamT).
      double D_e = (i == N - 2) ? 0.0 : Lam_e * (T[i + 1] - T[i]) / dx[i];
      double D_w = (i == 1) ? 0.0 : Lam_w * (T[i] - T[i - 1]) / dx[i];

      res_T[i] =
        Acc
        + (F_e - F_w) * A[i]
        - (D_e - D_w) * A[i]
        + Hw
        + Hp
        - work
        - Qe[i] * dV[i];

      if (exchanger.ua > 0.0 && !(i < (int)is_void.size() && is_void[i]))
        res_T[i] -= eb[i] * exchanger.ua * dV[i] * (exchanger.T_m - T[i]);
    }
  }

  void UpdateSolidTemperature(UserData* udata) {
    auto& state = udata->state;
    const auto& bed = udata->bed;
    const auto& fluid = udata->fluid;
    const auto thermal_mode = bed.equations_.thermal_equilibrium_mode_;

    auto& solid_temperature = state.solid_temperature;
    auto& res_Ts = solid_temperature.res_var;
    const auto& Ts = solid_temperature.var;
    const auto& Tf = state.fluid_temperature.var;

    int N = Ts.size();

    // LTE mode: algebraic constraint Ts = Tf
    if (thermal_mode == ThermalEquilibriumMode::LOCAL) {
      for (int i = 0; i < N; ++i) {
        res_Ts[i] = Ts[i] - Tf[i];
      }
      return;
    }

    // LTNE mode: full solid energy balance
    const auto& ddt_Ts = solid_temperature.ddt_var;
    const auto& particle_htc = state.particle_heat_transfer_coefficient;
    const auto& hp = particle_htc.var;

    const auto& internal_htc = state.internal_heat_transfer_coefficient;
    const auto& hi = internal_htc.var;

    const auto& Tw = state.wall_temperature.var;

    const auto& Cps = state.domain_vectors.Cps;
    const auto& ks = state.domain_vectors.ks;
    const auto& dV = state.domain_vectors.volume;
    const auto& dx = state.domain_vectors.cell_width;
    const auto& A = state.domain_vectors.area;
    const auto& rho_b = state.domain_vectors.bulk_density;
    const auto& eb = state.domain_vectors.eb;
    const auto& wall_gas_share = state.domain_vectors.wall_gas_share;
    const auto& lambda_s = state.domain_vectors.solid_conductivity;
    const auto& di = state.domain_vectors.di;
    const auto& cell_ads = state.domain_vectors.cell_adsorbents;

    // Adsorption energy source: heat = dH (ISOSTERIC) or dH - Cpa*Ts
    // (INTERNAL_ENERGY, Haghpanah/Ward). See HeatOfAdsorptionConvention in fluid.h.
    // In LTNE the adsorbed phase is at the solid temperature Ts.
    const HeatOfAdsorptionConvention hoa_convention = fluid.heat_of_adsorption_convention_;

    // Solid share (1 - eb) of the step's internal heat exchanger; the gas
    // takes the eb share in UpdateFluidTemperature.
    const InternalExchanger exchanger = ActiveInternalExchanger(udata);

    const auto& is_void = state.domain_vectors.is_void_cell;
    const int first_bed = bed.GetFirstBedCellIndex();
    const int last_bed = bed.GetLastBedCellIndex();

    // Axial conduction through the solid of a structured layer, whose channel
    // walls are continuous along the column (packed and void cells have
    // lambda_s = 0). The conductance of the face between cells j and j + 1 is
    // that of the two half cells in series,
    //   G = 1 / (dx_j / (2 lambda_j) + dx_{j+1} / (2 lambda_{j+1})),
    // which on a uniform grid is the harmonic face mean of lambda over dx; it
    // is 0 on any face touching a packed or void cell, and the bed-end faces
    // carry no flux. A bed without a structured layer skips the term, so its
    // residual is unchanged.
    bool any_solid_conduction = false;
    for (int i = first_bed; i <= last_bed && i < N; ++i) any_solid_conduction = any_solid_conduction || lambda_s[i] > 0.0;
    auto solid_face_conductance = [&](int j) {  // face between cells j and j + 1 [W/(m2-K)]
      if (j < first_bed || j + 1 > last_bed) return 0.0;
      const double lw = lambda_s[j], le = lambda_s[j + 1];
      if (!(lw > 0.0) || !(le > 0.0)) return 0.0;
      return 1.0 / (0.5 * dx[j] / lw + 0.5 * dx[j + 1] / le);
    };

    for (int i = 1; i < N - 1; ++i) {
      if (i < (int)is_void.size() && is_void[i]) continue;  // void cells: Ts is not part of the state vector

      // Solid heat capacity: rho_b*Cps plus the adsorbed-phase contribution,
      // summed per adsorbent: sum_k rho_b,k * sum_comp Cpa_comp * q_{k,comp}.
      double ads_heat_cap = 0.0;
      double Se = 0.0;
      for (const auto& ca : cell_ads[i]) {
        for (const auto& t : state.GetEnergyTerms(ca.id)) {
          ads_heat_cap += ca.bulk_density * t.heat_capacity * (*t.loading)[i];
          const double heat = AdsorptionHeat((*t.heat_of_adsorption)[i], t.heat_capacity, Ts[i], hoa_convention);
          Se += heat * ca.bulk_density * (*t.loading_rate)[i];
        }
      }
      double Csolid = rho_b[i] * Cps[i] + ads_heat_cap;
      Se *= dV[i];

      // Gas-to-solid heat transfer (hp is volumetric: W/(m3.K))
      double Hp = hp[i] * dV[i] * (Tf[i] - Ts[i]);

      double Acc = dV[i] * Csolid * ddt_Ts[i];

      // Solid share of the bed-to-wall exchange (1 - eb packed, all of it in a
      // structured layer); the gas carries the rest in UpdateFluidTemperature,
      // and the wall absorbs the sum.
      double Hi = 4 * hi[i] * dV[i] / di[i];
      double Hw = (1.0 - wall_gas_share[i]) * Hi * (Ts[i] - Tw[i]);

      res_Ts[i] =
        Acc
        - Hp
        + Hw
        - Se;

      if (any_solid_conduction) {
        const double q_e = solid_face_conductance(i) * (Ts[i + 1] - Ts[i]);      // [W/m2], into cell i
        const double q_w = solid_face_conductance(i - 1) * (Ts[i] - Ts[i - 1]);  // [W/m2], out of cell i
        res_Ts[i] -= (q_e - q_w) * A[i];
      }

      if (exchanger.ua > 0.0)
        res_Ts[i] -= (1.0 - eb[i]) * exchanger.ua * dV[i] * (exchanger.T_m - Ts[i]);
    }

    // Ghost cells mirror the first/last *bed* cell. When voids are disabled
    // first_bed = 1 and last_bed = N-2 so this matches the legacy behaviour.
    res_Ts[0] = Ts[0] - Ts[first_bed];
    res_Ts[N - 1] = Ts[N - 1] - Ts[last_bed];
  }

  void UpdateParticleHeatTransferCoefficient(UserData* udata) {
    auto& state = udata->state;
    auto& bed = udata->bed;

    auto& model = bed.equations_.particle_heat_transfer_model_;
    if (model) {
      for (auto& layer : bed.layers_.layers_) {
        model->Update(state, layer);
      }
    }

    ZeroVoidCellsAndMirrorGhosts(state.particle_heat_transfer_coefficient.var, bed);
  }

  void UpdateWallTemperature(UserData* udata) {
    auto& state = udata->state;
    const auto& bed = udata->bed;
    const auto thermal_mode = bed.equations_.thermal_equilibrium_mode_;

    const auto& Tf = state.fluid_temperature.var;
    const auto& Ts = state.solid_temperature.var;

    auto& wall_temperature = state.wall_temperature;
    auto& res_Tw = wall_temperature.res_var;
    const auto& Tw = wall_temperature.var;
    const auto& ddt_Tw = wall_temperature.ddt_var;

    const auto& internal_htc = state.internal_heat_transfer_coefficient;
    const auto& hi = internal_htc.var;

    const auto& external_htc = state.external_heat_transfer_coefficient;
    const auto& ho = external_htc.var;

    const auto& dV = state.domain_vectors.volume;
    const auto& dx = state.domain_vectors.cell_width;
    const auto& A = state.domain_vectors.area;
    const auto& wall_gas_share = state.domain_vectors.wall_gas_share;
    const auto& kw = state.domain_vectors.kw;
    const auto& Cpw = state.domain_vectors.Cpw;
    const auto& rho_w = state.domain_vectors.rho_w;

    const auto& d_i = state.domain_vectors.di;
    const auto& d_o = state.domain_vectors.de;
    // Outer-face sink: the step's heating medium (jacket) when it has one,
    // else the surroundings. The medium's h_m is already in ho (see
    // UpdateExternalHeatTransferCoefficient).
    double Tamb = bed.ambient_temperature_;
    if (const HeatingMedium* medium = ActiveHeatingMedium(state)) {
      Tamb = medium->Temperature(udata->residual_time_ - udata->step_start_time_);
    }

    int N = Tw.size();

    // Interior cells only (exclude ghost cells)
    for (int i = 1; i < N - 1; ++i) {
      const double Acc = rho_w[i] * Cpw[i] * ddt_Tw[i]; // [W/m3]

      // Every term in this residual is per unit WALL volume, so the flux terms
      // are (transfer area / wall cross-sectional area) * h * dT. For an annulus
      // of length dx the wall volume is (pi/4)(do^2 - di^2)dx, so
      //   inside:  (pi di dx) hi dT / ((pi/4)(do^2 - di^2) dx) = 4 di hi dT / (do^2 - di^2)
      //   outside: (pi do dx) ho dT / ((pi/4)(do^2 - di^2) dx) = 4 do ho dT / (do^2 - di^2)
      // Both fluxes must share this denominator, otherwise the heat the wall
      // absorbs does not equal the heat the bed gives up (see
      // ClosedBedWallAndAmbientEnergyBalanceCloses in test/validation).
      // Floor the annulus cross-section so a vanishingly thin wall (d_o -> d_i)
      // cannot drive this denominator to zero and spawn Inf/NaN. No-op for any
      // physical wall thickness (floor << real d_o^2 - d_i^2). See constants.h.
      const double denominator =
          detail::DenomFloor(d_o[i] * d_o[i] - d_i[i] * d_i[i],
                             constants::kMinDenomWallArea);

      // Internal heat flux: in LTNE (local thermal non-equilibrium), split between fluid (eb) and solid (1-eb)
      double flux_i;
      if (thermal_mode == ThermalEquilibriumMode::NON_LOCAL) {
        const double g = wall_gas_share[i];  // eb packed, 0 structured (UpdateFluidTemperature)
        flux_i = 4 * d_i[i] * hi[i] * (g * (Tf[i] - Tw[i]) + (1.0 - g) * (Ts[i] - Tw[i])) / denominator;
      } else {
        flux_i = 4 * d_i[i] * hi[i] * (Tf[i] - Tw[i]) / denominator;
      }

      const double flux_o = 4 * d_o[i] * ho[i] * (Tw[i] - Tamb) / denominator;

      const double conduction = kw[i] / (dx[i] * dx[i]) * (Tw[i + 1] - 2 * Tw[i] + Tw[i - 1]);

      res_Tw[i] = Acc - flux_i + flux_o - conduction;
    }

    bed.equations_.wall_temperature_boundary_condition_->Update(state);
  }

  void UpdateMolarConcentration(UserData* udata, Component component) {
    auto& state = udata->state;
    auto& bed = udata->bed;
    auto& fluid = udata->fluid;
    
    auto& C = state.concentration.var;

    auto& yi = state.mole_fraction.at(component).var;

    auto& molar_conc_eq = state.molar_concentration.at(component);
    auto& Ci = molar_conc_eq.var;
    
    int N = Ci.size();
    for (int i = 0; i < N; ++i) { // Interior cells
      double ci = C[i] * yi[i];
      Ci[i] = ci;
    }
  }

  void UpdateAdsorbedConcentration(UserData* udata, AdsorbentID ads_id, Component component) {
    auto& state = udata->state;
    const auto& bed = udata->bed;

    auto& qi = state.ads_concentration.at(ads_id).at(component).var;
    auto& ddt_qi = state.ads_concentration.at(ads_id).at(component).ddt_var;
    auto& res_qi = state.ads_concentration.at(ads_id).at(component).res_var;

    auto& qi_sat = state.eq_ads_concentration.at(ads_id).at(component).var;
    auto& Ki = state.mass_transfer_coefficient.at(ads_id).at(component).var;

    auto& domain = state.domain_vectors;
    auto& dx = domain.cell_width;
    const auto& is_void = domain.is_void_cell;
    const auto* mask = state.GetAdsorbentMask(ads_id);  // cells where this adsorbent is present

    // Kinetic rate model (LDF, LDF2, ...) for this adsorbent/component, resolved
    // once per call rather than per cell. AdsorptionSource::Initialise fills a
    // first-order LDF for every component, so a missing entry only happens for
    // an uninitialised source; fall back to first-order LDF in that case.
    const IMassTransferModel* rate_model = nullptr;
    {
      const auto& models = udata->adsorbent_directory.Get(ads_id).adsorption_source_.mass_transfer_models_;
      auto it = models.find(component);
      if (it != models.end()) rate_model = it->second.get();
    }

    // Manufactured rate source [mol/(kg s)] (verification hook).
    const auto* mms_loading = SourceOf(udata->manufactured_sources_.loading, component);

    int N = qi.size(); // Includes ghost cells
    for (int i = 1; i < N - 1; ++i) { // Interior cells
      if (i < (int)is_void.size() && is_void[i]) continue;  // void cells: qi is not part of the state vector
      if (mask && !(*mask)[i]) continue;  // adsorbent absent here: loading DOF is inactive
      const double rate = rate_model ? rate_model->Rate(qi_sat[i], qi[i], Ki[i])
                                     : Ki[i] * (qi_sat[i] - qi[i]);
      res_qi[i] = ddt_qi[i] * dx[i] - dx[i] * rate;
      if (mms_loading)
        res_qi[i] -= dx[i] * (*mms_loading)(CellCentre(domain, i), udata->residual_time_);
    }

    // Ghost cells mirror the first/last *bed* cell.
    const int first_bed = bed.GetFirstBedCellIndex();
    const int last_bed = bed.GetLastBedCellIndex();
    res_qi[0] = qi[0] - qi[first_bed];
    res_qi[N - 1] = qi[N - 1] - qi[last_bed];
  }

  void UpdateEquilibriumAdsorbedConcentration(UserData* udata, AdsorbentID ads_id, Component component) {
    auto& state = udata->state;
    auto& bed = udata->bed;
    auto& fluid = udata->fluid;
    auto& layers = bed.layers_.layers_;

    // Registered once per (adsorbent, component), but all of this adsorbent's
    // components are computed together on the first component's call (so
    // competitive/extended isotherms share their per-cell denominator); the
    // remaining per-component calls for this adsorbent are no-ops.
    if (component != fluid.components_.front())
      return;

    auto& adsorbent = udata->adsorbent_directory.Get(ads_id);
    auto& source = adsorbent.adsorption_source_;
    for (auto& layer : layers) {
      if (layer.composition_.find(ads_id) == layer.composition_.end()) continue;
      source.UpdateIsothermAll(state, ads_id, fluid.components_, layer.domain_.start_index, layer.domain_.end_index, &udata->variable_registry);
    }

    // Ghost cells mirror the first/last *bed* cell for this adsorbent's
    // components (skip voids if present).
    const int first_bed = bed.GetFirstBedCellIndex();
    const int last_bed = bed.GetLastBedCellIndex();
    for (auto& [comp, eq] : state.eq_ads_concentration.at(ads_id)) {
      auto& qi_sat = eq.var;
      const int N = (int)qi_sat.size();
      if (first_bed >= 0 && first_bed < N) qi_sat[0] = qi_sat[first_bed];
      if (last_bed >= 0 && last_bed < N) qi_sat[N - 1] = qi_sat[last_bed];
    }
  }

  void UpdateMoleFraction(UserData* udata, Component component) {
    auto& state = udata->state;
    auto& bed = udata->bed;
    auto& fluid = udata->fluid;
    auto& domain_vectors = state.domain_vectors;
    
    auto& mole_fraction = state.mole_fraction.at(component);
    auto& yi = mole_fraction.var;
    auto& ddt_yi = mole_fraction.ddt_var;
    auto& res_yi = mole_fraction.res_var;

    const auto& mole_fraction_face = state.mole_fraction_face.at(component);
    const auto& yi_face = mole_fraction_face.var;

    const auto& ads_conc = state.ads_concentration;  // nested: ads_id -> component -> equation
    const auto& cell_ads = domain_vectors.cell_adsorbents;

    const auto& pressure = state.pressure;
    const auto& P = pressure.var;
    const auto& ddt_P = pressure.ddt_var;

    const auto& temperature = state.fluid_temperature;
    const auto& T = temperature.var;
    const auto& ddt_T = temperature.ddt_var;

    auto& velocity = state.velocity;
    auto& u = velocity.var;

    auto& Dax = state.axial_dispersion_coefficient.var;

    const auto& dx = domain_vectors.cell_width;
    const auto& dV = domain_vectors.volume;
    const auto& A = domain_vectors.area;
    const auto& rho_b = domain_vectors.bulk_density;
    const auto& et = domain_vectors.et;
    const auto& eb = domain_vectors.eb;
    
    // Hoist the per-adsorbent lookup for THIS component out of the cell loop:
    // resolve &(dq/dt) once per adsorbent id (the map is the same object every
    // cell). A bed holds only a few adsorbents, so the per-cell lookup is a short
    // linear scan over a reused buffer (no hashing, no per-call allocation).
    static thread_local std::vector<std::pair<AdsorbentID, const std::vector<double>*>> ddt_q;
    ddt_q.clear();
    for (const auto& [ads_id, comp_map] : ads_conc) {
      auto it = comp_map.find(component);
      if (it != comp_map.end()) ddt_q.emplace_back(ads_id, &it->second.ddt_var);
    }
    auto ddt_q_of = [](AdsorbentID id) -> const std::vector<double>& {
      for (const auto& [ads_id, v] : ddt_q) if (ads_id == id) return *v;
      throw std::out_of_range("UpdateMoleFraction: unknown adsorbent");
    };

    // Manufactured component-balance source [mol/(m3 s)] (verification hook).
    const auto* mms_species = SourceOf(udata->manufactured_sources_.species, component);

    int N = yi.size();
    for (int n = 1; n < N - 1; ++n) { // Interior cells
      // Floor pressure/temperature denominators against a transient non-positive
      // state (no-op for physical P,T; floor << any real value). See constants.h.
      const double Pc = detail::DenomFloor(P[n], constants::kMinDenomPressure);
      const double Tc = detail::DenomFloor(T[n], constants::kMinDenomTemperature);

      double Pe = 0.5 * (P[n] + P[n + 1]);
      double Pw = 0.5 * (P[n] + P[n - 1]);
      double Te = 0.5 * (T[n] + T[n + 1]);
      double Tw = 0.5 * (T[n] + T[n - 1]);
      Te = detail::DenomFloor(Te, constants::kMinDenomTemperature);
      Tw = detail::DenomFloor(Tw, constants::kMinDenomTemperature);

      double u_e = u[n];
      double u_w = u[n - 1];

      double yie = yi_face[n];
      double yiw = yi_face[n - 1];

      double acc = dx[n] * et[n] * (-ddt_yi[n] - (yi[n] / Pc) * ddt_P[n] + (yi[n] / Tc) * ddt_T[n]);

      double fe = dx[n + 1] / (dx[n] + dx[n + 1]);
      double fw = dx[n] / (dx[n - 1] + dx[n]);

      double d_e = 1 / ((1 - fe) / Dax[n] + fe / Dax[n + 1]);
      double d_w = 1 / ((1 - fw) / Dax[n - 1] + fw / Dax[n]);

      // Danckwerts boundary faces: no dispersive flux crosses the two faces
      // of the bed (n == 1 west, n == N - 2 east), so the species flux through
      // a boundary face is the convective u C y_f alone, with y_f the boundary
      // value (the ghost cell the boundary condition sets). At an outlet or a
      // closed end the ghost mirrors the first cell and the term was zero
      // anyway; at an inlet the ghost holds the feed composition, and a
      // dispersive eb D C (y_1 - y_in) / dx would otherwise carry species into
      // the bed on top of the feed u C y_in, unseen by the stream logs (the
      // component balance then fails to close by an amount proportional to D).
      double D_e = (n == N - 2) ? 0.0
        : ((eb[n] * Tc * d_e) / (Pc * dx[n])) * (Pe / Te * (yi[n + 1] - yi[n]));
      double D_w = (n == 1) ? 0.0
        : ((eb[n] * Tc * d_w) / (Pc * dx[n])) * (Pw / Tw * (yi[n] - yi[n - 1]));

      double F_e = Tc / Pc * (yie * Pe * u_e / Te);
      double F_w = Tc / Pc * (yiw * Pw * u_w / Tw);

      // Adsorption source for this component, summed over the adsorbents present
      // in this cell, each weighted by its own bulk density: sum_k rho_b,k * dq_k/dt.
      double rho_b_ddt_q = 0.0;
      for (const auto& ca : cell_ads[n]) {
        rho_b_ddt_q += ca.bulk_density * ddt_q_of(ca.id)[n];
      }
      double adsorption_source = R * Tc * dx[n] * rho_b_ddt_q / Pc;

      res_yi[n] =
        acc                   // Accumulation
        + (D_e - D_w)         // Diffusion
        - (F_e - F_w)         // Convection
        - adsorption_source;  // Source

      // The residual is minus the balance times R T dx / P, as the
      // adsorption source above.
      if (mms_species)
        res_yi[n] += R * Tc * dx[n] * (*mms_species)(CellCentre(domain_vectors, n), udata->residual_time_) / Pc;
    }
  }

  void UpdateMoleFractionFaceUpwind(UserData* udata, Component component) {
    auto& state = udata->state;
    auto& yi_face = state.mole_fraction_face.at(component).var;
    const auto& yi = state.mole_fraction.at(component).var;
    const auto& u = state.velocity.var;

    int N = yi_face.size();
    for (int n = 0; n < N; ++n) { 
      double yi_f = UpwindFace(u[n], [&] { return yi[n]; }, [&] { return yi[n + 1]; });
      yi_face[n] = yi_f;
    }
  }

  namespace {
    // A cell whose value sits on a face rather than at a cell centre: the
    // ghost cells (Dirichlet inlet composition, or a copy of the end cell at
    // an outlet; either way the value at the bed-end face) and the lumped,
    // well-mixed void cells. A reconstruction may not interpolate through
    // such a cell as if it held a linear profile.
    inline bool IsLumpedCell(const DomainVectors& dv, int cell, int last_cell) {
      return cell <= 0 || cell >= last_cell || dv.IsVoidCell(static_cast<size_t>(cell));
    }

    // Bound a reconstructed interior face value y_f between the two cells
    // y_a, y_b that share the face: to [0, 1], widened to take in the two
    // cells' own values. A reconstruction can overshoot where a cell is a
    // local extremum (WENO's upwind extrapolation 3/2 y_c - 1/2 y_cu gives
    // -0.17 on the cells 0.45, 0, 1), and the bound stops that from carrying
    // a negative or > 1 composition. It must not be tighter than the cells
    // themselves, though: a component absent from the feed and the bed sits
    // at round-off around zero (and a pure one around one), and a plain
    // [0, 1] clamp zeroed the face of every cell whose value or stencil had
    // drifted marginally negative. Those cells could then not pass their
    // content on (the face, and its derivative, were zero), so negative
    // round-off was trapped while positive round-off left the bed, the
    // component drifted to -6e-10 (a few percent of its absolute tolerance),
    // and the face switching between clamped and free as Newton iterated
    // stalled the iteration: up to 25 times more steps than upwind, nearly
    // all after convergence failures (WP20-4). A face value between the two
    // cells it separates is never an overshoot, so allowing it keeps the
    // upwind behaviour of round-off (carried along and out) and leaves every
    // in-range profile exactly as [0, 1] did.
    inline double BoundFace(double y_f, double y_a, double y_b) {
      return std::clamp(y_f, std::min({0.0, y_a, y_b}), std::max({1.0, y_a, y_b}));
    }

    // Third-order WENO face value (Jiang and Shu 1996, k = 2; the scheme of
    // Haghpanah et al. 2013 and Yancy-Caballero et al. 2020, PEESEgroup/PSA
    // CycleSteps/WENO.m) at the face between cell c and its downwind
    // neighbour cd, reconstructed from the upwind side. cu is the cell
    // upwind of c. Two candidate face values, both linear in y:
    //   stencil {c, cd}: interpolation, y_c + (y_cd - y_c) dx_c / (dx_c + dx_cd)
    //   stencil {cu, c}: extrapolation, y_c + (y_c - y_cu) dx_c / (dx_cu + dx_c)
    // (uniform grid: (y_c + y_cd)/2 and 3/2 y_c - 1/2 y_cu), with linear
    // weights 2/3 and 1/3 and smoothness indicators beta = (slope * dx_c)^2
    // (uniform grid: (y_cd - y_c)^2 and (y_c - y_cu)^2). When cu is a lumped
    // cell its value sits on the face of c, half a cell away, so the
    // extrapolation becomes 2 y_c - y_cu with beta = (2 (y_c - y_cu))^2,
    // exactly the reference's first interior face.
    inline double WenoFace(double y_cu, double y_c, double y_cd,
                           double dx_cu, double dx_c, double dx_cd, bool cu_lumped) {
      constexpr double eps = 1e-12;
      const double h_up = cu_lumped ? 0.5 * dx_c : 0.5 * (dx_cu + dx_c);
      const double h_dn = 0.5 * (dx_c + dx_cd);
      const double s_up = (y_c - y_cu) / h_up;   // slope on the upwind stencil
      const double s_dn = (y_cd - y_c) / h_dn;   // slope on the central stencil
      const double cand_up = y_c + 0.5 * dx_c * s_up;
      const double cand_dn = y_c + 0.5 * dx_c * s_dn;
      const double beta_up = (s_up * dx_c) * (s_up * dx_c);
      const double beta_dn = (s_dn * dx_c) * (s_dn * dx_c);
      const double alpha_up = (1.0 / 3.0) / ((beta_up + eps) * (beta_up + eps));
      const double alpha_dn = (2.0 / 3.0) / ((beta_dn + eps) * (beta_dn + eps));
      return (alpha_up * cand_up + alpha_dn * cand_dn) / (alpha_up + alpha_dn);
    }
  }

  void UpdateMoleFractionFaceWENO(UserData* udata, Component component) {
    auto& state = udata->state;
    auto& yi_face = state.mole_fraction_face.at(component).var;
    const auto& yi = state.mole_fraction.at(component).var;
    const auto& u = state.velocity.var;
    const auto& dv = state.domain_vectors;
    const auto& dx = dv.cell_width;

    // Face n lies between cells n and n + 1; cells 0 and N are the ghosts.
    const int N = static_cast<int>(yi_face.size());
    for (int n = 0; n < N; ++n) {
      // Positive velocity: reconstruct from the left (cells n - 1, n, n + 1).
      // A reconstruction can overshoot where a cell is a local extremum and
      // produce a face mole fraction outside [0, 1]: the reconstructed value
      // is bounded to [0, 1], widened to the two cells' own values (see
      // BoundFace), so an overshoot cannot inject a negative or > 1
      // composition while round-off in the cells is still carried along as
      // under Upwind. Each face is single-valued, so every species stays
      // conserved. The upwind faces below are NOT bounded, as under Upwind:
      // the bed-end faces are not rescaled to the upwind cell's sum (see
      // RescaleInteriorMoleFractionFaces), and clamping an end cell's value
      // (a pure-N2 outlet at y = 1 + 1e-9) made the outlet face sum to one
      // instead of the cell's sum, so the cell's departure from one was no
      // longer carried out of the bed but grew with the outflow (WP20-1).
      auto left_biased = [&]() {
        // Upwind cell lumped (ghost: the inlet face value; void: well mixed,
        // its outflow has its composition) or downwind cell lumped (the
        // outlet face, or the inlet of a void): the upwind value, as the
        // reference's boundary faces. This keeps the boundary fluxes equal to
        // the stream integrals the cycle balance and KPIs are built from.
        if (IsLumpedCell(dv, n, N) || IsLumpedCell(dv, n + 1, N)) return yi[n];
        return BoundFace(WenoFace(yi[n - 1], yi[n], yi[n + 1], dx[n - 1], dx[n], dx[n + 1],
                                  IsLumpedCell(dv, n - 1, N)), yi[n], yi[n + 1]);
      };
      // Negative velocity: mirror image (cells n + 2, n + 1, n).
      auto right_biased = [&]() {
        if (IsLumpedCell(dv, n + 1, N) || IsLumpedCell(dv, n, N)) return yi[n + 1];
        return BoundFace(WenoFace(yi[n + 2], yi[n + 1], yi[n], dx[n + 2], dx[n + 1], dx[n],
                                  IsLumpedCell(dv, n + 2, N)), yi[n], yi[n + 1]);
      };
      yi_face[n] = UpwindFace(u[n], left_biased, right_biased);
    }
  }

  namespace {
    // Van Leer face value at the face between cell c and its downwind
    // neighbour cd, from the upwind side (cu upwind of c):
    // y_c + phi(r)/2 (y_cd - y_c), r = upwind slope / downwind slope,
    // phi(r) = (r + |r|)/(1 + |r|). r = a b / (b^2 + delta) is the ratio a/b
    // written so that it is the same for y and 1 - y (a two-component face
    // pair sums to one) and goes smoothly to 0 when b vanishes. A lumped cu
    // (ghost or void) sits on the face of c, half a cell away, which doubles
    // its slope.
    inline double VanLeerFace(double y_cu, double y_c, double y_cd, bool cu_lumped) {
      constexpr double delta = 1e-24;
      const double a = (cu_lumped ? 2.0 : 1.0) * (y_c - y_cu);
      const double b = y_cd - y_c;
      const double r = a * b / (b * b + delta);
      const double phi = (r + std::abs(r)) / (1.0 + std::abs(r));
      return y_c + 0.5 * phi * b;
    }

    // The convective face temperature of the energy balances, reconstructed
    // with the scheme the mole fractions use, so that T has the scheme's order
    // of accuracy too. Same stencils and the same faces left upwind: a face of
    // a lumped cell (a ghost or a void cell) takes the upwind value, and so
    // does every face under Upwind. Two differences from the mole fractions:
    //   * the temperature ghost cell does not hold the face value (it is set
    //     so that the mean of the ghost and the first cell is the boundary
    //     temperature), so the upstream value of a bed-end face, and the value
    //     a ghost contributes as the far cell of a reconstruction, is that
    //     mean (see BoundaryFaceUpstreamT);
    //   * the face value is not bounded: a temperature has no [0, 1] range,
    //     and the van Leer face lies between its two cells anyway;
    //   * the reconstruction acts on T / T_c, T_c the upwind cell's
    //     temperature, and the face value is scaled back. WENO's epsilon and
    //     van Leer's delta are absolute and sized for mole fractions of order
    //     one; on kelvin they put the noise-level differences of a nearly
    //     flat profile in the range where the weights switch, which made the
    //     residual noisy (more steps, run-to-run scatter). Relative to T_c a
    //     temperature has the scale of a mole fraction, and both
    //     reconstructions are otherwise homogeneous of degree one.
    // T_face[i] and T_face[i - 1] are the east and west faces of cell i, so
    // both cells of a face see the same value. The stencil reaches two cells
    // on either side, as the mole fractions' do, which the Jacobian band of
    // the high-order schemes already covers (BandReach in coupled_solver.cpp).
    void ConvectiveFaceTemperatures(const UserData* udata, const std::vector<double>& T,
                                    std::vector<double>& T_face) {
      const auto& u = udata->state.velocity.var;
      const auto& dv = udata->state.domain_vectors;
      const auto& dx = dv.cell_width;
      const AdvectionScheme scheme = udata->settings.advection_scheme_;
      const int last = static_cast<int>(T.size()) - 1;  // the second ghost; faces 0 .. last - 1
      T_face.resize(static_cast<size_t>(std::max(last, 0)));
      // Value of a lumped far cell at the face it shares with its bed neighbour.
      auto lumped_value = [&](int cell, int neighbour) {
        return (cell == 0 || cell == last) ? 0.5 * (T[cell] + T[neighbour]) : T[cell];
      };
      for (int n = 0; n < last; ++n) {
        // Positive velocity: from the left (cells n - 1, n, n + 1).
        auto left_biased = [&]() {
          const double upwind = (n == 0) ? 0.5 * (T[0] + T[1]) : T[n];
          if (scheme == AdvectionScheme::Upwind || IsLumpedCell(dv, n, last) || IsLumpedCell(dv, n + 1, last))
            return upwind;
          const bool cu_lumped = IsLumpedCell(dv, n - 1, last);
          const double T_cu = cu_lumped ? lumped_value(n - 1, n) : T[n - 1];
          const double s = T[n];
          return s * (scheme == AdvectionScheme::WENO
                          ? WenoFace(T_cu / s, 1.0, T[n + 1] / s, dx[n - 1], dx[n], dx[n + 1], cu_lumped)
                          : VanLeerFace(T_cu / s, 1.0, T[n + 1] / s, cu_lumped));
        };
        // Negative velocity: mirror image (cells n + 2, n + 1, n).
        auto right_biased = [&]() {
          const double upwind = (n + 1 == last) ? 0.5 * (T[n] + T[last]) : T[n + 1];
          if (scheme == AdvectionScheme::Upwind || IsLumpedCell(dv, n + 1, last) || IsLumpedCell(dv, n, last))
            return upwind;
          const bool cu_lumped = IsLumpedCell(dv, n + 2, last);
          const double T_cu = cu_lumped ? lumped_value(n + 2, n + 1) : T[n + 2];
          const double s = T[n + 1];
          return s * (scheme == AdvectionScheme::WENO
                          ? WenoFace(T_cu / s, 1.0, T[n] / s, dx[n + 2], dx[n + 1], dx[n], cu_lumped)
                          : VanLeerFace(T_cu / s, 1.0, T[n] / s, cu_lumped));
        };
        T_face[static_cast<size_t>(n)] = UpwindFace(u[static_cast<size_t>(n)], left_biased, right_biased);
      }
    }

    // The species fluxes add up to the total-mass flux only if the face mole
    // fractions sum to what the total balance carries. Upwind faces sum to the
    // upwind cell's sum, so a departure of the cell sum from one (which no
    // equation restores: every species and the total have their own balance)
    // is simply carried along and out of the bed. The nonlinear
    // reconstructions (WENO, van Leer), built component by component, do not
    // keep that sum, and worse, their linearisation about a steep front
    // weights the downwind cell, so a departure of the cell sum grows instead
    // of being carried: van Leer drove the cell sum to 90 in a CO2/N2
    // breakthrough. Rescaling each interior face so that its mole fractions
    // sum to the (smoothed) upwind cell's sum restores the upwind behaviour of
    // the sum while keeping the high-order profile of each component. Each
    // face stays single-valued, so every species stays conserved; the bed-end
    // faces are upwind already (and unclamped, so they keep the upwind cell's
    // sum exactly). Called after the last component's faces are
    // computed (the face equations run in the fluid's component order, see
    // Equations::Initialise); the coupling of each face to every component
    // over the stencil is why the high-order schemes widen the Jacobian band
    // (see BandReach in coupled_solver.cpp).
    void RescaleInteriorMoleFractionFaces(UserData* udata) {
      auto& state = udata->state;
      auto& faces = state.mole_fraction_face;
      if (faces.size() < 2) return;
      const auto& u = state.velocity.var;
      const int N = static_cast<int>(faces.begin()->second.var.size());
      for (int n = 1; n < N - 1; ++n) {
        double face_sum = 0.0, left_sum = 0.0, right_sum = 0.0;
        for (const auto& [comp, eq] : faces) {
          face_sum += eq.var[n];
          const auto& y = state.mole_fraction.at(comp).var;
          left_sum += y[n];
          right_sum += y[n + 1];
        }
        if (!(face_sum > 0.0)) continue;
        const double target = UpwindFace(u[n], [&] { return left_sum; }, [&] { return right_sum; });
        const double scale = target / face_sum;
        for (auto& [comp, eq] : faces) { (void)comp; eq.var[n] *= scale; }
      }
    }
  }

  void UpdateMoleFractionFaceVanLeer(UserData* udata, Component component) {
    auto& state = udata->state;
    auto& yi_face = state.mole_fraction_face.at(component).var;
    const auto& yi = state.mole_fraction.at(component).var;
    const auto& u = state.velocity.var;
    const auto& dv = state.domain_vectors;

    // Same stencil and boundary treatment as UpdateMoleFractionFaceWENO: the
    // bed-end faces and the faces of void cells take the upwind value,
    // unbounded. The limited face value lies between its two cells, so
    // BoundFace leaves it unchanged; a plain [0, 1] clamp acted only on
    // round-off in a cell marginally outside [0, 1], where it trapped that
    // round-off (a 4 to 15 times slower run with an absent component,
    // WP20-4).
    const int N = static_cast<int>(yi_face.size());
    for (int n = 0; n < N; ++n) {
      auto left_biased = [&]() {
        if (IsLumpedCell(dv, n, N) || IsLumpedCell(dv, n + 1, N)) return yi[n];
        return BoundFace(VanLeerFace(yi[n - 1], yi[n], yi[n + 1], IsLumpedCell(dv, n - 1, N)), yi[n], yi[n + 1]);
      };
      auto right_biased = [&]() {
        if (IsLumpedCell(dv, n + 1, N) || IsLumpedCell(dv, n, N)) return yi[n + 1];
        return BoundFace(VanLeerFace(yi[n + 2], yi[n + 1], yi[n], IsLumpedCell(dv, n + 2, N)), yi[n], yi[n + 1]);
      };
      yi_face[n] = UpwindFace(u[n], left_biased, right_biased);
    }
  }

  // Dispatch the convective face reconstruction to the scheme selected in the
  // solver settings. Upwind is the default; WENO and van Leer are higher-order
  // interior reconstructions. All three write the same face variable, so the
  // choice is a drop-in swap that leaves the conserved cell state untouched.
  // The high-order interior faces are rescaled to the upwind cell's sum of
  // mole fractions once the last component's are computed (see
  // RescaleInteriorMoleFractionFaces).
  void UpdateMoleFractionFace(UserData* udata, Component component) {
    switch (udata->settings.advection_scheme_) {
      case AdvectionScheme::WENO:
        UpdateMoleFractionFaceWENO(udata, component);
        break;
      case AdvectionScheme::VanLeer:
        UpdateMoleFractionFaceVanLeer(udata, component);
        break;
      case AdvectionScheme::Upwind:
      default:
        UpdateMoleFractionFaceUpwind(udata, component);
        return;
    }
    const auto& components = udata->fluid.components_;
    if (!components.empty() && component == components.back())
      RescaleInteriorMoleFractionFaces(udata);
  }

  void UpdateMassTransferCoefficient(UserData* udata, AdsorbentID ads_id, Component component){
    auto& state = udata->state;
    auto& bed = udata->bed;
    auto& fluid = udata->fluid;

    auto& adsorbent = udata->adsorbent_directory.Get(ads_id);
    auto& source = adsorbent.adsorption_source_;
    for (auto& layer : bed.layers_.layers_) {
      if (layer.composition_.find(ads_id) == layer.composition_.end()) continue;
      source.UpdateMassTransferCoefficient(state, ads_id, component, fluid, layer, &udata->variable_registry,
                                          adsorbent.solid_density_.Get());
    }

    ZeroVoidCellsAndMirrorGhosts(state.mass_transfer_coefficient.at(ads_id).at(component).var, bed);
  }

  static void BlanksLaw(UserData* udata, Component component) {
    auto& state = udata->state;
    auto& fluid = udata->fluid;

    const auto& P = state.pressure.var;
    const auto& T = state.fluid_temperature.var;
    auto& Dm = state.molecular_diffusion;
    
    auto& bed = udata->bed;
    // The pair table is already refreshed once by the caller
    // (UpdateMolecularDiffusionCoefficient) before this loop runs, so do NOT
    // recompute it here — that would repeat the expensive Neufield pow/exp pass
    // once per component per residual evaluation.
    auto& diffusion_coefficients = bed.equations_.binary_diffusion_model_.m_DiffusionCoefficients;

    const double eps = 1e-12;

    // Hoist the per-partner string-keyed lookups out of the cell loop: for this
    // component, resolve each partner's mole-fraction vector and binary-Dij
    // vector once (fluid.components_ order preserved → bit-identical sums).
    struct Partner { const std::vector<double>* yj; const std::vector<double>* Dij; };
    static thread_local std::vector<Partner> partners;
    partners.clear();
    partners.reserve(fluid.components_.size());
    for (const std::string& comp_j : fluid.components_) {
      if (comp_j == component) continue;
      partners.push_back({ &state.mole_fraction.at(comp_j).var,
                           &diffusion_coefficients[{component, comp_j}].Dm_ij });
    }
    auto& Dmi = Dm.at(component).var;

    for (int n = 0; n < P.size(); ++n) {
      double sum_yj = 0.0;
      for (const auto& p : partners) sum_yj += (*p.yj)[n];

      // Make positive
      double numerator = std::sqrt(sum_yj * sum_yj);

      double sum_yj_Dij = 0.0;
      for (const auto& p : partners) sum_yj_Dij += ((*p.yj)[n] / (*p.Dij)[n]);

      double denominator = std::sqrt(eps + std::pow(sum_yj_Dij, 2));

      Dmi[n] = numerator / denominator;
    }
  }

  static void MulticomponentDiffusionCorrection(UserData* udata, Component component) {
    auto& state = udata->state;
    auto& fluid = udata->fluid;

    const auto& P = state.pressure.var;
    const auto& T = state.fluid_temperature.var;
    auto& Dm = state.molecular_diffusion;
    
    auto& bed = udata->bed;
    bed.equations_.binary_diffusion_model_.Update(T, P);
    auto& diffusion_coefficients = bed.equations_.binary_diffusion_model_.m_DiffusionCoefficients;

    for (int n = 0; n < P.size(); ++n) {
      const std::string c1 = fluid.components_[0]; // primary reference
      const std::string c2 = fluid.components_[1]; // secondary reference

      for (const std::string& comp_i : fluid.components_) {
        std::string comp_j = c1;
        if (comp_i == c1) {
          comp_j = c2;
        }
        
        double D_ij = diffusion_coefficients[{comp_i, comp_j}].Dm_ij[n];

        double correction = 0.0;
        for (const std::string& comp_k : fluid.components_) {
          if (comp_k == comp_i || comp_k == comp_j) 
            continue;

          const auto& yk = state.mole_fraction.at(comp_k).var; 
          double D_ik = diffusion_coefficients[{comp_i, comp_k}].Dm_ij[n];
          correction += yk[n] * (D_ik - D_ij);
        }

        auto& Dmi = Dm.at(comp_i).var;
        Dmi[n] = D_ij + correction;
      }
    }
  }

  void UpdateMolecularDiffusionCoefficient(UserData* udata, Component component){
    auto& state = udata->state;
    auto& fluid = udata->fluid;
    auto& bed = udata->bed;

    size_t number_of_components = fluid.components_.size();

    const auto& P = state.pressure.var;
    const auto& T = state.fluid_temperature.var;
    auto& Dm = state.molecular_diffusion;
    
    if (number_of_components == 0)
      return;

    // This update is registered once per component, but a single evaluation
    // fills every component's coefficients (binary_diffusion_model_.Update
    // recomputes the whole pair table over all cells, with a pow + 3 exp per
    // cell). Do all the work on the first component's call and make the
    // remaining per-component calls no-ops — otherwise the expensive Neufield
    // omega recompute runs once per component per residual evaluation.
    if (component != fluid.components_.front())
      return;

    bed.equations_.binary_diffusion_model_.Update(T, P);
    auto& diffusion_coefficients = bed.equations_.binary_diffusion_model_.m_DiffusionCoefficients;

    if (number_of_components == 1) {
      // Pure fluid: D_m is the self-diffusion coefficient from the {i, i} pair
      // stored by BinaryDiffusionCoefficientModel::Initialise.
      const auto& component1 = fluid.components_[0];
      const auto& Dm_ii = diffusion_coefficients[{component1, component1}].Dm_ij;
      auto& Dm1 = Dm.at(component1).var;
      const int N = (int)P.size();
      for (int n = 0; n < N; ++n) {
        Dm1[n] = Dm_ii[n];
      }

      return;
    }

    if (number_of_components == 2) {
      const auto& component1 = fluid.components_[0];
      const auto& component2 = fluid.components_[1];
      // Hoist the string-pair map lookup and the destination vectors out of the
      // per-cell loop; the inner loop is then a plain copy.
      const auto& Dm_ij = diffusion_coefficients[{component1, component2}].Dm_ij;
      auto& Dm1 = Dm.at(component1).var;
      auto& Dm2 = Dm.at(component2).var;
      const int N = (int)P.size();
      for (int n = 0; n < N; ++n) {
        Dm1[n] = Dm_ij[n];
        Dm2[n] = Dm_ij[n];
      }

      return;
    }

    // Three or more components
    // MulticomponentDiffusionCorrection(udata, component);
    for (const auto& comp : fluid.components_) {
      BlanksLaw(udata, comp);
    }
  }

  void UpdateTotalMoleFraction(UserData* udata) {
    auto& state = udata->state;
    auto& y_total = state.total_mole_fraction.var;
    
    size_t N = y_total.size();
    for (size_t i = 0; i < N; ++i) {
      y_total[i] = 0.0;
    }
    
    for (const auto& [component, yi_eq] : state.mole_fraction) {
      const auto& yi = yi_eq.var;
      for (size_t i = 0; i < N; ++i) {
        y_total[i] += yi[i];
      }
    }
  }


  void UpdateFluidDensity(UserData* udata) {
    auto& state = udata->state;
    auto& fluid = udata->fluid;
    auto& rho = state.fluid_density.var;
    auto& molar_concentration = state.molar_concentration;

    // Resolve (C_i, M_i) once per component, in map order (bit-identical sum).
    static thread_local std::vector<std::pair<const std::vector<double>*, double>> terms;
    terms.clear();
    for (const auto& [component, Ci] : molar_concentration)
      terms.emplace_back(&Ci.var, fluid.molar_mass_.at(component));

    for (int n = 0; n < rho.size(); ++n) {
      rho[n] = 0;
      for (const auto& [ci, Mi] : terms) {
        rho[n] += (*ci)[n] * Mi;
      }
    }
  }
  
  void UpdateViscosity(UserData* udata) {
    auto& state = udata->state;
    auto& viscosity = state.viscosity.var;
    auto& temperature = state.fluid_temperature.var;
    auto& mole_fraction_map = state.mole_fraction;

    // CONSTANT viscosity mode: the user's mixture value everywhere.
    const double mu_const = udata->fluid.ConstantViscosity();
    if (mu_const > 0.0) {
      std::fill(viscosity.begin(), viscosity.end(), mu_const);
      return;
    }

    // Mole-fraction weighted mixing rule mu = sum_i y_i * mu_i(T). Resolve each
    // component's mole-fraction vector and its hard-coded mu_i(T) coefficients
    // ONCE (out of the cell loop) so the per-cell inner loop is pure arithmetic —
    // no CoolProp calls and no string/map lookups on the solver hot path.
    static thread_local std::vector<std::pair<const std::vector<double>*, const std::array<double, 5>*>> comps;
    comps.clear();
    for (auto& [component, mole_fraction] : mole_fraction_map)
      comps.emplace_back(&mole_fraction.var, Fluid::ViscosityCoeffs(component));

    for (int n = 0; n < viscosity.size(); ++n) {
      const double T = temperature[n];
      double mu = 0.0;
      for (const auto& [yi, coeffs] : comps)
        mu += (*yi)[n] * (coeffs ? Fluid::EvalViscosity(*coeffs, T) : 1.0e-5);
      viscosity[n] = mu;
    }
  }

  void UpdateHeatCapacity(UserData* udata) {
    auto& state = udata->state;
    auto& Cp = state.heat_capacity.var;
    auto& mole_fraction_map = state.mole_fraction;
    const auto& fluid = udata->fluid;

    // Mass-specific mixture heat capacity [J/(kg.K)], mixed on the molar basis:
    // Cp = sum_i y_i M_i Cp_i / sum_i y_i M_i (Fluid::MixtureHeatCapacityFromSums).
    // Hoist the string-keyed lookups out of the cell loop: resolve
    // (mole-fraction vector, M_i Cp_i, M_i) once per component (map iteration
    // order preserved → deterministic sums), then the cell loop is pure FMA.
    struct CpTerm { const std::vector<double>* yi; double molar_cp; double molar_mass; };
    static thread_local std::vector<CpTerm> terms;
    terms.clear();
    terms.reserve(mole_fraction_map.size());
    for (auto& [component, mole_fraction] : mole_fraction_map) {
      terms.push_back({ &mole_fraction.var, fluid.MolarHeatCapacity(component), fluid.molar_mass_.at(component) });
    }

    // Compute mixture heat capacity for each cell
    for (int n = 1; n < Cp.size() - 1; ++n) {
      double sum_cp = 0.0, sum_m = 0.0;
      for (const auto& t : terms) {
        const double y = (*t.yi)[n];
        sum_cp += y * t.molar_cp;
        sum_m += y * t.molar_mass;
      }
      Cp[n] = Fluid::MixtureHeatCapacityFromSums(sum_cp, sum_m);
    }

    Cp[0] = Cp[1];
    Cp[Cp.size() - 1] = Cp[Cp.size() - 2];
  }

  void UpdateSherwoodNumber(UserData* udata) {
    auto& state = udata->state;
    auto& Sh = state.sherwood_number.var;
    auto& Sc = state.schmidt_number.var;
    auto& Re = state.reynolds_number.var;
    const auto& dv = state.domain_vectors;

    size_t N = Sh.size();
    for (size_t n = 0; n < N; ++n) {
      if (n < dv.channel_nu.size() && dv.channel_nu[n] > 0.0) {
        // A structured channel: the fully developed laminar Sh (= Nu_T by the
        // heat-mass analogy; Shah and London 1978), or with entrance effects
        // Hawthorn's (1974) length-averaged developing-flow form
        //   Sh = Sh_inf (1 + 0.095 Re Sc D_h / L)^0.45.
        Sh[n] = ChannelFilmNumber(dv.channel_nu[n], std::max(Re[n], 0.0), std::max(Sc[n], 0.0),
                                  dv.hydraulic_diameter[n], dv.entrance_length[n]);
        continue;
      }
      // Packed bed: Wakao and Funazkri (1978), Sh = 2 + 1.1 Sc^(1/3) Re^0.6.
      Sh[n] = 2 + 1.1 * std::pow(Sc[n], 0.3333333) * std::pow(Re[n], 0.6);
    }
  }

  void UpdateReynoldsNumber(UserData* udata) {
    auto& state = udata->state;
    auto& domain = state.domain_vectors;
    // dp for a packed bed (particle Re on the superficial velocity); D_h / OFA
    // for a structured layer (channel Re on the channel velocity u / OFA).
    auto& dp = domain.re_length;

    auto& Re = state.reynolds_number.var;
    auto& rho = state.fluid_density.var;
    auto& u = state.velocity.var;
    auto& vis = state.viscosity.var;

    size_t N = Re.size();
    for (size_t n = 1; n < N - 1; ++n) {
      Re[n] = rho[n] * std::abs(u[n]) * dp[n] / vis[n];
    }

    Re[0] = Re[1];
    Re[N - 1] = Re[N - 2];
  }

  void UpdateSchmidtNumber(UserData* udata) {
    auto& state = udata->state;
    auto& Sc = state.schmidt_number.var;
    auto& vis = state.viscosity.var;
    auto& rho = state.fluid_density.var;
    auto& Dm = state.molecular_diffusion;
    auto& mole_frac = state.mole_fraction;
    
    constexpr double D_min = 1e-9;   // [m^2/s] choose something small but > 0

    // Resolve (y_i, Dm_i) once per component, in map order (bit-identical sum).
    static thread_local std::vector<std::pair<const std::vector<double>*, const std::vector<double>*>> terms;
    terms.clear();
    for (const auto& [component, yi] : mole_frac)
      terms.emplace_back(&yi.var, &Dm.at(component).var);

    size_t N = Sc.size();
    for (size_t n = 1; n < N - 1; ++n) {
      double D = 0;
      for (const auto& [yi, Dmi] : terms) {
        D += (*yi)[n] * (*Dmi)[n];
      }
      D = std::max(D, D_min);
      Sc[n] = vis[n] / (rho[n] * D);
    }

    Sc[0] = Sc[1];
    Sc[N - 1] = Sc[N - 2];
  }

  void UpdateFluidThermalConductivity(UserData* udata) {
    auto& state = udata->state;
    auto& kg = state.fluid_thermal_conductivity.var;
    auto& temperature = state.fluid_temperature.var;
    auto& mole_fraction_map = state.mole_fraction;

    // Mole-fraction weighted mixing rule k = sum_i y_i * k_i(T), using the hard-
    // coded per-component conductivity polynomials. Coefficient pointers are
    // resolved ONCE (out of the cell loop) so the per-cell inner loop is pure
    // arithmetic — no CoolProp calls / map lookups on the solver hot path.
    static thread_local std::vector<std::pair<const std::vector<double>*, const std::array<double, 5>*>> comps;
    comps.clear();
    for (auto& [component, mole_fraction] : mole_fraction_map)
      comps.emplace_back(&mole_fraction.var, Fluid::ThermalConductivityCoeffs(component));

    for (int n = 0; n < kg.size(); ++n) {
      const double T = temperature[n];
      double k = 0.0;
      for (const auto& [yi, coeffs] : comps)
        k += (*yi)[n] * (coeffs ? Fluid::EvalThermalConductivity(*coeffs, T) : 0.024);
      kg[n] = k;
    }
  }

  void UpdateInternalHeatTransferCoefficient(UserData* udata) {
    auto& state = udata->state;
    auto& bed = udata->bed;

    for (auto& layer : bed.layers_.layers_) {
      auto& eq = layer.equations_.internal_heat_transfer_model_;
      eq->Update(state, layer, bed.wall_);
    }

    ZeroVoidCellsAndMirrorGhosts(state.internal_heat_transfer_coefficient.var, bed);
  }

  void UpdateExternalHeatTransferCoefficient(UserData* udata) {
    auto& state = udata->state;
    auto& bed = udata->bed;

    // A step with a heating medium replaces every layer's outer-face
    // coefficient with the medium's h_m (one jacket around the whole bed).
    const HeatingMedium* medium = ActiveHeatingMedium(state);
    auto& ho = state.external_heat_transfer_coefficient.var;
    for (auto& layer : bed.layers_.layers_) {
      if (medium) {
        const double h_m = medium->heat_transfer_coefficient_;
        for (int n = layer.domain_.start_index; n <= layer.domain_.end_index; ++n) ho[n] = h_m;
      } else {
        auto& eq = layer.equations_.external_heat_transfer_model_;
        eq->Update(state, layer, bed.wall_);
      }
    }

    ZeroVoidCellsAndMirrorGhosts(state.external_heat_transfer_coefficient.var, bed);
  }

  void UpdateAxialThermalConductivity(UserData* udata) {
    auto& state = udata->state;
    auto& bed = udata->bed;

    for (auto& layer : bed.layers_.layers_) {
      auto& eq = layer.equations_.heat_dispersion_model_;
      eq->Update(state, layer);
    }

    ZeroVoidCellsAndMirrorGhosts(state.axial_thermal_conductivity.var, bed);
  }

  void UpdateAxialDispersionCoefficient(UserData* udata) {
    auto& state = udata->state;
    auto& bed = udata->bed;

    for (auto& layer : bed.layers_.layers_) {
      auto& eq = layer.equations_.dispersion_model_;
      eq->Update(state, layer);
    }

    ZeroVoidCellsAndMirrorGhosts(state.axial_dispersion_coefficient.var, bed);
  }

  void UpdateMolarMass(UserData* udata) {
    auto& state = udata->state;
    auto& fluid = udata->fluid;
    auto& Mr = state.molar_mass.var;
    auto& mole_fraction_map = state.mole_fraction;

    // Resolve (y_i, M_i) once per component, in map order (bit-identical sum).
    static thread_local std::vector<std::pair<const std::vector<double>*, double>> terms;
    terms.clear();
    for (auto& [component, mole_fraction] : mole_fraction_map)
      terms.emplace_back(&mole_fraction.var, fluid.molar_mass_.at(component));

    for (int n = 0; n < Mr.size(); ++n) {
      Mr[n] = 0;
      for (const auto& [yi, Mi] : terms) {
        Mr[n] += (*yi)[n] * Mi;
      }
    }
  }

  void UpdateVolumetricFlowRate(UserData* udata) {
    auto& state = udata->state;
    auto& bed = udata->bed;
    auto& domain_vectors = state.domain_vectors;

    auto& volumetric_flow_rate = state.volumetric_flow_rate.var;
    const auto& velocity = state.velocity.var;
    const auto& A = domain_vectors.area;

    int N = volumetric_flow_rate.size();
    for (int n = 0; n < N; ++n) {
      volumetric_flow_rate[n] = velocity[n] * A[n];
    }
  }

  void UpdateMolarFlowRate(UserData* udata) {
    auto& state = udata->state;
    auto& domain_vectors = state.domain_vectors;

    auto& molar_flow_rate = state.molar_flow_rate.var;
    const auto& volumetric_flow_rate = state.volumetric_flow_rate.var;
    const auto& pressure = state.pressure.var;
    const auto& temperature = state.fluid_temperature.var;

    const double R = constants::kIdealGasConstant;
    int N = molar_flow_rate.size();
    for (int n = 0; n < N; ++n) {
      // Pressure and temperature at face (average of adjacent cells)
      double Pf = 0.5 * (pressure[n] + pressure[n + 1]);
      double Tf = 0.5 * (temperature[n] + temperature[n + 1]);
      double Cf = Pf / (R * Tf); // Face concentration
      molar_flow_rate[n] = volumetric_flow_rate[n] * Cf;
    }
  }

  void UpdateMassSource(UserData* udata) {
    auto& state = udata->state;
    auto& mass_source = state.mass_source.var;
    const auto& component_mass_source = state.component_mass_source;

    int N = mass_source.size();
    for (int n = 0; n < N; ++n) {
      mass_source[n] = 0.0;
      for (const auto& [component, Smi] : component_mass_source) {
        mass_source[n] += Smi.var[n];
      }
    }
  }

  void UpdateEnergySource(UserData* udata) {
    auto& state = udata->state;
    auto& energy_source = state.energy_source.var;
    // Isothermal run: electric sources and heat of adsorption are inactive.
    if (state.IsIsothermal()) { std::fill(energy_source.begin(), energy_source.end(), 0.0); return; }
    const auto& component_energy_source = state.component_energy_source;
    const auto& Qe = state.domain_vectors.energy_source_power;

    int N = energy_source.size();
    for (int n = 0; n < N; ++n) {
      energy_source[n] = Qe[n];
      for (const auto& [component, Sei] : component_energy_source) {
        energy_source[n] += Sei.var[n];
      }
    }
  }

  void UpdateApproachToFluidisation(UserData* udata) {
    auto& state = udata->state;
    auto& bed = udata->bed;
    auto& domain_vectors = state.domain_vectors;
    auto& approach_to_fluidisation = state.approach_to_fluidisation.var;
    const auto& velocity = state.velocity.var;
    const auto& dp = domain_vectors.dp;
    const auto& rho = state.fluid_density.var;
    const auto& vis = state.viscosity.var;
    const auto& dPdx = state.pressure_gradient.var;

    const double g = constants::kGravity; 
    const auto& rho_b = domain_vectors.bulk_density;
    const auto& epsilon = domain_vectors.eb;
    const auto& dx = domain_vectors.cell_width;
    
    int N = approach_to_fluidisation.size();
    for (int n = 0; n < N; ++n) {
      double dPdx_max = -(rho_b[n] - rho[n]) * g;
      approach_to_fluidisation[n] = dPdx[n] / dPdx_max;    
    }
  }

  void UpdatePressureGradient(UserData* udata) {
    auto& state = udata->state;
    auto& domain_vectors = state.domain_vectors;
    auto& pressure_gradient = state.pressure_gradient.var;
    const auto& pressure = state.pressure.var;
    const auto& dx = domain_vectors.cell_width;

    int N = pressure_gradient.size();
    for (int n = 1; n < N - 1; ++n) {
      // Central difference for pressure gradient
      pressure_gradient[n] = (pressure[n + 1] - pressure[n]) / (dx[n]);
    }
    
    // Boundary conditions: extrapolate from interior
    if (N > 1) {
      pressure_gradient[0] = pressure_gradient[1];
      pressure_gradient[N - 1] = pressure_gradient[N - 2];
    }
  }

  void UpdateComponentMassSource(UserData* udata, Component component) {
    auto& state = udata->state;
    auto& Smi = state.component_mass_source.at(component).var;
    const auto& ads_conc = state.ads_concentration;
    const auto& cell_ads = state.domain_vectors.cell_adsorbents;

    // Resolve this component's dq/dt per adsorbent once (a bed holds only a
    // few adsorbents, so the per-cell lookup is a short linear scan).
    static thread_local std::vector<std::pair<AdsorbentID, const std::vector<double>*>> ddt_q;
    ddt_q.clear();
    for (const auto& [ads_id, comp_map] : ads_conc)
      ddt_q.emplace_back(ads_id, &comp_map.at(component).ddt_var);
    auto ddt_q_of = [](AdsorbentID id) -> const std::vector<double>& {
      for (const auto& [ads_id, v] : ddt_q) if (ads_id == id) return *v;
      throw std::out_of_range("UpdateComponentMassSource: unknown adsorbent");
    };

    // Component mass source from adsorption, summed per adsorbent: sum_k rho_b,k * dq_{k,comp}/dt
    int N = Smi.size();
    for (int n = 0; n < N; ++n) {
      double s = 0.0;
      for (const auto& ca : cell_ads[n]) {
        s += ca.bulk_density * ddt_q_of(ca.id)[n];
      }
      Smi[n] = s;
    }
  }

  void UpdateComponentEnergySource(UserData* udata, Component component) {
    auto& state = udata->state;
    auto& Sei = state.component_energy_source.at(component).var;
    // Isothermal run: heat of adsorption is inactive (no energy balance).
    if (state.IsIsothermal()) { std::fill(Sei.begin(), Sei.end(), 0.0); return; }
    const auto& ads_conc = state.ads_concentration;
    const auto& cell_ads = state.domain_vectors.cell_adsorbents;
    const auto& dH_ads = state.domain_vectors.dH_ads;

    // Use the same heat per mole adsorbed as the energy residuals
    // (AdsorptionHeat), so the logged source is the term the residual actually
    // uses: at the fluid temperature in LTE (UpdateTemperature) and at the
    // solid temperature in LTNE (UpdateSolidTemperature).
    const auto& fluid = udata->fluid;
    const HeatOfAdsorptionConvention hoa_convention = fluid.heat_of_adsorption_convention_;
    const double Cpa = hoa_convention == HeatOfAdsorptionConvention::INTERNAL_ENERGY
                           ? fluid.AdsorbedPhaseHeatCapacity(component) : 0.0;
    const auto& T_ads =
      (udata->bed.equations_.thermal_equilibrium_mode_ == ThermalEquilibriumMode::NON_LOCAL)
        ? state.solid_temperature.var
        : state.fluid_temperature.var;

    // Resolve this component's (dH, dq/dt) per adsorbent once (a bed holds only
    // a few adsorbents, so the per-cell lookup is a short linear scan).
    struct Term { AdsorbentID id; const std::vector<double>* dH; const std::vector<double>* ddt_q; };
    static thread_local std::vector<Term> terms;
    terms.clear();
    for (const auto& [ads_id, comp_map] : ads_conc)
      terms.push_back({ ads_id, &dH_ads.at(ads_id).at(component), &comp_map.at(component).ddt_var });
    auto term_of = [](AdsorbentID id) -> const Term& {
      for (const auto& t : terms) if (t.id == id) return t;
      throw std::out_of_range("UpdateComponentEnergySource: unknown adsorbent");
    };

    // Component energy source from adsorption, summed per adsorbent:
    // sum_k heat_{k,comp} * rho_b,k * dq_{k,comp}/dt
    int N = Sei.size();
    for (int n = 0; n < N; ++n) {
      double s = 0.0;
      for (const auto& ca : cell_ads[n]) {
        const Term& t = term_of(ca.id);
        const double heat = AdsorptionHeat((*t.dH)[n], Cpa, T_ads[n], hoa_convention);
        s += heat * ca.bulk_density * (*t.ddt_q)[n];
      }
      Sei[n] = s;
    }
  }

} // end adsim namespace
