/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Method of manufactured solutions (MMS): order of accuracy of the whole
// discretisation, with every term of the isothermal model active at once.
//
// The other order tests switch terms off to reach a closed form (an inert
// tracer for Ogata-Banks, no dispersion for a travelling step, a uniform bed
// for the energy ODEs). Here the solution is CHOSEN instead: smooth, time-
// dependent fields P(z,t), y(z,t), u(z,t), q(z,t), substituted into the
// continuous equations of the documented model (docs/mass-balance.md,
// docs/momentum-balance.md, docs/kinetic-models.md); what does not vanish is
// a source, which the test hands to the solver through the verification hook
// UserData::manufactured_sources_ (adsim/equations/manufactured_sources.h).
// The chosen fields are then the exact solution of the continuous problem,
// and the difference between them and the discrete solution is the
// discretisation error alone. A wrong coefficient anywhere (a missing eb, a
// T/P scaling, a face average in place of an upwind value) shows up as an
// error that does not converge -- the class of defect the closed-form tests,
// each with most terms switched off, cannot see.
//
// The model (isothermal, ideal gas C = P/(R T), superficial velocity u):
//   overall:   et dC/dt + d(u C)/dz + rho_b dq_A/dt                       = S_tot
//   species i: et d(C y_i)/dt + d(u C y_i)/dz - d(eb D C dy_i/dz)/dz
//              + rho_b dq_i/dt                                             = S_i
//   loading:   dq_A/dt = k (q_A*(P y_A) - q_A) + S_q,  q* single-site Langmuir
//   momentum:  -dP/dz = (mu / K) u - S_m,  K = eb^3 dp^2 / (180 (1 - eb)^2)
// with et = eb + (1 - eb) ep and rho_b = (1 - eb) rho_p. Component B is inert.
// The sources are derived below BY HAND from these equations (the product
// rule on C = P/(R T)); the test also re-derives every one by central finite
// differences of the manufactured fields and checks the two agree, so a slip
// in the hand derivation cannot hide behind the solver.
//
// The manufactured solution (zeta = z / L, g(t) = 1 - exp(-t / tau)):
//   P = P_L + g [dP (1 - zeta) + a sin(pi zeta)]
//   y = y0 + b g (1 - cos(pi zeta)) / 2         (y_A; y_B = 1 - y)
//   u = u_in [1 + c g sin(pi zeta / 2)]
//   q = q0 + d g cos(pi zeta / 2),  q0 = q*(P_L y0)
// It is chosen to meet the boundary and initial conditions the deck can
// express, which are constants:
//   * t = 0: uniform P_L, y0 and q0 = q*(P_L y0), the deck's initial state
//     (the loading starts at equilibrium with the initial gas);
//   * z = 0, a velocity inlet: u(0,t) = u_in, and y(0,t) = y0 with
//     dy/dz(0,t) = 0, so the Danckwerts inlet u C y_in = u C y - eb D C dy/dz
//     holds with the constant feed composition y_in = y0;
//   * z = L, a pressure outlet: P(L,t) = P_L, and dy/dz(L,t) = 0, the zero-
//     gradient (no dispersive flux) outlet the bed-end face imposes.
// No cosine or sine has a zero or an extremum inside the bed except where
// stated, the fields stay in their physical ranges (0 < y < 1, P, q, u > 0),
// and u > 0 everywhere so every face has a fixed upwind side.
//
// Expected orders (the assertions): the dispersion, pressure and Darcy terms
// are central differences, the face pressures and temperatures arithmetic
// means, and the sources and accumulations are sampled at the cell centre,
// all O(dz^2). The convective face value (of the mole fractions and of the
// temperature) is upwind (O(dz)), WENO (third-order reconstruction) or van
// Leer (second order, limited). So the error
// norm converges at order 1 for upwind and 2 for WENO and van Leer, which the
// O(dz^2) central terms cap. The bed-end faces take the upwind value in every
// scheme, which is first order in general; this solution has dy/dz = 0 at
// both ends (the Danckwerts conditions require it), so there they are second
// order too and the high-order schemes keep their order up to the bed ends.
// The tolerances on the orders come from the measured asymptotic range (see
// the numbers in the README).

#include "validation_helpers.h"

#include <cmath>
#include <functional>
#include <string>
#include <vector>

namespace adsim {

  using namespace adsim::validation;

  namespace {

    constexpr double kPi = 3.14159265358979323846;
    constexpr double kRgas = kR;  // [J/(mol K)], the CODATA value (validation_helpers.h)

    // ---- Parameters of the case and the manufactured solution ---------------
    struct MmsCase {
      // Bed and gas.
      double L = 0.5;            // [m]
      double diameter = 0.05;    // [m]
      double eb = 0.4;           // bed void fraction [-]
      double ep = 0.3;           // particle void fraction [-]
      double rho_p = 1000.0;     // particle density [kg/m3]
      double dp = 3.0e-4;        // particle diameter [m]
      double mu = 1.8e-5;        // constant gas viscosity [Pa s]
      double D = 5.0e-3;         // constant axial dispersion coefficient [m2/s]
      double k = 0.5;            // constant LDF coefficient [1/s]
      double T = 298.15;         // [K]
      // Single-site Langmuir of A: q* = m1 K p / (1 + K p), K = b1 exp(q1 / R T).
      double m1 = 3.0;           // [mol/kg]
      double b1 = 2.0e-8;        // [1/Pa]
      double q1 = 20000.0;       // [J/mol]
      // Boundary data.
      double P_L = 1.0e5;        // outlet pressure [Pa]
      double u_in = 0.05;        // inlet superficial velocity [m/s]
      double y0 = 0.3;           // feed and initial mole fraction of A [-]
      // Manufactured-solution amplitudes.
      double tau = 3.0;          // [s]
      double dP = 5000.0;        // [Pa], close to the Darcy drop of u_in (5.06 kPa)
      double a = 1000.0;         // [Pa]
      double b = 0.2;            // [-]
      double c = 0.3;            // [-]
      double d = 0.005;          // [mol/kg]
      double t_end = 6.0;        // [s], two time constants: well inside the transient

      // Non-isothermal case (energy balance in LTE, adiabatic: no wall
      // exchange), off by default.
      bool thermal = false;
      double e = 20.0;           // temperature amplitude [K]
      double Cp_A = 29.4;        // constant molar gas heat capacities [J/(mol K)]
      double Cp_B = 29.1;
      double Cps = 900.0;        // solid heat capacity [J/(kg K)]
      double Cpa = 40.0;         // adsorbed-phase heat capacity of A [J/(mol K)]
      double dH = 20000.0;       // isosteric heat of adsorption of A [J/mol]
      double lambda = 5.0;       // constant axial thermal conductivity [W/(m K)]

      double et() const { return eb + (1.0 - eb) * ep; }
      double rho_b() const { return (1.0 - eb) * rho_p; }
      double Kperm() const { return eb * eb * eb * dp * dp / (180.0 * (1.0 - eb) * (1.0 - eb)); }
      double LangmuirK(double TT) const { return b1 * std::exp(q1 / (kRgas * TT)); }
      double qstar(double pA, double TT) const {
        const double K = LangmuirK(TT);
        return m1 * K * pA / (1.0 + K * pA);
      }
      double q0() const { return qstar(P_L * y0, T); }
    };

    // Manufactured fields and the derivatives the sources need.
    struct Fields {
      double P, Pz, Pt;
      double y, yz, yzz, yt;
      double u, uz;
      double q, qt;
      double T, Tz, Tzz, Tt;
    };

    Fields Manufactured(const MmsCase& m, double z, double t) {
      const double zeta = z / m.L;
      const double g = 1.0 - std::exp(-t / m.tau);
      const double gt = std::exp(-t / m.tau) / m.tau;
      const double s1 = std::sin(kPi * zeta), c1 = std::cos(kPi * zeta);
      const double sh = std::sin(0.5 * kPi * zeta), ch = std::cos(0.5 * kPi * zeta);
      Fields f;
      const double Pshape = m.dP * (1.0 - zeta) + m.a * s1;
      f.P = m.P_L + g * Pshape;
      f.Pz = g * (-m.dP / m.L + m.a * kPi / m.L * c1);
      f.Pt = gt * Pshape;
      f.y = m.y0 + 0.5 * m.b * g * (1.0 - c1);
      f.yz = 0.5 * m.b * g * kPi / m.L * s1;
      f.yzz = 0.5 * m.b * g * (kPi / m.L) * (kPi / m.L) * c1;
      f.yt = 0.5 * m.b * gt * (1.0 - c1);
      f.u = m.u_in * (1.0 + m.c * g * sh);
      f.uz = m.u_in * m.c * g * 0.5 * kPi / m.L * ch;
      f.q = m.q0() + m.d * g * ch;
      f.qt = m.d * gt * ch;
      // T = T0 + e g (1 - cos(pi zeta)) / 2: the inlet temperature at z = 0 and
      // zero gradient at both ends (no conduction crosses a bed-end face).
      const double e = m.thermal ? m.e : 0.0;
      f.T = m.T + 0.5 * e * g * (1.0 - c1);
      f.Tz = 0.5 * e * g * kPi / m.L * s1;
      f.Tzz = 0.5 * e * g * (kPi / m.L) * (kPi / m.L) * c1;
      f.Tt = 0.5 * e * gt * (1.0 - c1);
      return f;
    }

    // The terms of each balance at (z, t), from the product rule on
    // C = P / (R T): C_t = C (P_t / P - T_t / T), C_z = C (P_z / P - T_z / T),
    // (u C)_z = u_z C + u C_z, (C y)_t = C_t y + C y_t,
    // (u C y)_z = u_z C y + u C_z y + u C y_z, (C y_z)_z = C_z y_z + C y_zz.
    // Kept as separate terms so the guards can check that each is large.
    struct Terms {
      // Overall balance.
      double tot_acc, tot_conv, tot_ads;
      // Species A (B is the same with y -> 1 - y and no adsorption).
      double A_acc, A_conv, A_disp, A_ads;
      double B_acc, B_conv, B_disp;
      // Loading A.
      double q_acc, q_ldf;
      // Momentum.
      double m_darcy, m_grad;
      // Energy (LTE, adiabatic): Cbed T_t + C Cp u T_z - lambda T_zz - et P_t
      // - dH rho_b q_t = S_T, Cbed = et C Cp + rho_b (Cps + Cpa q).
      double E_acc, E_conv, E_cond, E_work, E_ads;
      double Cbed;  // [J/(m3 K)]
      double Kp;  // K p_A: the Langmuir nonlinearity
    };

    Terms TermsAt(const MmsCase& m, double z, double t) {
      const Fields f = Manufactured(m, z, t);
      const double C = f.P / (kRgas * f.T);
      const double Ct = C * (f.Pt / f.P - f.Tt / f.T);
      const double Cz = C * (f.Pz / f.P - f.Tz / f.T);
      const double et = m.et(), eb = m.eb, D = m.D, rho_b = m.rho_b();
      Terms s;
      s.tot_acc = et * Ct;
      s.tot_conv = f.uz * C + f.u * Cz;
      s.tot_ads = rho_b * f.qt;
      s.A_acc = et * (Ct * f.y + C * f.yt);
      s.A_conv = f.uz * C * f.y + f.u * Cz * f.y + f.u * C * f.yz;
      s.A_disp = -eb * D * (Cz * f.yz + C * f.yzz);
      s.A_ads = rho_b * f.qt;
      const double yB = 1.0 - f.y;
      s.B_acc = et * (Ct * yB - C * f.yt);
      s.B_conv = f.uz * C * yB + f.u * Cz * yB - f.u * C * f.yz;
      s.B_disp = +eb * D * (Cz * f.yz + C * f.yzz);
      s.q_acc = f.qt;
      s.q_ldf = -m.k * (m.qstar(f.P * f.y, f.T) - f.q);
      s.m_darcy = m.mu / m.Kperm() * f.u;
      s.m_grad = f.Pz;
      const double Cp = f.y * m.Cp_A + yB * m.Cp_B;  // molar mixture heat capacity
      const double Cbed = et * C * Cp + rho_b * (m.Cps + m.Cpa * f.q);
      s.Cbed = Cbed;
      s.E_acc = Cbed * f.Tt;
      s.E_conv = C * Cp * f.u * f.Tz;
      s.E_cond = -m.lambda * f.Tzz;
      s.E_work = -et * f.Pt;
      s.E_ads = -m.dH * rho_b * f.qt;
      s.Kp = m.LangmuirK(f.T) * f.P * f.y;
      return s;
    }

    double SourceTotal(const MmsCase& m, double z, double t) {
      const Terms s = TermsAt(m, z, t);
      return s.tot_acc + s.tot_conv + s.tot_ads;
    }
    double SourceA(const MmsCase& m, double z, double t) {
      const Terms s = TermsAt(m, z, t);
      return s.A_acc + s.A_conv + s.A_disp + s.A_ads;
    }
    double SourceB(const MmsCase& m, double z, double t) {
      const Terms s = TermsAt(m, z, t);
      return s.B_acc + s.B_conv + s.B_disp;
    }
    double SourceLoading(const MmsCase& m, double z, double t) {
      const Terms s = TermsAt(m, z, t);
      return s.q_acc + s.q_ldf;
    }
    // -dP/dz = (mu/K) u - S_m  =>  S_m = (mu/K) u + dP/dz.
    double SourceMomentum(const MmsCase& m, double z, double t) {
      const Terms s = TermsAt(m, z, t);
      return s.m_darcy + s.m_grad;
    }
    double SourceEnergy(const MmsCase& m, double z, double t) {
      const Terms s = TermsAt(m, z, t);
      return s.E_acc + s.E_conv + s.E_cond + s.E_work + s.E_ads;
    }

    // The same sources by central finite differences of the manufactured
    // fields alone (no hand-derived derivative), to check the derivation.
    struct FdSources { double tot, A, B, q, mom, E; };
    FdSources SourcesByFiniteDifferences(const MmsCase& m, double z, double t) {
      const double hz = 1e-4 * m.L, ht = 1e-4 * m.tau;
      auto C = [&](double zz, double tt) {
        const Fields f = Manufactured(m, zz, tt);
        return f.P / (kRgas * f.T);
      };
      auto T = [&](double zz, double tt) { return Manufactured(m, zz, tt).T; };
      auto y = [&](double zz, double tt) { return Manufactured(m, zz, tt).y; };
      auto u = [&](double zz, double tt) { return Manufactured(m, zz, tt).u; };
      auto q = [&](double zz, double tt) { return Manufactured(m, zz, tt).q; };
      auto P = [&](double zz, double tt) { return Manufactured(m, zz, tt).P; };
      auto dz = [&](const std::function<double(double, double)>& f) {
        return (f(z + hz, t) - f(z - hz, t)) / (2 * hz);
      };
      auto dt = [&](const std::function<double(double, double)>& f) {
        return (f(z, t + ht) - f(z, t - ht)) / (2 * ht);
      };
      // Dispersive flux eb D C dy/dz at z +- hz/2, differenced over hz.
      auto disp_flux = [&](double zz, double sign) {
        const double yz = (y(zz + 0.5 * hz, t) - y(zz - 0.5 * hz, t)) / hz;
        return m.eb * m.D * C(zz, t) * sign * yz;
      };
      const double dqdt = dt(q);
      FdSources s;
      s.tot = m.et() * dt(C) + dz([&](double zz, double tt) { return u(zz, tt) * C(zz, tt); }) + m.rho_b() * dqdt;
      s.A = m.et() * dt([&](double zz, double tt) { return C(zz, tt) * y(zz, tt); })
          + dz([&](double zz, double tt) { return u(zz, tt) * C(zz, tt) * y(zz, tt); })
          - (disp_flux(z + 0.5 * hz, 1.0) - disp_flux(z - 0.5 * hz, 1.0)) / hz
          + m.rho_b() * dqdt;
      s.B = m.et() * dt([&](double zz, double tt) { return C(zz, tt) * (1.0 - y(zz, tt)); })
          + dz([&](double zz, double tt) { return u(zz, tt) * C(zz, tt) * (1.0 - y(zz, tt)); })
          - (disp_flux(z + 0.5 * hz, -1.0) - disp_flux(z - 0.5 * hz, -1.0)) / hz;
      s.q = dqdt - m.k * (m.qstar(P(z, t) * y(z, t), T(z, t)) - q(z, t));
      s.mom = m.mu / m.Kperm() * u(z, t) + dz(P);
      const double Cp = y(z, t) * m.Cp_A + (1.0 - y(z, t)) * m.Cp_B;
      const double Cbed = m.et() * C(z, t) * Cp + m.rho_b() * (m.Cps + m.Cpa * q(z, t));
      const double Tzz = (T(z + hz, t) - 2.0 * T(z, t) + T(z - hz, t)) / (hz * hz);
      s.E = Cbed * dt(T) + C(z, t) * Cp * u(z, t) * dz(T) - m.lambda * Tzz
          - m.et() * dt(P) - m.dH * m.rho_b() * dqdt;
      return s;
    }

    // ---- The deck -------------------------------------------------------------
    // One Kozeny-Carman layer, isothermal, constant viscosity, constant
    // dispersion and LDF coefficient; a velocity inlet at the bottom (z = 0)
    // and a pressure outlet at the top (z = L).
    void BuildDeck(UserData& data, const MmsCase& m, int n_cells, AdvectionScheme scheme) {
      data.fluid.AddComponent("Nitrogen");  // B, inert
      data.fluid.AddComponent("Oxygen");    // A, adsorbing
      data.fluid.viscosity_mode_ = ViscosityMode::CONSTANT;
      data.fluid.viscosity_ = m.mu;
      data.fluid.gas_heat_capacity_mode_ = GasHeatCapacityMode::CONSTANT;
      data.fluid.gas_heat_capacity_["Oxygen"] = m.Cp_A;
      data.fluid.gas_heat_capacity_["Nitrogen"] = m.Cp_B;
      data.fluid.adsorbed_heat_capacity_mode_ = AdsorbedHeatCapacityMode::CONSTANT;
      data.fluid.adsorbed_heat_capacity_["Oxygen"] = m.Cpa;
      data.fluid.adsorbed_heat_capacity_["Nitrogen"] = 0.0;
      data.fluid.heat_of_adsorption_convention_ = HeatOfAdsorptionConvention::ISOSTERIC;

      auto ads = data.adsorbent_directory.Add("Sorbent");
      auto& adsorbent = data.adsorbent_directory.Get(ads);
      adsorbent.ep_ = m.ep;
      adsorbent.solid_density_ = m.rho_p;
      adsorbent.heat_capacity_ = m.Cps;
      adsorbent.thermal_conductivity_ = 0.0;
      auto isotherm = std::make_unique<LangmuirIsotherm>();
      isotherm->params_["Oxygen"].m1 = m.m1;
      isotherm->params_["Oxygen"].b1 = m.b1;
      isotherm->params_["Oxygen"].q1 = m.q1;
      adsorbent.adsorption_source_.AddIsotherm({"Oxygen"}, std::move(isotherm));
      auto ldf = std::make_unique<LinearDrivingForce>();
      auto mtc = std::make_unique<MassTransferCoefficientConstant>();
      mtc->constant_mtc_ = m.k;
      ldf->mtc_model_ = std::move(mtc);
      adsorbent.adsorption_source_.mass_transfer_models_["Oxygen"] = std::move(ldf);
      adsorbent.adsorption_source_.heat_of_adsorption_["Oxygen"] = m.thermal ? m.dH : 0.0;

      std::get<CylinderBed>(data.bed.shape_).diameter = m.diameter;
      auto lid = data.bed.layers_.Add("Sorbent");
      auto& layer = data.bed.layers_.Get(lid);
      layer.momentum_model_ = MomentumModel::KarmenCozeny;
      layer.domain_.length = m.L;
      layer.domain_.n_cells = n_cells;
      layer.eb_ = m.eb;
      layer.particle_diameter_ = m.dp;
      layer.composition_[ads] = 1.0;
      auto disp = std::make_unique<DispersionModelConstant>();
      disp->dispersion_coefficient_ = m.D;
      layer.equations_.dispersion_model_ = std::move(disp);
      // Adiabatic bed: no exchange with the wall.
      auto ihtc = std::make_unique<ConstantInternalHeatTransferModel>();
      ihtc->heat_transfer_coefficient_ = 0.0;
      layer.equations_.internal_heat_transfer_model_ = std::move(ihtc);
      auto ehtc = std::make_unique<ConstantExternalHeatTransferModel>();
      ehtc->heat_transfer_coefficient_ = 0.0;
      layer.equations_.external_heat_transfer_model_ = std::move(ehtc);
      auto heat_disp = std::make_unique<HeatDispersionModelConstant>();
      heat_disp->dispersion_coefficient_ = m.lambda;
      layer.equations_.heat_dispersion_model_ = std::move(heat_disp);
      if (!m.thermal) data.bed.equations_.energy_balance_mode_ = EnergyBalanceMode::ISOTHERMAL;

      data.ic.pressure = m.P_L;
      data.ic.temperature = m.T;
      data.ic.wall_temperature = m.T;
      data.ic.velocity = m.u_in;
      data.ic.mole_fraction["Nitrogen"] = 1.0 - m.y0;
      data.ic.mole_fraction["Oxygen"] = m.y0;

      data.cycle.AddStep("Feed");
      Step& step = data.cycle.GetStep("Feed");
      step.step_time_ = m.t_end;

      auto inlet = std::make_unique<BoundaryStream>();
      inlet->name_ = "Feed";
      auto flow = std::make_unique<VelocityBoundaryValue>();
      static_cast<ConstantBoundaryValue*>(flow->profile_.get())->constant_value_ = m.u_in;
      inlet->boundary_flow_ = std::move(flow);
      inlet->mole_fraction["Nitrogen"] = 1.0 - m.y0;
      inlet->mole_fraction["Oxygen"] = m.y0;
      inlet->temperature_ = m.T;
      auto inlet_bc = std::make_unique<FlowBoundaryCondition>();
      inlet_bc->direction_ = BoundaryDirection::INLET;
      inlet_bc->location_ = BoundaryLocation::BOTTOM;
      inlet_bc->streams_.push_back(std::move(inlet));
      step.bottom_boundary_condition_ = std::move(inlet_bc);

      auto outlet_bc = std::make_unique<PressureBoundaryCondition>();
      outlet_bc->direction_ = BoundaryDirection::OUTLET;
      outlet_bc->location_ = BoundaryLocation::TOP;
      outlet_bc->pressure_ = std::make_unique<ConstantBoundaryValue>(m.P_L);
      outlet_bc->mole_fraction_["Nitrogen"] = 1.0 - m.y0;
      outlet_bc->mole_fraction_["Oxygen"] = m.y0;
      step.top_boundary_condition_ = std::move(outlet_bc);

      data.cycle.number_of_cycles_ = 1;
      data.settings.reltol = 1e-10;
      data.settings.max_step = 0.02;
      data.settings.max_log_points_per_step_ = 4;
      data.settings.advection_scheme_ = scheme;

      // The verification hook: the manufactured sources (copies of m).
      auto& src = data.manufactured_sources_;
      src.total = [m](double z, double t) { return SourceTotal(m, z, t); };
      src.species["Oxygen"] = [m](double z, double t) { return SourceA(m, z, t); };
      src.species["Nitrogen"] = [m](double z, double t) { return SourceB(m, z, t); };
      src.loading["Oxygen"] = [m](double z, double t) { return SourceLoading(m, z, t); };
      src.momentum = [m](double z, double t) { return SourceMomentum(m, z, t); };
      if (m.thermal) src.energy = [m](double z, double t) { return SourceEnergy(m, z, t); };
    }

    // ---- Errors -----------------------------------------------------------------
    struct Errors {
      double y = 0, P = 0, q = 0, u = 0, T = 0;    // discrete L2 over the bed at t_end
      double y_max = 0;                            // max norm of y
      double y_end = 0;                            // |error| of y in the last cell (outlet end)
      double y_interior = 0;                       // L2 of y over 0.25 L < z < 0.75 L
      double t_logged = 0;
    };

    Errors Run(const MmsCase& m, int n_cells, AdvectionScheme scheme) {
      UserData data;
      BuildDeck(data, m, n_cells, scheme);
      std::string err;
      EXPECT_TRUE(data.PerformChecks(err)) << err;
      CoupledSolver solver{data};
      EXPECT_GE(solver.IntegrateCycle(data.cycle, false), 0)
        << "solver failed (scheme " << static_cast<int>(scheme) << ", " << n_cells << " cells)";

      const auto& log = solver.GetLog();
      const auto Plog = log.GetVariable(0, 0, VariableType::PRESSURE);
      const auto ulog = log.GetVariable(0, 0, VariableType::VELOCITY);
      const auto ylog = log.GetComponentVariable(0, 0, ComponentVariableType::MOLE_FRACTION, "Oxygen");
      const auto qlog = log.GetComponentVariable(0, 0, ComponentVariableType::ADSORBED_CONCENTRATION, "Oxygen");
      const auto Tlog = log.GetVariable(0, 0, VariableType::TEMPERATURE);
      Errors e;
      if (Plog.values.empty()) return e;
      const double t = Plog.times.back();
      e.t_logged = t;
      const auto& P = Plog.values.back();
      const auto& y = ylog.values.back();
      const auto& q = qlog.values.back();
      const auto& u = ulog.values.back();
      const auto& T = Tlog.values.back();
      const double dz = m.L / n_cells;
      // Cell n (1..N) has its centre at (n - 1/2) dz; index 0 and N + 1 are the ghosts.
      double sy = 0, sP = 0, sq = 0, su = 0, sT = 0, sy_int = 0;
      int n_int = 0;
      for (int n = 1; n <= n_cells; ++n) {
        const double z = (n - 0.5) * dz;
        const Fields f = Manufactured(m, z, t);
        const double ey = y[n] - f.y;
        sy += ey * ey;
        sP += (P[n] - f.P) * (P[n] - f.P);
        sq += (q[n] - f.q) * (q[n] - f.q);
        sT += (T[n] - f.T) * (T[n] - f.T);
        e.y_max = std::max(e.y_max, std::abs(ey));
        if (z > 0.25 * m.L && z < 0.75 * m.L) { sy_int += ey * ey; ++n_int; }
      }
      // Faces 0..N at z = i dz (face 0 is the inlet, N the outlet).
      for (int i = 0; i <= n_cells; ++i) {
        const Fields f = Manufactured(m, i * dz, t);
        su += (u[i] - f.u) * (u[i] - f.u);
      }
      e.y = std::sqrt(sy / n_cells);
      e.P = std::sqrt(sP / n_cells);
      e.q = std::sqrt(sq / n_cells);
      e.u = std::sqrt(su / (n_cells + 1));
      e.T = std::sqrt(sT / n_cells);
      e.y_interior = std::sqrt(sy_int / std::max(n_int, 1));
      e.y_end = std::abs(y[n_cells] - Manufactured(m, (n_cells - 0.5) * dz, t).y);
      return e;
    }

    double Order(double coarse, double fine) { return std::log2(coarse / fine); }

    const char* SchemeName(AdvectionScheme s) {
      switch (s) {
        case AdvectionScheme::WENO: return "WENO";
        case AdvectionScheme::VanLeer: return "van Leer";
        default: return "upwind";
      }
    }

    void MeasureOrders(const MmsCase& m, AdvectionScheme scheme, std::vector<Errors>& out) {
      for (int n : {20, 40, 80, 160}) out.push_back(Run(m, n, scheme));
      for (size_t k = 0; k < out.size(); ++k) {
        const auto& e = out[k];
        std::printf("[MMS %s%s] N=%3d  y %.3e (int %.3e, max %.3e, end %.3e)  P %.3e  q %.3e  u %.3e  T %.3e",
                    SchemeName(scheme), m.thermal ? " thermal" : "", 20 << k, e.y, e.y_interior,
                    e.y_max, e.y_end, e.P, e.q, e.u, e.T);
        if (k > 0) {
          const auto& c = out[k - 1];
          std::printf("  | order y %.2f int %.2f max %.2f end %.2f P %.2f q %.2f u %.2f T %.2f",
                      Order(c.y, e.y), Order(c.y_interior, e.y_interior), Order(c.y_max, e.y_max),
                      Order(c.y_end, e.y_end), Order(c.P, e.P), Order(c.q, e.q), Order(c.u, e.u),
                      m.thermal ? Order(c.T, e.T) : 0.0);
        }
        std::printf("\n");
      }
    }

  } // namespace

  // ===========================================================================
  // The derivation: the hand-derived sources equal finite differences of the
  // manufactured fields, and every term is large.
  // ===========================================================================

  TEST(ManufacturedSolutions, SourcesMatchFiniteDifferencesAndEveryTermIsActive) {
    for (bool thermal : {false, true}) {
      MmsCase m;
      m.thermal = thermal;
      double max_rel = 0.0;
      Terms peak{};
      double Cbed_max = 0.0;
      auto track = [](double& p, double v) { p = std::max(p, std::abs(v)); };
      for (double t : {0.5, 1.5, 3.0, 4.5, 6.0}) {
        for (int k = 1; k < 20; ++k) {
          const double z = m.L * k / 20.0;
          const FdSources fd = SourcesByFiniteDifferences(m, z, t);
          const double S[6] = {SourceTotal(m, z, t), SourceA(m, z, t), SourceB(m, z, t),
                               SourceLoading(m, z, t), SourceMomentum(m, z, t), SourceEnergy(m, z, t)};
          const double F[6] = {fd.tot, fd.A, fd.B, fd.q, fd.mom, fd.E};
          // Relative to the largest term of that balance, the scale a wrong
          // term would be measured against.
          const Terms s = TermsAt(m, z, t);
          const double scale[6] = {
            std::max({std::abs(s.tot_acc), std::abs(s.tot_conv), std::abs(s.tot_ads)}),
            std::max({std::abs(s.A_acc), std::abs(s.A_conv), std::abs(s.A_disp), std::abs(s.A_ads)}),
            std::max({std::abs(s.B_acc), std::abs(s.B_conv), std::abs(s.B_disp)}),
            std::max(std::abs(s.q_acc), std::abs(s.q_ldf)),
            std::max(std::abs(s.m_darcy), std::abs(s.m_grad)),
            std::max({std::abs(s.E_acc), std::abs(s.E_conv), std::abs(s.E_cond), std::abs(s.E_work),
                      std::abs(s.E_ads), 1.0})};
          for (int i = 0; i < 6; ++i)
            max_rel = std::max(max_rel, std::abs(S[i] - F[i]) / scale[i]);
          track(peak.tot_acc, s.tot_acc); track(peak.tot_conv, s.tot_conv); track(peak.tot_ads, s.tot_ads);
          track(peak.A_acc, s.A_acc); track(peak.A_conv, s.A_conv); track(peak.A_disp, s.A_disp);
          track(peak.A_ads, s.A_ads); track(peak.q_acc, s.q_acc); track(peak.q_ldf, s.q_ldf);
          track(peak.m_darcy, s.m_darcy); track(peak.m_grad, s.m_grad); track(peak.Kp, s.Kp);
          track(peak.E_acc, s.E_acc); track(peak.E_conv, s.E_conv); track(peak.E_cond, s.E_cond);
          track(peak.E_work, s.E_work); track(peak.E_ads, s.E_ads); track(Cbed_max, s.Cbed);
          // The overall balance is the sum of the component balances.
          EXPECT_NEAR(SourceA(m, z, t) + SourceB(m, z, t), SourceTotal(m, z, t),
                      1e-12 * std::abs(SourceTotal(m, z, t)) + 1e-12);
        }
      }
      // Central differences with step 1e-4 (of L or tau) are exact to
      // O(1e-8) relative and lose ~1e-16 / 1e-4 to rounding: a slip in the
      // hand derivation (a dropped or mis-signed term, each at least 15 % of
      // its balance below) would be 1e-1.
      EXPECT_LT(max_rel, 1e-5) << "hand-derived source disagrees with finite differences";

      // Guards: in each gas balance every term is at least 15 % of the
      // largest, so none of them can be wrong without moving the solution.
      const double gas = std::max({peak.A_acc, peak.A_conv, peak.A_disp, peak.A_ads});
      for (double term : {peak.A_acc, peak.A_conv, peak.A_disp, peak.A_ads})
        EXPECT_GT(term, 0.15 * gas);
      const double tot = std::max({peak.tot_acc, peak.tot_conv, peak.tot_ads});
      for (double term : {peak.tot_acc, peak.tot_conv, peak.tot_ads})
        EXPECT_GT(term, 0.15 * tot);
      EXPECT_GT(peak.q_ldf, peak.q_acc);  // the isotherm drives the loading source
      EXPECT_GT(peak.m_darcy, 1e4);       // [Pa/m]: both momentum terms large...
      EXPECT_GT(peak.m_grad, 1e4);
      EXPECT_GT(peak.Kp, 1.5);            // ...and the Langmuir isotherm well past its Henry limit
      if (thermal) {
        // The solid heat capacity dominates the energy balance, as in a real
        // bed, so its terms are compared by the temperature each would move
        // over the run if it were dropped, term * t_end / Cbed: at least
        // 5e-3 K (measured: work 8e-3, conduction 1.7e-2, convection 3.9e-2,
        // adsorption 0.17), a hundred times the temperature error ThermalOrder
        // allows the high-order schemes at 160 cells (5e-5 K).
        for (double term : {peak.E_conv, peak.E_cond, peak.E_work, peak.E_ads})
          EXPECT_GT(term * m.t_end / Cbed_max, 5e-3);
      }
      std::printf("[MMS terms%s] tot %.3g %.3g %.3g | A %.3g %.3g %.3g %.3g | q %.3g %.3g | mom %.3g %.3g"
                  " | E %.3g %.3g %.3g %.3g %.3g (Cbed %.3g) | Kp %.2f | FD %.1e\n",
                  thermal ? " thermal" : "", peak.tot_acc, peak.tot_conv, peak.tot_ads, peak.A_acc,
                  peak.A_conv, peak.A_disp, peak.A_ads, peak.q_acc, peak.q_ldf, peak.m_darcy, peak.m_grad,
                  peak.E_acc, peak.E_conv, peak.E_cond, peak.E_work, peak.E_ads, Cbed_max, peak.Kp, max_rel);
    }
  }

  // The hook changes nothing when its sources are absent: a run with every
  // source set to return zero is bit-identical to a run without the hook,
  // which is the code path of every deck (the residuals only branch on an
  // empty source).
  TEST(ManufacturedSolutions, ZeroSourcesAreBitIdentical) {
    MmsCase m;
    m.thermal = true;
    auto run = [&](bool zero_hook) {
      UserData data;
      BuildDeck(data, m, 20, AdvectionScheme::Upwind);
      data.manufactured_sources_ = ManufacturedSources{};
      if (zero_hook) {
        auto zero = [](double, double) { return 0.0; };
        auto& src = data.manufactured_sources_;
        src.total = zero;
        src.species["Oxygen"] = zero;
        src.species["Nitrogen"] = zero;
        src.loading["Oxygen"] = zero;
        src.momentum = zero;
        src.energy = zero;
        EXPECT_TRUE(src.Any());
      } else {
        EXPECT_FALSE(data.manufactured_sources_.Any());
      }
      std::string err;
      EXPECT_TRUE(data.PerformChecks(err)) << err;
      CoupledSolver solver{data};
      EXPECT_GE(solver.IntegrateCycle(data.cycle, false), 0);
      const auto& log = solver.GetLog();
      std::vector<std::vector<double>> out;
      for (auto type : {VariableType::PRESSURE, VariableType::VELOCITY, VariableType::TEMPERATURE})
        for (const auto& v : log.GetVariable(0, 0, type).values) out.push_back(v);
      for (const auto& v : log.GetComponentVariable(0, 0, ComponentVariableType::MOLE_FRACTION, "Oxygen").values)
        out.push_back(v);
      return out;
    };
    const auto plain = run(false);
    const auto zero = run(true);
    ASSERT_FALSE(plain.empty());
    EXPECT_EQ(plain, zero);
  }

  // ===========================================================================
  // Observed orders
  // ===========================================================================

  // Orders on 20/40/80/160 cells; the bounds come from the measured
  // asymptotic range (README; the same case on 320 and 640 cells gives
  // upwind 1.00, van Leer 1.97-2.03 and WENO 1.88-2.23 in every field).
  namespace {
    void ExpectOrders(const std::vector<Errors>& e, double lo, double hi, double lo_max, const char* what) {
      ASSERT_EQ(e.size(), 4u);
      for (size_t k = 2; k < e.size(); ++k) {  // the 20 -> 40 -> 80 and 40 -> 80 -> 160 pairs
        const auto& c = e[k - 1];
        const auto& f = e[k];
        SCOPED_TRACE(std::string(what) + ", " + std::to_string(20 << (k - 1)) + " -> "
                     + std::to_string(20 << k) + " cells");
        for (auto [name, p] : {std::pair{"y", Order(c.y, f.y)}, std::pair{"P", Order(c.P, f.P)},
                               std::pair{"q", Order(c.q, f.q)}, std::pair{"u", Order(c.u, f.u)}}) {
          EXPECT_GT(p, lo) << name;
          EXPECT_LT(p, hi) << name;
        }
        EXPECT_GT(Order(c.y_max, f.y_max), lo_max) << "max norm of y";
        EXPECT_GT(Order(c.y_end, f.y_end), lo_max - 0.3) << "y in the outlet cell";
      }
    }
  }

  TEST(ManufacturedSolutions, UpwindOrder) {
    std::vector<Errors> e;
    MeasureOrders(MmsCase{}, AdvectionScheme::Upwind, e);
    // Measured 0.95-1.00.
    ExpectOrders(e, 0.93, 1.07, 0.93, "upwind");
  }

  TEST(ManufacturedSolutions, WenoOrder) {
    std::vector<Errors> e;
    MeasureOrders(MmsCase{}, AdvectionScheme::WENO, e);
    // Measured: L2 1.81-2.51 (P lowest, u highest), max 1.90-1.94, outlet
    // cell 1.59-1.72 (rising to 1.88 on 640 cells).
    ExpectOrders(e, 1.7, 2.7, 1.8, "WENO");
    // The third-order interior reconstruction: on these grids the interior
    // error converges faster than second order (2.55-2.69, falling to 2.24
    // on 640 cells as the O(dz^2) central terms take over) and, from 40
    // cells on, below van Leer's (1.7 to 4 times; equal at 20 cells).
    for (size_t k = 2; k < e.size(); ++k)
      EXPECT_GT(Order(e[k - 1].y_interior, e[k].y_interior), 2.3);
    std::vector<Errors> vl;
    MeasureOrders(MmsCase{}, AdvectionScheme::VanLeer, vl);
    for (size_t k = 1; k < e.size(); ++k) EXPECT_LT(e[k].y_interior, vl[k].y_interior);
  }

  TEST(ManufacturedSolutions, VanLeerOrder) {
    std::vector<Errors> e;
    MeasureOrders(MmsCase{}, AdvectionScheme::VanLeer, e);
    // Measured: L2 1.86-2.10, max and outlet cell 1.83-1.90.
    ExpectOrders(e, 1.75, 2.3, 1.75, "van Leer");
  }

  // Non-isothermal: the energy balance's convective face temperature is
  // reconstructed with the same scheme as the mole fractions, so the
  // temperature converges at order 1 with upwind and 2 with WENO and van Leer,
  // and the other fields keep the orders of the isothermal case.
  //
  // GUARD: when the face temperature was always the upwind value (before the
  // fix), T converged at 0.97-1.05 for all three schemes and its error pulled
  // the high-order schemes' y order down to 1.3-1.5 between 80 and 160 cells
  // (about 1 on 640), far below the bounds here.
  TEST(ManufacturedSolutions, ThermalOrder) {
    MmsCase m;
    m.thermal = true;
    std::vector<Errors> up, weno, vl;
    MeasureOrders(m, AdvectionScheme::Upwind, up);
    MeasureOrders(m, AdvectionScheme::WENO, weno);
    MeasureOrders(m, AdvectionScheme::VanLeer, vl);
    // Every field: measured upwind 0.95-1.00; WENO L2 2.08-2.69, max
    // 1.72-2.03, outlet cell 1.46-1.72; van Leer L2 1.86-2.43, max and outlet
    // cell 1.71-1.90 (the isothermal bounds hold).
    ExpectOrders(up, 0.93, 1.07, 0.93, "upwind, thermal");
    ExpectOrders(weno, 1.7, 2.7, 1.7, "WENO, thermal");
    ExpectOrders(vl, 1.75, 2.5, 1.7, "van Leer, thermal");
    // The temperature: measured upwind 0.97-0.99, WENO 2.06-2.25, van Leer
    // 1.95-2.03 (on 160 -> 320 -> 640 cells: upwind 1.00, WENO 2.19 and
    // 1.99, van Leer 2.03 and 1.99).
    auto expect_T_order = [](const std::vector<Errors>& e, double lo, double hi, const char* what) {
      for (size_t k = 2; k < e.size(); ++k) {
        const double p = Order(e[k - 1].T, e[k].T);
        EXPECT_GT(p, lo) << what << ", T, " << (20 << (k - 1)) << " -> " << (20 << k) << " cells";
        EXPECT_LT(p, hi) << what << ", T, " << (20 << (k - 1)) << " -> " << (20 << k) << " cells";
      }
    };
    expect_T_order(up, 0.9, 1.1, "upwind");
    expect_T_order(weno, 1.85, 2.5, "WENO");
    expect_T_order(vl, 1.85, 2.3, "van Leer");
    // The high-order schemes' errors stay far below upwind's (measured at
    // 160 cells: T 1.6-2.1e-5 against 1.7e-3 K, y 0.7-1.0e-6 against 9.0e-5).
    for (const auto* e : {&weno, &vl}) {
      EXPECT_LT(e->back().T, 5e-5);
      EXPECT_LT(e->back().T, 0.03 * up.back().T);
      EXPECT_LT(e->back().y, 0.05 * up.back().y);
    }
  }

} // namespace adsim
