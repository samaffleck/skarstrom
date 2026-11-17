/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Axial dispersion with a coefficient per component (the Ruthven and
// Taylor-Aris models).
//
// Each component i disperses with its own D_ax,i, built on its mixture-averaged
// molecular diffusivity D_m,i (Chapman-Enskog binary D_ij with Blanc's law).
// Fickian fluxes -eb D_ax,i C dy_i/dz with different coefficients do not sum
// to zero, and the overall balance has no dispersion term, so the component
// balance takes the drift-corrected flux (Ramshaw, J. Non-Equilib. Thermodyn.
// 15 (1990) 295, in the molar-average frame):
//
//   J_i = -eb C (D_ax,i dy_i/dz - y_i sum_j D_ax,j dy_j/dz),  sum_i J_i = 0.
//
// The checks:
//  1. Residual level, exact: in a still column (no flow, no accumulation, no
//     sources) with a ternary composition profile across a layer interface,
//     the residuals hold the face-flux differences alone. Summed over the
//     components in each cell they vanish (the fluxes sum to zero), and
//     summed over the cells for each component they vanish (each face flux is
//     single valued). Also when the mole fractions do not sum exactly to 1.
//  2. Residual level, binary: with two components D_m,1 = D_m,2 = D_12, so
//     the two components share one coefficient and the residuals equal those
//     of the Constant model at that coefficient (the correction is then zero
//     to round-off).
//  3. Solved, trace limit: H2 and CO2 at 1 % each in N2 through an inert bed.
//     The correction is second order in the trace fraction, so each tracer's
//     outlet curve is the Ogata-Banks curve of its own D_ax,i; H2 spreads more
//     than CO2.
//  4. Solved, comparable fractions: a 40 % H2, 30 % CO2 feed into N2. The
//     mole fractions keep summing to 1 and the logged streams account for
//     every component's inventory change.
//  5. Registry and deck: the Ruthven and Taylor-Aris models have no
//     molecular-diffusivity input, and a deck that stored one (cereal
//     versions 1 and 2) loads.

#include "validation_helpers.h"

#include "adsim/parameters/parameter_registry.h"
#include "adsim/solver/kpis.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <numeric>
#include <string>
#include <vector>

namespace adsim {

  using namespace adsim::validation;

  namespace {

    const std::string kN2 = "Nitrogen";
    const std::string kH2 = "Hydrogen";
    const std::string kCO2 = "CarbonDioxide";

    constexpr double kT0 = 298.0;      // [K]
    constexpr double kP0 = 1.0e5;      // [Pa]
    constexpr double kDiameter = 0.05; // [m]

    struct LayerSpec {
      double length;  // [m]
      int cells;
      double eb;      // [-]
    };

    struct Column {
      std::vector<std::string> components{kN2, kH2, kCO2};
      std::vector<LayerSpec> layers{{0.15, 6, 0.35}, {0.15, 15, 0.48}};
      // A Constant model of this coefficient instead of `kind` (0: `kind`).
      double constant_dispersion = 0.0;  // [m2/s]
      AxialDispersionType kind = AxialDispersionType::RUTHVEN;
      double particle_diameter = 0.002;
    };

    // An isothermal column of inert packing (no isotherm), each layer with
    // Ruthven dispersion (or a constant coefficient).
    void BuildColumn(UserData& data, const Column& c) {
      for (const auto& comp : c.components) data.fluid.AddComponent(comp);
      data.fluid.viscosity_mode_ = ViscosityMode::CONSTANT;
      data.fluid.viscosity_ = 1.8e-5;
      data.bed.equations_.energy_balance_mode_ = EnergyBalanceMode::ISOTHERMAL;
      std::get<CylinderBed>(data.bed.shape_).diameter = kDiameter;
      data.bed.wall_.thickness_ = 0.004;
      data.bed.wall_.density_ = 7900.0;
      data.bed.wall_.heat_capacity_ = 500.0;
      data.bed.wall_.thermal_conductivity_ = 16.0;

      for (size_t k = 0; k < c.layers.size(); ++k) {
        const auto& s = c.layers[k];
        const auto ads_id = data.adsorbent_directory.Add("Packing " + std::to_string(k + 1));
        auto& ads = data.adsorbent_directory.Get(ads_id);
        ads.solid_density_ = 1000.0;
        ads.heat_capacity_ = 900.0;
        ads.thermal_conductivity_ = 0.0;
        ads.ep_ = 0.0;

        const auto lid = data.bed.layers_.Add("Layer " + std::to_string(k + 1));
        auto& layer = data.bed.layers_.Get(lid);
        layer.momentum_model_ = MomentumModel::KarmenCozeny;
        layer.domain_.length = s.length;
        layer.domain_.n_cells = s.cells;
        layer.eb_ = s.eb;
        layer.particle_diameter_ = c.particle_diameter;
        layer.composition_[ads_id] = 1.0;
        if (c.constant_dispersion > 0.0) {
          auto disp = std::make_unique<DispersionModelConstant>();
          disp->dispersion_coefficient_ = c.constant_dispersion;
          layer.equations_.dispersion_model_ = std::move(disp);
        } else {
          layer.equations_.dispersion_model_ = CreateAxialDispersionModel(c.kind);
        }
        auto heat_disp = std::make_unique<HeatDispersionModelConstant>();
        heat_disp->dispersion_coefficient_ = 0.0;
        layer.equations_.heat_dispersion_model_ = std::move(heat_disp);
        auto ihtc = std::make_unique<ConstantInternalHeatTransferModel>();
        ihtc->heat_transfer_coefficient_ = 0.0;
        layer.equations_.internal_heat_transfer_model_ = std::move(ihtc);
        auto ehtc = std::make_unique<ConstantExternalHeatTransferModel>();
        ehtc->heat_transfer_coefficient_ = 0.0;
        layer.equations_.external_heat_transfer_model_ = std::move(ehtc);
      }

      data.ic.pressure = kP0;
      data.ic.temperature = kT0;
      data.ic.wall_temperature = kT0;
      data.ic.velocity = 0.0;
      for (const auto& comp : c.components) data.ic.mole_fraction[comp] = (comp == kN2) ? 1.0 : 0.0;
      data.cycle.number_of_cycles_ = 1;
    }

    // =========================================================================
    // Residual level
    // =========================================================================

    // The species residuals of a still column (u = 0, every time derivative 0,
    // no sources) whose composition varies along it, in mol/(m2 s) per cell:
    // residual[component][cell]. Every component's profile is a front centred
    // on the layer interface; `sum_error` scales the composition so that the
    // mole fractions sum to 1 + sum_error, as the solver allows.
    struct StillColumn {
      std::vector<std::vector<double>> residual;  // per component, per cell
      std::vector<double> width;
      int interface_cell = 0;  // first cell of the second layer
      double max_drift = 0.0;  // max_f |S_f| / sum_i |k_i,f dy_i| over the interior faces
      std::vector<std::vector<double>> dax;  // D_ax,i per component, per cell
    };

    StillColumn StillColumnResiduals(const Column& c, double sum_error = 0.0) {
      UserData data;
      BuildColumn(data, c);
      data.cycle.AddStep("Hold");
      Step& step = data.cycle.GetStep("Hold");
      step.step_time_ = 1.0;
      auto bottom = std::make_unique<ClosedBoundaryCondition>();
      bottom->location_ = BoundaryLocation::BOTTOM;
      step.bottom_boundary_condition_ = std::move(bottom);
      auto top = std::make_unique<ClosedBoundaryCondition>();
      top->location_ = BoundaryLocation::TOP;
      step.top_boundary_condition_ = std::move(top);
      std::string error;
      EXPECT_TRUE(data.EvaluateAllExpressions(error) && data.PerformChecks(error)) << error;
      data.Initialise();
      data.cycle.step_list_.front().Initialise(data.state);

      auto& s = data.state;
      const auto& dom = s.domain_vectors;
      const int N = static_cast<int>(s.pressure.var.size());
      const size_t nc = c.components.size();

      StillColumn out;
      out.interface_cell = 1 + c.layers.front().cells;
      double z_if = 0.0;
      for (int i = 1; i < out.interface_cell; ++i) z_if += dom.cell_width[i];

      // Fronts of different widths and centres, so every gradient differs.
      double z = 0.0;
      for (int i = 0; i < N; ++i) {
        const double zc = (i == 0) ? 0.0 : (i == N - 1 ? z : z + 0.5 * dom.cell_width[i]);
        if (i > 0 && i < N - 1) z += dom.cell_width[i];
        std::vector<double> y(nc);
        double rest = 1.0;
        for (size_t k = 1; k < nc; ++k) {
          const double w = 0.03 + 0.02 * k;
          y[k] = (0.45 / k) * (0.5 + 0.5 * std::tanh((zc - z_if + 0.01 * k) / w));
          rest -= y[k];
        }
        y[0] = rest;
        for (size_t k = 0; k < nc; ++k)
          s.mole_fraction.at(c.components[k]).var[i] = y[k] * (1.0 + sum_error);
        s.pressure.var[i] = kP0;
        s.fluid_temperature.var[i] = kT0;
      }
      for (const auto& comp : c.components) {
        auto& v = s.mole_fraction.at(comp).var;
        v[0] = v[1];
        v[N - 1] = v[N - 2];
      }
      std::fill(s.velocity.var.begin(), s.velocity.var.end(), 0.0);

      s.UpdateIndependentEquations(&data);
      for (const auto& comp : c.components) {
        UpdateMoleFractionFace(&data, comp);
        UpdateMoleFraction(&data, comp);
      }

      out.width.assign(dom.cell_width.begin(), dom.cell_width.end());
      for (const auto& comp : c.components) {
        std::vector<double> r(N, 0.0);
        for (int i = 1; i < N - 1; ++i)
          r[i] = s.mole_fraction.at(comp).res_var[i] * s.pressure.var[i] / (kR * s.fluid_temperature.var[i]);
        out.residual.push_back(r);
        out.dax.push_back(s.PerComponentDispersion() ? s.component_axial_dispersion.at(comp).var
                                                     : s.axial_dispersion_coefficient.var);
      }
      if (s.PerComponentDispersion()) {
        // The size of the drift the correction removes, relative to the
        // uncorrected fluxes through the same face.
        for (int f = 1; f < N - 2; ++f) {
          double mag = 0.0;
          for (const auto& comp : c.components) {
            const auto& y = s.mole_fraction.at(comp).var;
            mag += std::abs(s.component_axial_dispersion.at(comp).var[f] * (y[f + 1] - y[f]));
          }
          if (mag > 0.0) out.max_drift = std::max(out.max_drift, std::abs(s.dispersion_drift.var[f]) / mag);
        }
      }
      return out;
    }

  } // namespace

  // 1. The dispersive fluxes of the components sum to zero through every face,
  //    so in each cell the species residuals sum to zero; and each face flux is
  //    single valued, so each component's residuals sum to zero over the
  //    column. Both to round-off, also across a layer interface where the cell
  //    width and eb jump and when the mole fractions sum to 1 + 1e-6.
  //
  // Measured: cell sums at most 1.6e-15, species sums at most 1.1e-16 of the
  // summed magnitudes.
  //
  // GUARD: the coefficients really differ (D_ax,H2 is more than three times
  // D_ax,CO2 in a still column, where D_ax,i = 0.7 D_m,i), and the drift
  // sum_j k_j,f dy_j that the correction removes is a sizeable fraction of the
  // uncorrected fluxes, so a missing correction would fail the cell sums.
  TEST(PerComponentDispersion, StillColumnFluxesSumToZeroInEveryCell) {
    for (double sum_error : {0.0, 1e-6}) {
      Column c;
      const StillColumn r = StillColumnResiduals(c, sum_error);
      const int N = static_cast<int>(r.width.size());
      const int mid = r.interface_cell + 2;
      EXPECT_GT(r.dax[1][mid], 3.0 * r.dax[2][mid]) << "H2 and CO2 disperse alike";
      // Measured 0.14: without the correction a cell would gain or lose 14 %
      // of its dispersive flux magnitude, eleven orders above the bound below.
      EXPECT_GT(r.max_drift, 0.1) << "the correction has nothing to correct";

      double worst_cell = 0.0;
      for (int i = 1; i < N - 1; ++i) {
        double net = 0.0, mag = 0.0;
        for (const auto& res : r.residual) {
          net += res[i];
          mag += std::abs(res[i]);
        }
        if (mag > 0.0) worst_cell = std::max(worst_cell, std::abs(net) / mag);
      }
      std::vector<double> worst_species;
      for (const auto& res : r.residual) {
        double net = 0.0, mag = 0.0;
        for (int i = 1; i < N - 1; ++i) {
          net += res[i];
          mag += std::abs(res[i]);
        }
        EXPECT_GT(mag, 0.0) << "a component does not disperse";
        worst_species.push_back(std::abs(net) / mag);
      }
      std::printf("[still column, sum error %.0e] cell sum %.2e, species sums N2 %.2e H2 %.2e CO2 %.2e, drift %.2f\n",
                  sum_error, worst_cell, worst_species[0], worst_species[1], worst_species[2], r.max_drift);
      EXPECT_LT(worst_cell, 1e-12) << "the dispersive fluxes do not sum to zero in a cell";
      for (size_t k = 0; k < worst_species.size(); ++k)
        EXPECT_LT(worst_species[k], 1e-12) << c.components[k] << ": a face flux is not single valued";
    }
  }

  // 2. Binary: Blanc's law gives D_m,1 = D_m,2 = D_12, so Ruthven's residuals
  //    are those of a constant coefficient 0.7 D_12 (no flow in a still
  //    column). Uniform T and P make D_12 uniform; it is read from the run.
  TEST(PerComponentDispersion, BinaryEqualsOneCoefficientAtD12) {
    Column mixed;
    mixed.components = {kN2, kCO2};
    const StillColumn a = StillColumnResiduals(mixed);
    const int mid = a.interface_cell + 2;
    // D_ax = 0.7 D_m in a still column.
    const double d12 = a.dax[0][mid] / 0.7;
    EXPECT_NEAR(a.dax[0][mid], a.dax[1][mid], 1e-15 * a.dax[0][mid]);

    Column shared = mixed;
    shared.constant_dispersion = 0.7 * d12;
    const StillColumn b = StillColumnResiduals(shared);

    double worst = 0.0, scale = 0.0;
    for (size_t k = 0; k < a.residual.size(); ++k)
      for (size_t i = 0; i < a.residual[k].size(); ++i) {
        worst = std::max(worst, std::abs(a.residual[k][i] - b.residual[k][i]));
        scale = std::max(scale, std::abs(b.residual[k][i]));
      }
    std::printf("[binary] D_12 %.4e m2/s, max residual difference %.2e of the largest\n", d12, worst / scale);
    EXPECT_GT(scale, 0.0);
    EXPECT_LT(worst, 1e-12 * scale);
  }

  // =========================================================================
  // Solved breakthroughs
  // =========================================================================

  namespace {

    constexpr double kU = 0.005;       // [m/s] interstitial velocity
    constexpr double kLength = 0.3;    // [m]
    constexpr int kCells = 300;
    constexpr double kEb = 0.4;
    constexpr double kDp = 0.001;      // [m]

    // One inert layer, a constant molar feed at the bottom, a pressure outlet
    // at the top. Upwind, so the numerical dispersion is the known u dx / 2.
    // u_i is the interstitial velocity of the feed.
    void BuildBreakthrough(UserData& data, const std::vector<std::pair<std::string, double>>& feed,
                           double run_time, int cells, AxialDispersionType kind = AxialDispersionType::RUTHVEN,
                           double u_i = kU, double length = kLength) {
      Column c;
      c.layers = {{length, cells, kEb}};
      c.particle_diameter = kDp;
      c.kind = kind;
      BuildColumn(data, c);
      data.settings.advection_scheme_ = AdvectionScheme::Upwind;

      const double area = 0.25 * std::acos(-1.0) * kDiameter * kDiameter;
      const double molar_flow = kEb * u_i * area * kP0 / (kR * kT0);

      data.cycle.AddStep("Feed");
      Step& step = data.cycle.GetStep("Feed");
      step.step_time_ = run_time;
      auto stream = std::make_unique<BoundaryStream>();
      auto flow = std::make_unique<MolarFlowBoundaryValue>();
      static_cast<ConstantBoundaryValue*>(flow->profile_.get())->constant_value_ = molar_flow;
      stream->boundary_flow_ = std::move(flow);
      for (const auto& [comp, y] : feed) stream->mole_fraction[comp] = y;
      stream->temperature_ = kT0;
      auto inlet = std::make_unique<FlowBoundaryCondition>();
      inlet->direction_ = BoundaryDirection::INLET;
      inlet->location_ = BoundaryLocation::BOTTOM;
      inlet->streams_.push_back(std::move(stream));
      step.bottom_boundary_condition_ = std::move(inlet);

      auto outlet = std::make_unique<PressureBoundaryCondition>();
      outlet->direction_ = BoundaryDirection::OUTLET;
      outlet->location_ = BoundaryLocation::TOP;
      auto p_out = std::make_unique<ConstantBoundaryValue>();
      p_out->constant_value_ = kP0;
      outlet->pressure_ = std::move(p_out);
      for (const auto& [comp, y] : feed) outlet->mole_fraction_[comp] = y;
      step.top_boundary_condition_ = std::move(outlet);

      data.settings.reltol = 1e-9;
      data.settings.max_step = 0.2;
    }

    // Outlet mole fraction of `component` over the step (the top face).
    std::pair<std::vector<double>, std::vector<double>> Outlet(const EquationLog& log, const std::string& component) {
      const auto y = log.GetComponentFaceVariable(0, 0, ComponentVariableType::MOLE_FRACTION, component);
      std::vector<double> series;
      for (const auto& snapshot : y.values) series.push_back(snapshot.empty() ? 0.0 : snapshot.back());
      return {y.times, series};
    }

    double CurveMean(const std::vector<double>& t, const std::vector<double>& F) {
      double m1 = 0.0;
      for (size_t i = 1; i < t.size(); ++i) m1 += 0.5 * ((1.0 - F[i - 1]) + (1.0 - F[i])) * (t[i] - t[i - 1]);
      return m1;
    }

    // RMS of F - OgataBanks(D) over [t_lo, t_hi], the analytic curve shifted
    // by t_shift (the finite-column offset of the front position), as in
    // test_analytic_dispersion.cpp.
    double ShiftedRms(const std::vector<double>& t, const std::vector<double>& F, double length, double u, double D,
                      double t_shift, double t_lo, double t_hi) {
      double sum = 0.0;
      int n = 0;
      for (size_t i = 0; i < t.size(); ++i) {
        if (t[i] < t_lo || t[i] > t_hi) continue;
        const double d = F[i] - OgataBanks(length, u, D, t[i] - t_shift);
        sum += d * d;
        ++n;
      }
      return n > 0 ? std::sqrt(sum / n) : 1.0;
    }


    // A trace breakthrough: H2 and CO2 at 1 % each in N2 through an inert bed
    // with `kind` dispersion. Each tracer's outlet curve against Ogata-Banks
    // at its own D_ax,i + u dx / 2 (upwind), and at the other tracer's.
    struct TraceResult {
      double u, d_num, dm_h2, dm_co2, d_h2, d_co2;
      double own_h2, other_h2, own_co2, other_co2;
    };

    TraceResult RunTrace(AxialDispersionType kind, double u_i, double length, int cells, UserData& data) {
      TraceResult r{};
      const double t_front = length / u_i;
      BuildBreakthrough(data, {{kN2, 0.98}, {kH2, 0.01}, {kCO2, 0.01}}, 2.5 * t_front, cells, kind, u_i, length);
      CoupledSolver solver{data};
      std::string error;
      EXPECT_TRUE(data.PerformChecks(error)) << error;
      EXPECT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";
      EXPECT_TRUE(data.state.PerComponentDispersion());

      const auto& log = solver.GetLog();
      const auto u_face = log.GetVariable(0, 0, VariableType::VELOCITY).values.back();  // a copy: the entry is a temporary
      double u_sum = 0.0;
      int u_n = 0;
      for (size_t i = 1; i + 1 < u_face.size(); ++i) { u_sum += u_face[i]; ++u_n; }
      r.u = (u_sum / u_n) / kEb;
      r.d_num = 0.5 * r.u * length / cells;

      // The bed ends at the feed composition, which is the composition the
      // fronts travelled through to within the trace fractions.
      const int mid = cells / 2;
      r.dm_h2 = data.state.molecular_diffusion.at(kH2).var[mid];
      r.dm_co2 = data.state.molecular_diffusion.at(kCO2).var[mid];
      r.d_h2 = data.state.component_axial_dispersion.at(kH2).var[mid];
      r.d_co2 = data.state.component_axial_dispersion.at(kCO2).var[mid];

      auto fit = [&](const std::string& comp, double d_own, double d_other, double& own, double& other) {
        auto [t, y] = Outlet(log, comp);
        std::vector<double> F(y.size());
        for (size_t i = 0; i < y.size(); ++i) F[i] = y[i] / 0.01;
        EXPECT_NEAR(F.back(), 1.0, 0.01) << comp << ": breakthrough incomplete";
        const double shift = CurveMean(t, F) - length / r.u;
        own = ShiftedRms(t, F, length, r.u, d_own + r.d_num, shift, 0.3 * t_front, 2.2 * t_front);
        other = ShiftedRms(t, F, length, r.u, d_other + r.d_num, shift, 0.3 * t_front, 2.2 * t_front);
      };
      fit(kH2, r.d_h2, r.d_co2, r.own_h2, r.other_h2);
      fit(kCO2, r.d_co2, r.d_h2, r.own_co2, r.other_co2);
      std::printf("[trace, %s] D_ax H2 %.3e, CO2 %.3e, numerical %.3e m2/s; RMS H2 %.4f (at CO2's D %.4f), "
                  "CO2 %.4f (at H2's D %.4f)\n", AxialDispersionModeList[static_cast<int>(kind)], r.d_h2, r.d_co2,
                  r.d_num, r.own_h2, r.other_h2, r.own_co2, r.other_co2);
      return r;
    }

  } // namespace

  // 3. Trace limit. The correction to a tracer's flux is y_i S_f, second order
  //    in the trace fraction (S_f is itself of the order of the trace
  //    gradients), so each tracer disperses with its own D_ax,i: its outlet
  //    curve is the Ogata-Banks curve at D_ax,i + u dx / 2 (upwind).
  //
  // GUARD: D_ax,H2 is more than twice D_ax,CO2, and the CO2 curve fits its own
  // coefficient much better than H2's (and the other way round), so a run
  // with one shared coefficient would fail. The molecular diffusivities are
  // checked against measured binary values at 298 K, 1 atm (Marrero and
  // Mason, J. Phys. Chem. Ref. Data 1 (1972) 3: H2-N2 7.79e-5, CO2-N2
  // 1.65e-5 m2/s), scaled to 1 bar.
  TEST(PerComponentDispersion, TraceComponentsSpreadWithTheirOwnCoefficient) {
    UserData data;
    const TraceResult r = RunTrace(AxialDispersionType::RUTHVEN, kU, kLength, kCells, data);
    const double scale_1bar = 101325.0 / kP0;
    EXPECT_NEAR(r.dm_h2, 7.79e-5 * scale_1bar, 0.1 * 7.79e-5 * scale_1bar);
    EXPECT_NEAR(r.dm_co2, 1.65e-5 * scale_1bar, 0.1 * 1.65e-5 * scale_1bar);
    EXPECT_NEAR(r.d_h2, 0.7 * r.dm_h2 + 0.5 * kDp * r.u, 1e-3 * r.d_h2);
    ASSERT_GT(r.d_h2, 2.0 * r.d_co2) << "the tracers disperse alike";
    ASSERT_GT(r.d_co2, 3.0 * r.d_num) << "numerical dispersion dominates";
    EXPECT_LT(r.own_h2, 0.01);
    EXPECT_LT(r.own_co2, 0.01);
    EXPECT_GT(r.other_h2, 4.0 * r.own_h2);
    EXPECT_GT(r.other_co2, 4.0 * r.own_co2);
  }

  // 4. Comparable fractions: the correction is first order here. The mole
  //    fractions still sum to 1 (to the solver tolerance), and the moles of
  //    each component the logged streams carry in and out equal the change of
  //    the bed's inventory.
  TEST(PerComponentDispersion, ComparableFractionsKeepEveryComponentBalance) {
    UserData data;
    const double t_front = kLength / kU;
    BuildBreakthrough(data, {{kN2, 0.3}, {kH2, 0.4}, {kCO2, 0.3}}, 1.5 * t_front, 60);

    CoupledSolver solver{data};
    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;
    ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";

    const auto& log = solver.GetLog();
    const auto& blog = solver.GetBoundaryLog();
    double worst_sum = 0.0;
    const auto y_n2 = log.GetComponentVariable(0, 0, ComponentVariableType::MOLE_FRACTION, kN2);
    const auto y_h2 = log.GetComponentVariable(0, 0, ComponentVariableType::MOLE_FRACTION, kH2);
    const auto y_co2 = log.GetComponentVariable(0, 0, ComponentVariableType::MOLE_FRACTION, kCO2);
    for (size_t s = 0; s < y_n2.values.size(); ++s)
      for (size_t i = 1; i + 1 < y_n2.values[s].size(); ++i)
        worst_sum = std::max(worst_sum, std::abs(y_n2.values[s][i] + y_h2.values[s][i] + y_co2.values[s][i] - 1.0));

    const Step& step = data.cycle.step_list_.front();
    const std::string feed_name = static_cast<const FlowBoundaryCondition*>(step.bottom_boundary_condition_.get())
                                    ->streams_.front()->name_;
    const std::string out_name = step.step_name_ + " Pressure Boundary " + BoundaryLocationList[(int)BoundaryLocation::TOP];
    for (const auto& comp : {kN2, kH2, kCO2}) {
      const auto in = blog.GetComponentVariableData(0, 0, feed_name, BoundaryComponentVariableType::MOLAR_FLOW, comp);
      const auto out = blog.GetComponentVariableData(0, 0, out_name, BoundaryComponentVariableType::MOLAR_FLOW, comp);
      const double moles_in = IntegrateVariableInTime(in.variable_timeseries, in.time);
      const double moles_out = IntegrateVariableInTime(out.variable_timeseries, out.time);  // velocity sign: + leaves the top
      const double inventory =
        CalculateComponentAccumulation_Volume(log, 0, 0, ComponentVariableType::COMPONENT_MOLAR_CONCENTRATION, comp);
      const double scale = std::max({std::abs(moles_in), std::abs(moles_out), std::abs(inventory)});
      std::printf("[comparable] %s: in %.6e, out %.6e, inventory change %.6e mol, closure %.2e\n", comp.c_str(),
                  moles_in, moles_out, inventory, std::abs(moles_in - moles_out - inventory) / scale);
      // Measured 6e-8 to 1.7e-7 (the trapezoidal quadrature of the stream logs).
      EXPECT_LT(std::abs(moles_in - moles_out - inventory), 1e-5 * scale) << comp;
    }
    std::printf("[comparable] max |sum y - 1| %.2e\n", worst_sum);
    EXPECT_LT(worst_sum, 1e-6);
  }

  // 5. The Ruthven and Taylor-Aris models have no molecular-diffusivity input
  //    (the registry lists none for the layer), and a deck saved when they had
  //    one loads: the stored value is dropped.
  TEST(PerComponentDispersion, RegistryHasNoDiffusivityAndOldDecksLoad) {
    UserData data;
    Column c;
    // One layer: cereal reads a class version once per archive and applies
    // it to every object of the class, so the spliced version-2 object must
    // be the only Ruthven model.
    c.layers = {{0.3, 10, 0.4}};
    BuildColumn(data, c);
    const auto first = data.bed.layers_.layers_.front().id_;
    ParameterRegistry reg(data);
    const std::string lp = "layer[" + std::to_string(first) + "].";
    EXPECT_EQ(reg.Find(lp + "molecular_diffusion_coefficient"), nullptr);
    EXPECT_EQ(reg.Find(lp + "component_diffusivity"), nullptr);
    auto* ruthven = dynamic_cast<DispersionModelRuthven*>(data.bed.layers_.layers_.front().equations_.dispersion_model_.get());
    ASSERT_NE(ruthven, nullptr);
    EXPECT_TRUE(ruthven->PerComponent());
    EXPECT_EQ(ruthven->ReadsTransportProperties(), static_cast<unsigned>(kMolecularDiffusion));

    // A version-2 deck: the class version and the two dropped fields.
    std::string json;
    SaveToString(data, json);
    const std::string v3 = "\"polymorphic_name\": \"adsim::DispersionModelRuthven\"";
    size_t at = json.find(v3);
    ASSERT_NE(at, std::string::npos);
    const size_t version_at = json.find("\"cereal_class_version\": 3", at);
    ASSERT_NE(version_at, std::string::npos);
    std::string old = json;
    old.replace(version_at, std::string("\"cereal_class_version\": 3").size(),
                "\"cereal_class_version\": 2, \"Molecular_Diffusion\": {\"cereal_class_version\": 2, \"Value\": 2.5e-05, \"Expression\": "
                "\"2.5e-05\", \"DisplayUnitIndex\": 0}, \"Component_Diffusivity\": 0");
    // (cereal looks the fields up by name, so their place in the object does
    // not matter. The ExpressionParam carries its class version, as the first
    // one an archive loads must.)
    UserData loaded;
    std::string error;
    try {
      LoadFromString(loaded, old);
    } catch (const std::exception& e) {
      error = e.what();
    }
    ASSERT_TRUE(error.empty()) << error;
    EXPECT_NE(dynamic_cast<DispersionModelRuthven*>(loaded.bed.layers_.layers_.front().equations_.dispersion_model_.get()),
              nullptr);
  }

  // =========================================================================
  // Packed-bed correlations: Edwards-Richardson (1968), Wakao-Funazkri (1978)
  // =========================================================================
  //
  // With Pe = u_i d_p / D_ax (u_i interstitial) and Re Sc = u d_p / D_m (u
  // superficial), so eb / (Re Sc) = D_m / (u_i d_p):
  //   Edwards-Richardson: 1/Pe = 0.73 eb/(Re Sc) + 0.5 / (1 + 13 * 0.73 eb/(Re Sc))
  //   Wakao-Funazkri:     1/Pe = 20 / (Re Sc) + 1/2   (eb D_ax / D_m = 20 + 0.5 Sc Re)
  // as Knox et al., Ind. Eng. Chem. Res. 55 (2016) 4734, eqs 10a, 10b and 11
  // print them. The references below are written from these forms.

  namespace {
    double EdwardsRichardsonReference(double dm, double ui, double dp) {
      const double inv_pe_sc = dm / (ui * dp);  // eb / (Re Sc)
      return ui * dp * (0.73 * inv_pe_sc + 0.5 / (1.0 + 13.0 * 0.73 * inv_pe_sc));
    }
    double WakaoFunazkriReference(double dm, double eb, double ui, double dp) {
      const double re_sc = eb * ui * dp / dm;  // superficial u d_p / D_m
      return ui * dp * (20.0 / re_sc + 0.5);
    }
  } // namespace

  // Limits: Edwards-Richardson gives 0.73 D_m at u -> 0 and Pe = 2 at u -> inf;
  // Wakao-Funazkri gives 20 D_m / eb at u -> 0 and Pe = 2 at high Re. Hand
  // values at intermediate points from the dimensionless forms above.
  TEST(DispersionCorrelations, LimitsAndHandValues) {
    const double dm = 1.6e-5, dp = 2.0e-3, eb = 0.4;
    EXPECT_DOUBLE_EQ(DispersionModelEdwardsRichardson::Coefficient(dm, 0.0, dp), 0.73 * dm);
    EXPECT_DOUBLE_EQ(DispersionModelWakaoFunazkri::Coefficient(dm, eb, 0.0, dp), 20.0 * dm / eb);
    for (double ui : {1e4, 1e5}) {  // Pe -> 2 as 20 D_m / (eb u_i d_p) -> 0
      EXPECT_NEAR(ui * dp / DispersionModelEdwardsRichardson::Coefficient(dm, ui, dp), 2.0, 1e-3) << ui;
      EXPECT_NEAR(ui * dp / DispersionModelWakaoFunazkri::Coefficient(dm, eb, ui, dp), 2.0, 1e-3) << ui;
    }
    for (double ui : {1e-4, 3e-3, 0.05, 0.4}) {
      EXPECT_NEAR(DispersionModelEdwardsRichardson::Coefficient(dm, ui, dp), EdwardsRichardsonReference(dm, ui, dp),
                  1e-13 * EdwardsRichardsonReference(dm, ui, dp)) << ui;
      EXPECT_NEAR(DispersionModelWakaoFunazkri::Coefficient(dm, eb, ui, dp), WakaoFunazkriReference(dm, eb, ui, dp),
                  1e-13 * WakaoFunazkriReference(dm, eb, ui, dp)) << ui;
    }
    // A hand calculation: u_i 0.05 m/s, d_p 2 mm: u_i d_p = 1e-4 m2/s,
    // 0.73 * 1.6e-5 = 1.168e-5; 0.5e-4 / (1 + 9.49 * 0.16) = 1.98541e-5.
    EXPECT_NEAR(DispersionModelEdwardsRichardson::Coefficient(dm, 0.05, dp), 3.15341e-5, 1e-9);
    // 20 * 1.6e-5 / 0.4 + 0.5e-4 = 8.5e-4.
    EXPECT_NEAR(DispersionModelWakaoFunazkri::Coefficient(dm, eb, 0.05, dp), 8.5e-4, 1e-12);
  }

  // Metamorphic: the coefficient grows with the velocity, and the
  // Wakao-Funazkri value lies above Edwards-Richardson's everywhere.
  TEST(DispersionCorrelations, GrowWithVelocity) {
    const double dm = 1.6e-5, dp = 2.0e-3, eb = 0.4;
    double prev_er = 0.0, prev_wf = 0.0;
    for (double ui = 1e-5; ui < 10.0; ui *= 1.5) {
      const double er = DispersionModelEdwardsRichardson::Coefficient(dm, ui, dp);
      const double wf = DispersionModelWakaoFunazkri::Coefficient(dm, eb, ui, dp);
      EXPECT_GT(er, prev_er) << ui;
      EXPECT_GT(wf, prev_wf) << ui;
      EXPECT_GT(wf, er) << ui;
      prev_er = er;
      prev_wf = wf;
    }
  }

  // Solved, trace limit: each tracer's outlet curve is the Ogata-Banks curve
  // of its own D_ax,i (+ u dx / 2), and the state's coefficient is the
  // correlation's at the state's D_m,i. Edwards-Richardson at u_i 5 mm/s
  // (molecular regime, D_ax close to 0.73 D_m); Wakao-Funazkri at u_i 0.1 m/s
  // over 1 m (its 20 D_m / eb term would otherwise mix the column).
  TEST(DispersionCorrelations, EdwardsRichardsonTraceBreakthrough) {
    UserData data;
    const TraceResult r = RunTrace(AxialDispersionType::EDWARDS_RICHARDSON, kU, kLength, kCells, data);
    EXPECT_NEAR(r.d_h2, EdwardsRichardsonReference(r.dm_h2, r.u, kDp), 1e-3 * r.d_h2);
    EXPECT_NEAR(r.d_co2, EdwardsRichardsonReference(r.dm_co2, r.u, kDp), 1e-3 * r.d_co2);
    ASSERT_GT(r.d_h2, 2.0 * r.d_co2) << "the tracers disperse alike";
    ASSERT_GT(r.d_co2, 3.0 * r.d_num) << "numerical dispersion dominates";
    EXPECT_LT(r.own_h2, 0.01);
    EXPECT_LT(r.own_co2, 0.01);
    EXPECT_GT(r.other_h2, 4.0 * r.own_h2);
    EXPECT_GT(r.other_co2, 4.0 * r.own_co2);
  }

  TEST(DispersionCorrelations, WakaoFunazkriTraceBreakthrough) {
    UserData data;
    const double u_i = 0.1, length = 1.0;
    const TraceResult r = RunTrace(AxialDispersionType::WAKAO_FUNAZKRI, u_i, length, kCells, data);
    EXPECT_NEAR(r.d_h2, WakaoFunazkriReference(r.dm_h2, kEb, r.u, kDp), 1e-3 * r.d_h2);
    EXPECT_NEAR(r.d_co2, WakaoFunazkriReference(r.dm_co2, kEb, r.u, kDp), 1e-3 * r.d_co2);
    ASSERT_GT(r.d_h2, 2.0 * r.d_co2) << "the tracers disperse alike";
    ASSERT_GT(r.d_co2, 3.0 * r.d_num) << "numerical dispersion dominates";
    const double pe_h2 = r.u * length / r.d_h2;
    ASSERT_GT(pe_h2, 5.0) << "the column is mixed: Ogata-Banks does not apply";
    EXPECT_LT(r.own_h2, 0.01);
    EXPECT_LT(r.own_co2, 0.01);
    EXPECT_GT(r.other_h2, 4.0 * r.own_h2);
    EXPECT_GT(r.other_co2, 4.0 * r.own_co2);
  }

  // The deck checks refuse the packed-bed correlations on a structured layer.
  TEST(DispersionCorrelations, RefusedOnAStructuredLayer) {
    for (AxialDispersionType kind : {AxialDispersionType::EDWARDS_RICHARDSON, AxialDispersionType::WAKAO_FUNAZKRI}) {
      UserData data;
      Column c;
      c.layers = {{0.3, 10, 0.4}};
      c.kind = kind;
      BuildColumn(data, c);
      std::string error;
      EXPECT_TRUE(data.bed.PerformChecks(data.adsorbent_directory, error)) << error;  // fine on a packed layer
      auto& layer = data.bed.layers_.layers_.front();
      layer.kind_ = LayerKind::Laminate;
      layer.equations_.dispersion_model_ = CreateAxialDispersionModel(kind);
      error.clear();
      EXPECT_FALSE(data.bed.PerformChecks(data.adsorbent_directory, error)) << AxialDispersionModeList[static_cast<int>(kind)];
      EXPECT_NE(error.find("packed-bed dispersion"), std::string::npos) << error;
    }
  }

} // namespace adsim
