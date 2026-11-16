/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Axial dispersion with a coefficient per component (component_diffusivity =
// mixture_averaged on the Ruthven and Taylor-Aris models).
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
//  2. Residual level, binary: with two components D_m,1 = D_m,2 = D_12, so the
//     mixture-averaged residuals equal the shared-mode residuals with
//     D_m = D_12 (the correction is then zero to round-off).
//  3. Solved, trace limit: H2 and CO2 at 1 % each in N2 through an inert bed.
//     The correction is second order in the trace fraction, so each tracer's
//     outlet curve is the Ogata-Banks curve of its own D_ax,i; H2 spreads more
//     than CO2.
//  4. Solved, comparable fractions: a 40 % H2, 30 % CO2 feed into N2. The
//     mole fractions keep summing to 1 and the logged streams account for
//     every component's inventory change.
//  5. Registry and deck: the choice is a registry entry of the Ruthven and
//     Taylor-Aris models only, and it survives a save and a load.

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
      ComponentDiffusivity mode = ComponentDiffusivity::MIXTURE_AVERAGED;
      double shared_dm = 1.6e-5;  // [m2/s], the shared mode's D_m
      double particle_diameter = 0.002;
    };

    // An isothermal column of inert packing (no isotherm), each layer with
    // Ruthven dispersion in the given mode.
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
        auto disp = std::make_unique<DispersionModelRuthven>();
        disp->Dm = c.shared_dm;
        disp->component_diffusivity_ = c.mode;
        layer.equations_.dispersion_model_ = std::move(disp);
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

  // 2. Binary: Blanc's law gives D_m,1 = D_m,2 = D_12, so the mixture-averaged
  //    residuals are the shared-mode residuals with D_m = D_12. Uniform T and P
  //    make D_12 uniform; it is read from the mixture-averaged run.
  TEST(PerComponentDispersion, BinaryMixtureAveragedEqualsSharedAtD12) {
    Column mixed;
    mixed.components = {kN2, kCO2};
    const StillColumn a = StillColumnResiduals(mixed);
    const int mid = a.interface_cell + 2;
    // D_ax = 0.7 D_m in a still column.
    const double d12 = a.dax[0][mid] / 0.7;
    EXPECT_NEAR(a.dax[0][mid], a.dax[1][mid], 1e-15 * a.dax[0][mid]);

    Column shared = mixed;
    shared.mode = ComponentDiffusivity::SHARED;
    shared.shared_dm = d12;
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
    void BuildBreakthrough(UserData& data, const std::vector<std::pair<std::string, double>>& feed,
                           double run_time, int cells) {
      Column c;
      c.layers = {{kLength, cells, kEb}};
      c.particle_diameter = kDp;
      BuildColumn(data, c);
      data.settings.advection_scheme_ = AdvectionScheme::Upwind;

      const double area = 0.25 * std::acos(-1.0) * kDiameter * kDiameter;
      const double molar_flow = kEb * kU * area * kP0 / (kR * kT0);

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
    double ShiftedRms(const std::vector<double>& t, const std::vector<double>& F, double u, double D,
                      double t_shift, double t_lo, double t_hi) {
      double sum = 0.0;
      int n = 0;
      for (size_t i = 0; i < t.size(); ++i) {
        if (t[i] < t_lo || t[i] > t_hi) continue;
        const double d = F[i] - OgataBanks(kLength, u, D, t[i] - t_shift);
        sum += d * d;
        ++n;
      }
      return n > 0 ? std::sqrt(sum / n) : 1.0;
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
    const double t_front = kLength / kU;
    BuildBreakthrough(data, {{kN2, 0.98}, {kH2, 0.01}, {kCO2, 0.01}}, 2.5 * t_front, kCells);

    CoupledSolver solver{data};
    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;
    ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";
    ASSERT_TRUE(data.state.PerComponentDispersion());

    const auto& log = solver.GetLog();
    const auto vel = log.GetVariable(0, 0, VariableType::VELOCITY);
    const auto& u_face = vel.values.back();
    double u_sum = 0.0;
    int u_n = 0;
    for (size_t i = 1; i + 1 < u_face.size(); ++i) { u_sum += u_face[i]; ++u_n; }
    const double u = (u_sum / u_n) / kEb;
    const double dx = kLength / kCells;
    const double d_num = 0.5 * u * dx;

    // The bed ends at the feed composition, which is the composition the
    // fronts travelled through to within the trace fractions.
    const int mid = kCells / 2;
    const double scale_1bar = 101325.0 / kP0;
    const double dm_h2 = data.state.molecular_diffusion.at(kH2).var[mid];
    const double dm_co2 = data.state.molecular_diffusion.at(kCO2).var[mid];
    EXPECT_NEAR(dm_h2, 7.79e-5 * scale_1bar, 0.1 * 7.79e-5 * scale_1bar);
    EXPECT_NEAR(dm_co2, 1.65e-5 * scale_1bar, 0.1 * 1.65e-5 * scale_1bar);
    const double d_h2 = data.state.component_axial_dispersion.at(kH2).var[mid];
    const double d_co2 = data.state.component_axial_dispersion.at(kCO2).var[mid];
    EXPECT_NEAR(d_h2, 0.7 * dm_h2 + 0.5 * kDp * u, 1e-3 * d_h2);
    ASSERT_GT(d_h2, 2.0 * d_co2) << "the tracers disperse alike";
    ASSERT_GT(d_co2, 3.0 * d_num) << "numerical dispersion dominates";

    struct Fit { double own, other; };
    auto fit = [&](const std::string& comp, double d_own, double d_other) {
      auto [t, y] = Outlet(log, comp);
      std::vector<double> F(y.size());
      for (size_t i = 0; i < y.size(); ++i) F[i] = y[i] / 0.01;
      EXPECT_NEAR(F.back(), 1.0, 0.01) << comp << ": breakthrough incomplete";
      const double shift = CurveMean(t, F) - kLength / u;
      return Fit{ShiftedRms(t, F, u, d_own + d_num, shift, 0.3 * t_front, 2.2 * t_front),
                 ShiftedRms(t, F, u, d_other + d_num, shift, 0.3 * t_front, 2.2 * t_front)};
    };
    const Fit h2 = fit(kH2, d_h2, d_co2);
    const Fit co2 = fit(kCO2, d_co2, d_h2);
    std::printf("[trace] D_ax H2 %.3e, CO2 %.3e, numerical %.3e m2/s; RMS H2 %.4f (at CO2's D %.4f), "
                "CO2 %.4f (at H2's D %.4f)\n", d_h2, d_co2, d_num, h2.own, h2.other, co2.own, co2.other);
    EXPECT_LT(h2.own, 0.01);
    EXPECT_LT(co2.own, 0.01);
    EXPECT_GT(h2.other, 4.0 * h2.own);
    EXPECT_GT(co2.other, 4.0 * co2.own);
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

  // 5. The choice is a registry entry of the Ruthven and Taylor-Aris models
  //    (not of the Constant model), and a saved deck keeps it.
  TEST(PerComponentDispersion, RegistryEntryAndDeckRoundTrip) {
    UserData data;
    Column c;
    c.mode = ComponentDiffusivity::SHARED;
    BuildColumn(data, c);
    auto constant = std::make_unique<DispersionModelConstant>();
    data.bed.layers_.layers_.back().equations_.dispersion_model_ = std::move(constant);
    const auto first = data.bed.layers_.layers_.front().id_;
    const auto second = data.bed.layers_.layers_.back().id_;

    ParameterRegistry reg(data);
    const std::string path = "layer[" + std::to_string(first) + "].component_diffusivity";
    Parameter* p = reg.Find(path);
    ASSERT_NE(p, nullptr);
    EXPECT_EQ(p->Info().choices, (std::vector<std::string>{"shared", "mixture_averaged"}));
    EXPECT_EQ(p->GetValue(), 0.0);
    EXPECT_EQ(reg.Find("layer[" + std::to_string(second) + "].component_diffusivity"), nullptr);
    p->SetText("mixture_averaged");
    auto* ruthven = dynamic_cast<DispersionModelRuthven*>(data.bed.layers_.layers_.front().equations_.dispersion_model_.get());
    ASSERT_NE(ruthven, nullptr);
    EXPECT_EQ(ruthven->component_diffusivity_, ComponentDiffusivity::MIXTURE_AVERAGED);
    EXPECT_TRUE(ruthven->PerComponent());
    EXPECT_EQ(ruthven->ReadsTransportProperties(), static_cast<unsigned>(kMolecularDiffusion));

    std::string json;
    SaveToString(data, json);
    UserData loaded;
    LoadFromString(loaded, json);
    auto* back = dynamic_cast<DispersionModelRuthven*>(loaded.bed.layers_.layers_.front().equations_.dispersion_model_.get());
    ASSERT_NE(back, nullptr);
    EXPECT_EQ(back->component_diffusivity_, ComponentDiffusivity::MIXTURE_AVERAGED);
  }

} // namespace adsim
