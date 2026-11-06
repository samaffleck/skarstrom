/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Solved tests for structured layers (MixedAdsorbentLayer with a monolith or
// laminate kind): laminar flow through straight channels.
//
// In a channel of hydraulic diameter D_h, fully developed laminar flow obeys
// Hagen-Poiseuille,  -dP/dz = 2 (f Re) mu u_ch / D_h^2,  with the Fanning
// friction factor times Reynolds number of the cross-section (Shah and
// London, Laminar Flow Forced Convection in Ducts, 1978): 14.23 for a square,
// 16 for a circle, 24 for parallel plates. The superficial velocity is
// u = OFA u_ch. For an isothermal ideal gas with the mass flux G = rho u
// constant and rho = P M / (R T) this integrates to
//
//   P_in^2 - P_out^2 = 2 L (R T / M) (2 fRe / (OFA D_h^2)) mu G.
//
// The geometry is re-derived here from its definition (pitch p = n^-1/2,
// open width a = p - w - 2 t_wc; laminate open gap g - 2 t_wc and period
// g + t), not taken from the solver.

#include "breakthrough_case.h"
#include "adsim/solver/operating_limits.h"

#include <cmath>
#include <functional>
#include <string>

namespace adsim {

  using namespace adsim::validation;

  namespace {

    constexpr double kLength = 0.2;       // [m]
    constexpr double kDiameter = 0.05;    // [m]
    constexpr double kT = 300.0;          // [K]
    constexpr double kPout = 1.0e5;       // [Pa]
    constexpr double kMolarFlow = 2.0e-2; // [mol/s]: laminar (Re ~ 25) in every channel below

    struct Channel {
      std::string name;
      LayerKind kind;
      ChannelShape shape;
      double cell_density, wall, sheet, gap, washcoat;
    };

    // Open fraction, hydraulic diameter and fRe from the definitions.
    void Geometry(const Channel& c, double& ofa, double& dh, double& fre) {
      if (c.kind == LayerKind::Laminate) {
        const double open = c.gap - 2.0 * c.washcoat;
        ofa = open / (c.gap + c.sheet);
        dh = 2.0 * open;  // 4 (open x width) / (2 width)
        fre = 24.0;
        return;
      }
      const double p = 1.0 / std::sqrt(c.cell_density);
      const double a = p - c.wall - 2.0 * c.washcoat;
      dh = a;  // square: 4 a^2 / (4 a); circle: its diameter
      if (c.shape == ChannelShape::Circular) {
        ofa = M_PI * a * a / (4.0 * p * p);
        fre = 16.0;
      } else {
        ofa = a * a / (p * p);
        fre = 14.23;
      }
    }

    struct Result {
      double dP = 0.0, mu = 0.0, M = 0.0, u_mid = 0.0;
    };

    using Inspect = std::function<void(const UserData&, const CoupledSolver&)>;

    Result RunChannelFlow(const Channel& c, int n_cells = 40, double molar_flow = kMolarFlow,
                          const Inspect& inspect = nullptr) {
      UserData data;
      data.fluid.AddComponent("Nitrogen");
      const auto ads = data.adsorbent_directory.Add("Inert");
      data.adsorbent_directory.Get(ads).ep_ = 0.0;

      const auto lid = data.bed.layers_.Add("Structured");
      auto& layer = data.bed.layers_.Get(lid);
      layer.domain_.length = kLength;
      layer.domain_.n_cells = n_cells;
      layer.composition_[ads] = 1.0;
      layer.kind_ = c.kind;
      layer.structured_.channel_shape_ = c.shape;
      layer.structured_.cell_density_ = c.cell_density;
      layer.structured_.wall_thickness_ = c.wall;
      layer.structured_.sheet_thickness_ = c.sheet;
      layer.structured_.gap_ = c.gap;
      layer.structured_.washcoat_thickness_ = c.washcoat;
      auto disp = std::make_unique<DispersionModelConstant>();
      disp->dispersion_coefficient_ = 0.0;
      layer.equations_.dispersion_model_ = std::move(disp);
      std::get<CylinderBed>(data.bed.shape_).diameter = kDiameter;
      data.bed.equations_.energy_balance_mode_ = EnergyBalanceMode::ISOTHERMAL;

      data.ic.pressure = kPout;
      data.ic.temperature = kT;
      data.ic.wall_temperature = kT;
      data.ic.velocity = 0.0;
      data.ic.mole_fraction["Nitrogen"] = 1.0;

      data.cycle.AddStep("Flow");
      Step& step = data.cycle.GetStep("Flow");
      step.step_time_ = 10.0;
      auto inlet = std::make_unique<BoundaryStream>();
      inlet->name_ = "Feed";
      auto flow = std::make_unique<MolarFlowBoundaryValue>();
      static_cast<ConstantBoundaryValue*>(flow->profile_.get())->constant_value_ = molar_flow;
      inlet->boundary_flow_ = std::move(flow);
      inlet->mole_fraction["Nitrogen"] = 1.0;
      inlet->temperature_ = kT;
      auto inlet_bc = std::make_unique<FlowBoundaryCondition>();
      inlet_bc->direction_ = BoundaryDirection::INLET;
      inlet_bc->location_ = BoundaryLocation::BOTTOM;
      inlet_bc->streams_.push_back(std::move(inlet));
      step.bottom_boundary_condition_ = std::move(inlet_bc);
      auto outlet_bc = std::make_unique<PressureBoundaryCondition>();
      outlet_bc->direction_ = BoundaryDirection::OUTLET;
      outlet_bc->location_ = BoundaryLocation::TOP;
      outlet_bc->pressure_ = std::make_unique<ConstantBoundaryValue>(kPout);
      outlet_bc->mole_fraction_["Nitrogen"] = 1.0;
      step.top_boundary_condition_ = std::move(outlet_bc);
      data.cycle.number_of_cycles_ = 1;
      data.settings.reltol = 1e-10;
      data.settings.max_step = 1.0;
      data.settings.advection_scheme_ = AdvectionScheme::Upwind;

      std::string err;
      EXPECT_TRUE(data.EvaluateAllExpressions(err) && data.PerformChecks(err)) << c.name << ": " << err;
      CoupledSolver solver{data};
      EXPECT_GE(solver.IntegrateCycle(data.cycle, false), 0) << c.name << ": solver did not complete";

      if (inspect) inspect(data, solver);
      const auto& log = solver.GetLog();
      const auto P = log.GetVariable(0, 0, VariableType::PRESSURE).values.back();
      const auto mu = log.GetVariable(0, 0, VariableType::VISCOSITY).values.back();
      const auto u = log.GetVariable(0, 0, VariableType::VELOCITY).values.back();
      const size_t N = P.size();
      Result r;
      // Face pressures at z = 0 and z = L (ghost / end-cell averages).
      r.dP = 0.5 * (P[0] + P[1]) - 0.5 * (P[N - 2] + P[N - 1]);
      r.mu = mu[N / 2];
      r.M = data.fluid.molar_mass_.at("Nitrogen");
      r.u_mid = u[u.size() / 2];
      return r;
    }

  } // namespace

  // Each shape at a typical size (400 cpsi monoliths, one with a washcoat; a
  // 0.2 mm laminate): the simulated pressure drop is the closed form.
  // TOLERANCE 1e-5: the steady state is reached (10 s against a gas residence
  // time of ~0.1 s) and the momentum balance is linear, so only the
  // integrator (reltol 1e-10) and the viscosity's variation with the
  // pressure (none for the ideal-gas viscosity model) remain.
  TEST(StructuredLayers, HagenPoiseuillePressureDropForEachShape) {
    const Channel cases[] = {
      {"square monolith", LayerKind::Monolith, ChannelShape::Square, 6.2e5, 1.65e-4, 0.0, 0.0, 0.0},
      {"circular monolith", LayerKind::Monolith, ChannelShape::Circular, 6.2e5, 1.65e-4, 0.0, 0.0, 0.0},
      {"square monolith with washcoat", LayerKind::Monolith, ChannelShape::Square, 6.2e5, 1.0e-4, 0.0, 0.0, 5.0e-5},
      {"laminate", LayerKind::Laminate, ChannelShape::Square, 0.0, 0.0, 2.0e-4, 2.0e-4, 0.0},
      {"laminate with washcoat", LayerKind::Laminate, ChannelShape::Square, 0.0, 0.0, 2.0e-4, 3.0e-4, 5.0e-5},
    };
    const double area = 0.25 * M_PI * kDiameter * kDiameter;
    for (const auto& c : cases) {
      double ofa, dh, fre;
      Geometry(c, ofa, dh, fre);
      const Result r = RunChannelFlow(c);
      const double G = kMolarFlow * r.M / area;  // [kg/(m2 s)]
      const double rhs = 2.0 * kLength * (kR * kT / r.M) * (2.0 * fre / (ofa * dh * dh)) * r.mu * G;
      const double dP = std::sqrt(kPout * kPout + rhs) - kPout;
      // Channel Reynolds number rho u_ch D_h / mu = G D_h / (OFA mu): laminar.
      const double re = G * dh / (ofa * r.mu);
      std::cout << "[" << c.name << "] OFA " << ofa << ", D_h " << dh << " m, Re " << re << ", dP simulated "
                << r.dP << " Pa, analytic " << dP << " Pa\n";
      ASSERT_LT(re, 2000.0) << c.name;
      ASSERT_GT(dP, 10.0) << c.name << ": pressure drop too small to test";
      EXPECT_NEAR(r.dP, dP, 1e-5 * dP) << c.name;
    }
  }

  // The three shapes at the same open width compare as their fRe: at equal
  // OFA D_h^2 / fRe the drops agree, so a wrong fRe in any one shows. Here,
  // the circular monolith at the square's pitch has OFA pi/4 times smaller and
  // fRe 16 / 14.23 larger, so its P_in^2 - P_out^2 is (4/pi)(16/14.23) times
  // the square's.
  TEST(StructuredLayers, CircularOverSquareIsTheRatioOfOpenAreaAndFriction) {
    const Channel square{"square", LayerKind::Monolith, ChannelShape::Square, 6.2e5, 1.65e-4, 0.0, 0.0, 0.0};
    Channel circle = square;
    circle.name = "circle";
    circle.shape = ChannelShape::Circular;
    const Result rs = RunChannelFlow(square);
    const Result rc = RunChannelFlow(circle);
    auto sq = [](const Result& r) { return (kPout + r.dP) * (kPout + r.dP) - kPout * kPout; };
    ASSERT_NEAR(rs.mu, rc.mu, 1e-12 * rs.mu);
    EXPECT_NEAR(sq(rc) / sq(rs), (4.0 / M_PI) * (16.0 / 14.23), 1e-5);
  }

  // A structured layer is a fixed solid: an upward pressure gradient above
  // the weight of the solid does not lift it, so the operating limits give no
  // fluidisation warning for its faces. Guard: the logged approach to
  // fluidisation (computed as for a packed bed) exceeds 1.
  TEST(StructuredLayers, NoFluidisationWarningInsideAStructuredLayer) {
    const Channel thin{"thin laminate", LayerKind::Laminate, ChannelShape::Square, 0.0, 0.0, 2.0e-4, 5.0e-5, 0.0};
    double max_approach = 0.0;
    size_t fluidisation_warnings = 0;
    RunChannelFlow(thin, 40, 2.0e-2, [&](const UserData& data, const CoupledSolver& solver) {
      const auto a = solver.GetLog().GetVariable(0, 0, VariableType::APPROACH_TO_FLUIDISATION).values;
      for (const auto& row : a)
        for (size_t i = 1; i + 2 < row.size(); ++i) max_approach = std::max(max_approach, row[i]);
      for (const auto& w : CheckOperatingLimits(data, solver.GetLog()).warnings)
        fluidisation_warnings += w.limit == OperatingLimit::Fluidisation;
    });
    ASSERT_GT(max_approach, 1.0) << "the guard needs a pressure gradient above the solid's weight";
    EXPECT_EQ(fluidisation_warnings, 0u);
  }

} // namespace adsim
