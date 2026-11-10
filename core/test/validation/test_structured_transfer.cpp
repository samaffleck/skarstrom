/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Heat and mass transfer in structured layers (monolith, laminate), each
// closure checked against its source's own result, written out here
// independently of core/src:
//
//  1. Fully developed laminar Nu (constant wall temperature) and Sh by the
//     heat-mass analogy, Shah and London, Laminar Flow Forced Convection in
//     Ducts (Academic Press, 1978): 2.976 square, 3.657 circle, 7.541
//     parallel plates.
//  2. Hawthorn (1974, AIChE Symp. Ser. 70(137)) length-averaged developing
//     flow: Nu = Nu_inf (1 + 0.095 Re Pr D_h / L)^0.45.
//  3. Taylor-Aris dispersion: Taylor (1953, Proc. R. Soc. A 219, 186),
//     K = a^2 U^2 / (48 D) in a tube of radius a; Aris (1956, Proc. R. Soc.
//     A 235, 67), K = W^2 U^2 / (210 D) between plates W apart; both with the
//     molecular term D added.
//  4. Glueckauf (1955, Trans. Faraday Soc. 51, 1540) linear driving force:
//     k = 15 D / r^2 for a sphere of radius r, 3 D / l^2 for a slab of
//     half-thickness l.
//  5. Wakao and Funazkri (1978, Chem. Eng. Sci. 33, 1375), packed beds:
//     Nu = 2 + 1.1 Pr^(1/3) Re^0.6; with Gunn (1978, Int. J. Heat Mass
//     Transfer 21, 467) and Ranz and Marshall (1952, Chem. Eng. Prog. 48,
//     141) as the other packed-layer film correlations (5b).
//  6. A linear-isotherm breakthrough through a parallel-passage contactor
//     (laminate) against the moment analysis of the LDF model (derived
//     below), which exercises the slab LDF, Taylor-Aris dispersion and the
//     open-area holdup together.
//  7. The channel Reynolds number warning of the operating limits.
//
// The equation numbers of these results in their sources are to be added by
// whoever checks the texts; the results themselves are stated above.

#include "breakthrough_case.h"
#include "adsim/equations/heat_transfer_model.h"
#include "adsim/solver/operating_limits.h"

#include <cmath>
#include <functional>
#include <string>

namespace adsim {

  using namespace adsim::validation;

  namespace {

    constexpr double kT = 300.0;     // [K]
    constexpr double kPout = 1.0e5;  // [Pa]
    const std::string kCarrier = "Nitrogen";
    const std::string kTracer = "CarbonDioxide";

    // A one-layer structured bed with an N2 carrier and a dilute CO2 tracer
    // on a Henry isotherm with the macropore LDF (constant D_m), Taylor-Aris
    // dispersion and the film gas-solid model, isothermal.
    struct ChannelDeck {
      LayerKind kind = LayerKind::Laminate;
      ChannelShape shape = ChannelShape::Square;
      double cell_density = 6.2e5, wall = 1.65e-4;  // monolith
      double sheet = 2.0e-3, gap = 2.0e-4;           // laminate
      double washcoat = 0.0;
      double length = 0.2, diameter = 0.05;          // [m]
      int cells = 40;
      double molar_flow = 2.0e-3;                    // [mol/s]
      double y_feed = 1.0e-3;                        // tracer [-]
      double henry_k = 5.0e-8;                       // q = k p [mol/(kg Pa)]
      double rho_s = 1100.0, ep = 0.35, tortuosity = 3.0;
      double dm_ldf = 1.0e-5;                        // [m2/s], the LDF's constant D_m
      double dm_dispersion = 1.6e-5;                 // [m2/s], Taylor-Aris D_m
      bool ltne = false;
      bool packed = false;                           // a packed bed instead (for the packed closures)
      FilmCorrelation film = FilmCorrelation::WakaoFunazkri;  // packed layer's film correlation
      bool haynes = false;                           // Haynes rate coefficient (film + pore) instead
      double run_time = 600.0, reltol = 1e-8, max_step = 0.5;
    };

    void Build(UserData& data, const ChannelDeck& d) {
      data.fluid.AddComponent(kCarrier);
      data.fluid.AddComponent(kTracer);
      const auto id = data.adsorbent_directory.Add("Sheet");
      auto& ads = data.adsorbent_directory.Get(id);
      ads.solid_density_ = d.rho_s;
      ads.ep_ = d.ep;
      ads.heat_capacity_ = 900.0;
      ads.thermal_conductivity_ = 0.5;
      ads.adsorption_source_.Initialise(data.fluid);
      auto henry = std::make_unique<HenryIsotherm>();
      henry->params_[kTracer].k = d.henry_k;
      henry->params_[kTracer].Q = 0.0;
      ads.adsorption_source_.AddIsotherm({kTracer}, std::move(henry));
      auto ldf = std::make_unique<LinearDrivingForce>();
      if (d.haynes) {
        auto haynes = std::make_unique<MassTransferCoefficient_Haynes>();
        haynes->macropore_void_fraction_ = d.ep;
        haynes->tortuosity_ = d.tortuosity;
        ldf->mtc_model_ = ToResistances(std::move(haynes));  // generic model, as the example decks store it
      } else {
        auto macro = std::make_unique<MassTransferCoefficientMacropore>();
        macro->macropore_void_fraction_ = d.ep;
        macro->tortuosity_ = d.tortuosity;
        macro->molecular_diffusion_coefficient_ = d.dm_ldf;
        ldf->mtc_model_ = ToResistances(std::move(macro));  // generic model, as the example decks store it
      }
      ads.adsorption_source_.mass_transfer_models_[kTracer] = std::move(ldf);
      ads.adsorption_source_.heat_of_adsorption_[kTracer] = 0.0;

      const auto lid = data.bed.layers_.Add("Channels");
      auto& layer = data.bed.layers_.Get(lid);
      layer.domain_.length = d.length;
      layer.domain_.n_cells = d.cells;
      layer.composition_[id] = 1.0;
      if (d.packed) {
        layer.kind_ = LayerKind::Packed;
        layer.eb_ = 0.4;
        layer.particle_diameter_ = 2.0e-3;
        layer.momentum_model_ = MomentumModel::BlakeKozeny;
        layer.film_correlation_ = d.film;
        auto disp = std::make_unique<DispersionModelConstant>();
        disp->dispersion_coefficient_ = 1e-5;
        layer.equations_.dispersion_model_ = std::move(disp);
      } else {
        layer.kind_ = d.kind;
        layer.structured_.channel_shape_ = d.shape;
        layer.structured_.cell_density_ = d.cell_density;
        layer.structured_.wall_thickness_ = d.wall;
        layer.structured_.sheet_thickness_ = d.sheet;
        layer.structured_.gap_ = d.gap;
        layer.structured_.washcoat_thickness_ = d.washcoat;
        auto ta = std::make_unique<DispersionModelTaylorAris>();
        ta->Dm = d.dm_dispersion;
        layer.equations_.dispersion_model_ = std::move(ta);
      }
      std::get<CylinderBed>(data.bed.shape_).diameter = d.diameter;
      data.bed.equations_.particle_heat_transfer_model_ =
          CreateParticleHeatTransferModel(ParticleHeatTransferType::FILM);
      if (d.ltne) {
        data.bed.equations_.thermal_equilibrium_mode_ = ThermalEquilibriumMode::NON_LOCAL;
      } else {
        data.bed.equations_.energy_balance_mode_ = EnergyBalanceMode::ISOTHERMAL;
      }

      data.ic.pressure = kPout;
      data.ic.temperature = kT;
      data.ic.wall_temperature = kT;
      data.ic.velocity = 0.0;
      data.ic.mole_fraction[kCarrier] = 1.0;
      data.ic.mole_fraction[kTracer] = 0.0;

      data.cycle.AddStep("Feed");
      Step& step = data.cycle.GetStep("Feed");
      step.step_time_ = d.run_time;
      auto inlet = std::make_unique<BoundaryStream>();
      inlet->name_ = "Feed inlet";
      auto flow = std::make_unique<MolarFlowBoundaryValue>();
      static_cast<ConstantBoundaryValue*>(flow->profile_.get())->constant_value_ = d.molar_flow;
      inlet->boundary_flow_ = std::move(flow);
      inlet->mole_fraction[kCarrier] = 1.0 - d.y_feed;
      inlet->mole_fraction[kTracer] = d.y_feed;
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
      outlet_bc->mole_fraction_[kCarrier] = 1.0 - d.y_feed;
      outlet_bc->mole_fraction_[kTracer] = d.y_feed;
      step.top_boundary_condition_ = std::move(outlet_bc);
      data.cycle.number_of_cycles_ = 1;
      data.settings.reltol = d.reltol;
      data.settings.max_step = d.max_step;
    }

    using Inspect = std::function<void(const UserData&, const CoupledSolver&)>;

    void RunChannelDeck(const ChannelDeck& d, const Inspect& inspect) {
      UserData data;
      Build(data, d);
      std::string err;
      ASSERT_TRUE(data.EvaluateAllExpressions(err) && data.PerformChecks(err)) << err;
      CoupledSolver solver{data};
      ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";
      inspect(data, solver);
    }

    // The value at the middle interior index at the end of the run.
    double Mid(const std::vector<Var1D>& rows) { return rows.back()[rows.back().size() / 2]; }

  } // namespace

  // 1. Shah and London's tabulated Nu_T, and Sh = Nu by the analogy.
  TEST(StructuredTransfer, FullyDevelopedNusseltIsShahAndLondonsTable) {
    MixedAdsorbentLayer l;
    EXPECT_EQ(l.GetChannelNusselt(), 0.0);  // packed: no channel value
    l.kind_ = LayerKind::Monolith;
    l.structured_.channel_shape_ = ChannelShape::Square;
    EXPECT_DOUBLE_EQ(l.GetChannelNusselt(), 2.976);
    l.structured_.channel_shape_ = ChannelShape::Circular;
    EXPECT_DOUBLE_EQ(l.GetChannelNusselt(), 3.657);
    l.kind_ = LayerKind::Laminate;
    EXPECT_DOUBLE_EQ(l.GetChannelNusselt(), 7.541);
  }

  // 2. Hawthorn's developing-flow form: the fully developed value for a long
  //    channel, and the stated power law otherwise.
  TEST(StructuredTransfer, HawthornDevelopingFlowForm) {
    const double nu = 2.976, dh = 1.1e-3;
    EXPECT_DOUBLE_EQ(ChannelFilmNumber(nu, 50.0, 0.7, dh, 0.0), nu);  // switched off
    EXPECT_NEAR(ChannelFilmNumber(nu, 50.0, 0.7, dh, 1e6), nu, 1e-6 * nu);  // L >> Re Pr D_h
    for (double re : {1.0, 30.0, 500.0}) {
      for (double L : {0.05, 0.3}) {
        const double gz = re * 0.7 * dh / L;  // Graetz-type group
        EXPECT_NEAR(ChannelFilmNumber(nu, re, 0.7, dh, L), nu * std::pow(1.0 + 0.095 * gz, 0.45), 1e-13 * nu);
      }
    }
  }

  // 3. Taylor-Aris: the logged axial dispersion coefficient against Taylor's
  //    radius form (circular channels) and Aris's gap form (plates), with the
  //    channel velocity u / OFA from the logged superficial velocity.
  TEST(StructuredTransfer, TaylorArisCoefficientIsTaylorsAndArissResult) {
    for (LayerKind kind : {LayerKind::Monolith, LayerKind::Laminate}) {
      ChannelDeck d;
      d.kind = kind;
      d.shape = ChannelShape::Circular;
      d.gap = 4.0e-4;
      d.sheet = 2.0e-4;
      d.molar_flow = kind == LayerKind::Laminate ? 1.0e-1 : 2.0e-2;
      d.run_time = 5.0;
      d.cells = 20;
      RunChannelDeck(d, [&](const UserData& data, const CoupledSolver& solver) {
        const auto& log = solver.GetLog();
        // The coefficient and the velocity at the same index (the closure
        // reads u[n] for D_ax[n]).
        const auto Dx = log.GetVariable(0, 0, VariableType::AXIAL_DISPERSION_COEFFICIENT).values.back();
        const auto uu = log.GetVariable(0, 0, VariableType::VELOCITY).values.back();
        const size_t n = Dx.size() / 2;
        ASSERT_LT(n, uu.size());
        const double u = std::abs(uu[n]);
        const double D_logged = Dx[n];
        const auto& layer = data.bed.layers_.layers_.front();
        const double U = u / layer.GetOpenFraction();  // mean channel velocity
        const double D = d.dm_dispersion;
        double expected;
        if (kind == LayerKind::Monolith) {
          const double a = 0.5 * layer.GetHydraulicDiameter();  // tube radius
          expected = D + a * a * U * U / (48.0 * D);             // Taylor (1953)
        } else {
          const double W = d.gap;                                // plate spacing
          expected = D + W * W * U * U / (210.0 * D);            // Aris (1956)
        }
        std::cout << "[Taylor-Aris " << LayerKindList[static_cast<int>(kind)] << "] U " << U << " m/s, D_ax "
                  << D_logged << " m2/s, analytic " << expected << "\n";
        ASSERT_GT(expected, 2.0 * D) << "the convective term must matter";
        EXPECT_NEAR(D_logged, expected, 1e-12 * expected);
      });
    }
  }

  // 4. Glueckauf's LDF: 3 D_eff / l^2 for the slab of a laminate sheet (l =
  //    t / 2) and of a washcoat (l = t_wc); 15 D_eff / r_p^2 for packed
  //    spheres, as before. D_eff = eps_p D_m / tau.
  TEST(StructuredTransfer, SlabLdfIsGlueckaufsThreeOverHalfThicknessSquared) {
    struct Case { const char* name; bool packed; double washcoat; double l; double c; };
    const ChannelDeck base;
    const Case cases[] = {
      {"laminate sheet", false, 0.0, 0.5 * base.sheet, 3.0},
      {"laminate washcoat", false, 5.0e-5, 5.0e-5, 3.0},
      {"packed spheres", true, 0.0, 1.0e-3, 15.0},  // r_p = d_p / 2
    };
    for (const auto& c : cases) {
      ChannelDeck d;
      d.packed = c.packed;
      d.washcoat = c.washcoat;
      d.gap = 4.0e-4;
      d.run_time = 2.0;
      d.cells = 10;
      RunChannelDeck(d, [&](const UserData&, const CoupledSolver& solver) {
        const double k = Mid(solver.GetLog().GetComponentVariable(0, 0, ComponentVariableType::MASS_TRANSFER_COEFFICIENT, kTracer).values);
        const double D_eff = d.ep * d.dm_ldf / d.tortuosity;
        const double expected = c.c * D_eff / (c.l * c.l);
        std::cout << "[LDF " << c.name << "] k " << k << " 1/s, Glueckauf " << expected << " 1/s\n";
        EXPECT_NEAR(k, expected, 1e-12 * expected) << c.name;
      });
    }
  }

  // 5. The film gas-solid model: Shah and London's Nu on the channel wall
  //    area for a laminate, Wakao and Funazkri on 6 (1 - eb) / d_p for packed
  //    spheres, from the logged gas properties.
  TEST(StructuredTransfer, FilmGasSolidCoefficient) {
    for (bool packed : {false, true}) {
      ChannelDeck d;
      d.packed = packed;
      d.ltne = true;
      d.gap = 4.0e-4;
      d.molar_flow = 1.0e-2;
      d.run_time = 2.0;
      d.cells = 10;
      RunChannelDeck(d, [&](const UserData& data, const CoupledSolver& solver) {
        const auto& log = solver.GetLog();
        const double hp = Mid(log.GetVariable(0, 0, VariableType::PARTICLE_HEAT_TRANSFER_COEFFICIENT).values);
        const double kg = Mid(log.GetVariable(0, 0, VariableType::THERMAL_CONDUCTIVITY).values);
        const double cp = Mid(log.GetVariable(0, 0, VariableType::HEAT_CAPACITY).values);
        const double mu = Mid(log.GetVariable(0, 0, VariableType::VISCOSITY).values);
        const double re = Mid(log.GetVariable(0, 0, VariableType::REYNOLDS_NUMBER).values);
        const auto& layer = data.bed.layers_.layers_.front();
        double expected;
        if (packed) {
          const double eb = 0.4, dp = 2.0e-3;
          // The cube root as the engine has always evaluated it for Sh
          // (film_correlation.h), so the packed Nu and Sh are one function.
          const double nu = 2.0 + 1.1 * std::pow(cp * mu / kg, 0.3333333) * std::pow(re, 0.6);
          expected = nu * kg / dp * 6.0 * (1.0 - eb) / dp;
        } else {
          const double open = d.gap, period = d.gap + d.sheet;
          const double dh = 2.0 * open, av = 2.0 / period;  // plates
          expected = 7.541 * kg / dh * av;
          EXPECT_NEAR(layer.GetHydraulicDiameter(), dh, 1e-18);
        }
        std::cout << "[film " << (packed ? "packed" : "laminate") << "] h a_v " << hp << " W/(m3 K), expected "
                  << expected << "\n";
        EXPECT_NEAR(hp, expected, 1e-10 * expected);
      });
    }
  }

  // 5b. Packed-bed film correlations (Wakao-Funazkri, Gunn, Ranz-Marshall;
  //     film_correlation.h): the layer's one choice sets both the particle Nu
  //     of the film heat model and the Sh of the Haynes coefficient, each
  //     written out here from its source and the logged gas properties, with
  //     Re on the superficial velocity and d_p.
  TEST(FilmCorrelation, HeatAndMassFollowTheLayerChoice) {
    const double eb = 0.4, dp = 2.0e-3;  // ChannelDeck's packed layer
    struct Case { FilmCorrelation film; const char* name; };
    const Case cases[] = {{FilmCorrelation::WakaoFunazkri, "Wakao-Funazkri"},
                          {FilmCorrelation::Gunn, "Gunn"},
                          {FilmCorrelation::RanzMarshall, "Ranz-Marshall"}};
    auto number = [&](FilmCorrelation c, double re, double x) {
      switch (c) {
      case FilmCorrelation::Gunn:
        return (7.0 - 10.0 * eb + 5.0 * eb * eb) * (1.0 + 0.7 * std::pow(re, 0.2) * std::cbrt(x)) +
               (1.33 - 2.4 * eb + 1.2 * eb * eb) * std::pow(re, 0.7) * std::cbrt(x);
      case FilmCorrelation::RanzMarshall: return 2.0 + 0.6 * std::sqrt(re) * std::cbrt(x);
      default: return 2.0 + 1.1 * std::pow(x, 0.3333333) * std::pow(re, 0.6);
      }
    };
    double hp_by[3] = {}, k_by[3] = {};
    for (const auto& c : cases) {
      ChannelDeck d;
      d.packed = true;
      d.ltne = true;
      d.haynes = true;
      d.film = c.film;
      d.molar_flow = 1.0e-2;
      d.run_time = 2.0;
      d.cells = 10;
      RunChannelDeck(d, [&](const UserData&, const CoupledSolver& solver) {
        const auto& log = solver.GetLog();
        auto mid = [&](VariableType v) { return Mid(log.GetVariable(0, 0, v).values); };
        const double hp = mid(VariableType::PARTICLE_HEAT_TRANSFER_COEFFICIENT);
        const double kg = mid(VariableType::THERMAL_CONDUCTIVITY), cp = mid(VariableType::HEAT_CAPACITY);
        const double mu = mid(VariableType::VISCOSITY), re = mid(VariableType::REYNOLDS_NUMBER);
        const double sc = mid(VariableType::SCHMIDT_NUMBER), sh = mid(VariableType::SHERWOOD_NUMBER);
        const double dm = Mid(log.GetComponentVariable(0, 0, ComponentVariableType::MOLECULAR_DIFFUSION, kTracer).values);
        const double k = Mid(log.GetComponentVariable(0, 0, ComponentVariableType::MASS_TRANSFER_COEFFICIENT, kTracer).values);
        ASSERT_GT(re, 1.0) << "the case should have a non-trivial Reynolds number";

        const double nu_expected = number(c.film, re, cp * mu / kg);
        const double hp_expected = nu_expected * kg / dp * 6.0 * (1.0 - eb) / dp;
        const double sh_expected = number(c.film, re, sc);
        // Haynes: 1/k = r_p / (3 k_f) + tau r_p^2 / (15 eps_p D_m), k_f = Sh D_m / d_p.
        const double rp = 0.5 * dp, kf = sh_expected * dm / dp;
        const double k_expected = 1.0 / (rp / (3.0 * kf) + d.tortuosity * rp * rp / (15.0 * d.ep * dm));
        std::cout << "[film " << c.name << "] Re " << re << " Nu " << hp / (kg / dp * 6.0 * (1.0 - eb) / dp)
                  << " (expected " << nu_expected << "), Sh " << sh << " (expected " << sh_expected << ")\n";
        EXPECT_NEAR(hp, hp_expected, 1e-10 * hp_expected) << c.name;
        EXPECT_NEAR(sh, sh_expected, 1e-12 * sh_expected) << c.name;
        EXPECT_NEAR(k, k_expected, 1e-10 * k_expected) << c.name;
        if (c.film == FilmCorrelation::WakaoFunazkri)  // the engine's fixed Sh before the choice existed
          EXPECT_EQ(sh, 2 + 1.1 * std::pow(sc, 0.3333333) * std::pow(re, 0.6));
        hp_by[static_cast<int>(c.film)] = hp;
        k_by[static_cast<int>(c.film)] = k;
      });
    }
    // The choice matters: the three correlations give different coefficients,
    // and the single-sphere Ranz-Marshall gives the smallest film coefficient.
    EXPECT_LT(hp_by[2], hp_by[0]);
    EXPECT_LT(hp_by[2], hp_by[1]);
    EXPECT_LT(k_by[2], k_by[0]);
    EXPECT_NE(hp_by[0], hp_by[1]);
  }

  // 6. A linear-isotherm breakthrough through a parallel-passage contactor.
  //
  // For a dilute tracer with holdup e_t, interstitial velocity v = u / e_t,
  // dispersion D_e = e_b D_ax / e_t (the dispersive flux is e_b D_ax dc/dz),
  // a linear isotherm with the dimensionless capacity
  // delta = rho_b K R T / e_t (q = K p) and the LDF dq/dt = k (q* - q), the
  // Laplace transform of the column's response is exp(-(L/v) G(s) +
  // (D_e L / v^3) G(s)^2 + ...) with G(s) = s (1 + delta k / (s + k)), so
  //   mean      t_1   = (L / v) (1 + delta),
  //   variance  sigma^2 = 2 delta L / (v k) + 2 D_e L (1 + delta)^2 / v^3,
  // to first order in 1 / Pe (the moment results of Kucera 1965 and Suzuki
  // and Smith 1971 for the LDF model). The step response's moments are
  // t_1 = int (1 - F) dt and E[t^2] = int 2 t (1 - F) dt with F = y / y_feed.
  //
  // k is Glueckauf's slab value 3 eps_p D_m / (tau l^2) and D_ax Taylor-Aris,
  // so the test exercises both with the open-area holdup and bulk density.
  // The run uses enough cells that the scheme's numerical dispersion is a
  // small part of sigma^2; it is checked to fall when the grid is refined.
  TEST(StructuredTransfer, ParallelPassageBreakthroughMomentsMatchTheLdfModel) {
    struct Moments { double t1, var; };
    auto measure = [](int cells) {
      ChannelDeck d;
      d.cells = cells;
      d.sheet = 1.0e-3;        // l = 0.5 mm: a slow LDF, k = 14 1/s
      d.gap = 3.0e-4;
      d.molar_flow = 4.0e-3;
      d.henry_k = 2.4e-6;      // delta ~ 10
      d.run_time = 80.0;
      d.max_step = 0.05;
      Moments m{0.0, 0.0};
      RunChannelDeck(d, [&](const UserData&, const CoupledSolver& solver) {
        const auto& log = solver.GetLog();
        auto [times, y_out] = TopFaceMoleFraction(log, 0, 0, kTracer);
        ASSERT_NEAR(y_out.back(), d.y_feed, 1e-4 * d.y_feed) << "the tracer has not broken through";
        std::vector<double> unbroken(times.size()), weighted(times.size());
        for (size_t n = 0; n < times.size(); ++n) {
          unbroken[n] = 1.0 - y_out[n] / d.y_feed;
          weighted[n] = 2.0 * times[n] * unbroken[n];
        }
        m.t1 = TimeIntegral(times, unbroken);
        m.var = TimeIntegral(times, weighted) - m.t1 * m.t1;
      });
      return std::make_pair(m, d);
    };

    const auto [coarse, d] = measure(120);
    const auto fine = measure(240).first;

    // The analytic moments from the geometry (re-derived here).
    const double period = d.gap + d.sheet;
    const double ofa = d.gap / period;                   // open frontal area
    const double et = ofa + (1.0 - ofa) * d.ep;          // gas holdup
    const double rho_b = (1.0 - ofa) * d.rho_s;          // bulk density
    const double area = 0.25 * M_PI * d.diameter * d.diameter;
    const double u = d.molar_flow * kR * kT / (kPout * area);  // superficial velocity
    const double v = u / et;
    const double delta = rho_b * d.henry_k * kR * kT / et;
    const double l = 0.5 * d.sheet;
    const double k = 3.0 * d.ep * d.dm_ldf / (d.tortuosity * l * l);  // Glueckauf (1955) slab
    const double U = u / ofa, W = d.gap;
    const double D_ax = d.dm_dispersion + W * W * U * U / (210.0 * d.dm_dispersion);  // Aris (1956)
    const double De = ofa * D_ax / et;
    const double t1 = d.length / v * (1.0 + delta);
    const double var_ldf = 2.0 * delta * d.length / (v * k);
    const double var_disp = 2.0 * De * d.length * (1.0 + delta) * (1.0 + delta) / (v * v * v);
    const double var = var_ldf + var_disp;

    std::cout << "[parallel passage] delta " << delta << ", k " << k << " 1/s, D_ax " << D_ax << " m2/s\n"
              << "  t1: analytic " << t1 << " s, 120 cells " << coarse.t1 << ", 240 cells " << fine.t1 << "\n"
              << "  sigma^2: analytic " << var << " s^2 (LDF " << var_ldf << ", dispersion " << var_disp
              << "), 120 cells " << coarse.var << ", 240 cells " << fine.var << "\n";
    // Guards: the LDF and dispersion terms both matter; the front is retained.
    ASSERT_GT(delta, 5.0);
    ASSERT_GT(var_ldf, 0.5 * var);
    ASSERT_GT(var_disp, 0.02 * var);
    // The mean is the holdup and capacity: exact up to the quadrature.
    EXPECT_NEAR(fine.t1, t1, 2e-3 * t1);
    // The variance converges to the analytic value as the grid is refined
    // (numerical dispersion falls), and the finer grid is within 5 %.
    EXPECT_LT(std::abs(fine.var - var), std::abs(coarse.var - var));
    EXPECT_NEAR(fine.var, var, 0.05 * var);
  }

  // 7. The operating limits warn when a structured layer's channel Reynolds
  //    number passes 2000, and not at laminar flow.
  TEST(StructuredTransfer, ChannelReynoldsWarning) {
    for (double flow : {2.0e-3, 2.0}) {
      ChannelDeck d;
      d.kind = LayerKind::Monolith;
      d.molar_flow = flow;
      d.run_time = 1.0;
      d.cells = 10;
      d.max_step = 0.05;
      RunChannelDeck(d, [&](const UserData& data, const CoupledSolver& solver) {
        const auto report = CheckOperatingLimits(data, solver.GetLog());
        size_t warnings = 0;
        for (const auto& w : report.warnings) warnings += w.limit == OperatingLimit::ChannelReynolds;
        ASSERT_FALSE(report.steps.empty());
        const double re = report.steps.front().max_channel_reynolds;
        std::cout << "[channel Re] flow " << flow << " mol/s: max Re " << re << ", warnings " << warnings << "\n";
        if (re > kLaminarChannelReynolds) EXPECT_EQ(warnings, 1u);
        else EXPECT_EQ(warnings, 0u);
        if (flow > 1.0) EXPECT_GT(re, kLaminarChannelReynolds) << "the high-flow case must pass the limit";
        else EXPECT_LT(re, 100.0);
      });
    }
  }

} // namespace adsim
