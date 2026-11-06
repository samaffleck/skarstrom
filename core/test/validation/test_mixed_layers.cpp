/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Solved tests for mixed adsorbent layers: one layer packed with several
// solids at set mass fractions w_k.
//
// The model (porous_media.cpp, domain.cpp, adsorbent_equations.cpp): within a
// cell all solids share the interparticle gas and one solid temperature; each
// solid k carries its own loading q_k at its own bulk density
//
//     rho_b,k = (1 - eb) rho_mix w_k,   rho_mix = 1 / sum_k (w_k / rho_k)
//
// (rho_mix is the particle density of the blend: mass fractions turn into
// solid-volume fractions phi_k = (w_k / rho_k) rho_mix). The particle voidage
// and the conductivity mix on phi_k, the heat capacity on w_k. The logged
// loading is the mass-weighted sum sum_k w_k q_k [mol per kg of solid].
//
// The existing coverage stops at the property mixing (MixedAdsorbentProperties
// unit tests) and one integration test comparing an outlet curve. Here the
// solved model is checked against references written from the physics:
//
//  1. A layer of A and an identical copy A' (50/50 and 25/75) gives the same
//     solution as a pure-A layer, in every field and cell, in LTE and LTNE.
//  2. A mixed with an inert solid of another density, voidage, heat capacity
//     and conductivity gives the same solution as pure A whose density,
//     voidage, heat capacity and conductivity are set to the blend's by hand
//     (the bulk density of A is unchanged, and the inert's heat capacity is
//     added to A's per kg of A), in LTE and LTNE.
//  3. The linear-isotherm stoichiometric breakthrough time of a two-adsorbent
//     layer is the column's mass-weighted capacity over the feed rate.
//  4. Mole balance (open, non-isothermal breakthrough) and energy balance
//     (closed bed heated to desorb from two different adsorbents) close.
//
// None of this code shares anything with core/src: the mixing rules are
// re-derived here from the definitions above.

#include "validation_helpers.h"
#include "adsim/equations/heat_transfer_model.h"

#include <algorithm>
#include <array>
#include <functional>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace adsim {

  using namespace adsim::validation;

  namespace {

    const std::string kO2 = "Oxygen";
    const std::string kN2 = "Nitrogen";

    // One solid of the blend. O2 follows a single-site Langmuir isotherm
    // q = m b e^{Q/RT} p / (1 + b e^{Q/RT} p), N2 a Henry isotherm
    // q = k e^{Q/RT} p; both with constant-coefficient LDF kinetics.
    struct Solid {
      double w = 1.0;       // mass fraction in the layer [-]
      double rho = 1100.0;  // particle density [kg/m3]
      double cp = 900.0;    // [J/(kg.K)]
      double k = 0.5;       // [W/(m.K)]
      double ep = 0.35;     // particle voidage [-]
      bool adsorbs = true;  // false: an inert solid (no isotherm at all)
      double m_O2 = 1.0, b_O2 = 1.2e-8, Q_O2 = 15000.0, dH_O2 = 15000.0, kldf_O2 = 0.3;
      double k_N2 = 1.0e-9, Q_N2 = 10000.0, dH_N2 = 10000.0, kldf_N2 = 1.0;
    };

    // A non-isothermal O2/N2 breakthrough on one (possibly mixed) layer:
    // nonlinear isotherm, heats of adsorption, wall heat transfer, axial
    // dispersion and a solid-conductivity-dependent heat dispersion
    // (Yagi-Kunii, which reads the layer's mixed ks).
    struct MixedDeck {
      std::vector<Solid> solids{Solid{}};
      bool ltne = false;
      int n_cells = 30;
      double length = 0.5;        // [m]
      double diameter = 0.05;     // [m]
      double eb = 0.4;            // [-]
      double molar_flow = 0.005;  // [mol/s]
      double y_feed = 0.2;        // O2 [-]
      double run_time = 40.0;     // [s]: the O2 front is mid-column at the end
      double reltol = 1e-6;
      double max_step = 0.5;
    };

    constexpr double kT0 = 298.0;  // [K]
    constexpr double kP0 = 1.0e5;  // [Pa]

    void BuildMixedDeck(UserData& data, const MixedDeck& d) {
      data.bed.wall_.thickness_ = 0.005;
      data.bed.wall_.heat_capacity_ = 500.0;
      data.bed.wall_.density_ = 7800.0;
      data.bed.wall_.thermal_conductivity_ = 16.0;
      data.bed.ambient_temperature_ = kT0;
      std::get<CylinderBed>(data.bed.shape_).diameter = d.diameter;

      data.fluid.AddComponent(kN2);
      data.fluid.AddComponent(kO2);

      const auto lid = data.bed.layers_.Add("Blend");
      auto& layer = data.bed.layers_.Get(lid);
      layer.momentum_model_ = MomentumModel::KarmenCozeny;
      layer.domain_.length = d.length;
      layer.domain_.n_cells = d.n_cells;
      layer.eb_ = d.eb;
      layer.particle_diameter_ = 0.002;

      for (size_t n = 0; n < d.solids.size(); ++n) {
        const Solid& s = d.solids[n];
        const auto id = data.adsorbent_directory.Add("Solid " + std::to_string(n + 1));
        auto& ads = data.adsorbent_directory.Get(id);
        ads.solid_density_ = s.rho;
        ads.heat_capacity_ = s.cp;
        ads.thermal_conductivity_ = s.k;
        ads.ep_ = s.ep;
        if (s.adsorbs) {
          auto lang = std::make_unique<LangmuirIsotherm>();
          lang->params_[kO2].m1 = s.m_O2;
          lang->params_[kO2].b1 = s.b_O2;
          lang->params_[kO2].q1 = s.Q_O2;
          ads.adsorption_source_.AddIsotherm({kO2}, std::move(lang));
          auto henry = std::make_unique<HenryIsotherm>();
          henry->params_[kN2].k = s.k_N2;
          henry->params_[kN2].Q = s.Q_N2;
          ads.adsorption_source_.AddIsotherm({kN2}, std::move(henry));
          for (const auto& [c, kldf] : {std::pair{kO2, s.kldf_O2}, std::pair{kN2, s.kldf_N2}}) {
            auto ldf = std::make_unique<LinearDrivingForce>();
            auto mtc = std::make_unique<MassTransferCoefficientConstant>();
            mtc->constant_mtc_ = kldf;
            ldf->mtc_model_ = std::move(mtc);
            ads.adsorption_source_.mass_transfer_models_[c] = std::move(ldf);
          }
          ads.adsorption_source_.heat_of_adsorption_[kO2] = s.dH_O2;
          ads.adsorption_source_.heat_of_adsorption_[kN2] = s.dH_N2;
        }
        layer.composition_[id] = s.w;
      }

      auto disp = std::make_unique<DispersionModelConstant>();
      disp->dispersion_coefficient_ = 1e-4;
      layer.equations_.dispersion_model_ = std::move(disp);
      layer.equations_.heat_dispersion_model_ =
        CreateAxialHeatDispersionModel(AxialHeatDispersionType::YAGI_AND_KUNII);
      auto ihtc = std::make_unique<ConstantInternalHeatTransferModel>();
      ihtc->heat_transfer_coefficient_ = 20.0;
      layer.equations_.internal_heat_transfer_model_ = std::move(ihtc);
      auto ehtc = std::make_unique<ConstantExternalHeatTransferModel>();
      ehtc->heat_transfer_coefficient_ = 5.0;
      layer.equations_.external_heat_transfer_model_ = std::move(ehtc);

      if (d.ltne) {
        data.bed.equations_.thermal_equilibrium_mode_ = ThermalEquilibriumMode::NON_LOCAL;
        auto php = std::make_unique<ConstantParticleHeatTransferModel>();
        php->heat_transfer_coefficient_ = 2.0e3;  // [W/(m3.K)]: Ts lags T by ~1 K
        data.bed.equations_.particle_heat_transfer_model_ = std::move(php);
      }

      data.ic.pressure = kP0;
      data.ic.temperature = kT0;
      data.ic.wall_temperature = kT0;
      data.ic.velocity = 0.0;
      data.ic.mole_fraction[kN2] = 1.0;
      data.ic.mole_fraction[kO2] = 0.0;

      data.cycle.AddStep("Feed");
      Step& step = data.cycle.GetStep("Feed");
      step.step_time_ = d.run_time;

      auto inlet = std::make_unique<BoundaryStream>();
      inlet->name_ = "Feed inlet";
      auto flow = std::make_unique<MolarFlowBoundaryValue>();
      static_cast<ConstantBoundaryValue*>(flow->profile_.get())->constant_value_ = d.molar_flow;
      inlet->boundary_flow_ = std::move(flow);
      inlet->mole_fraction[kN2] = 1.0 - d.y_feed;
      inlet->mole_fraction[kO2] = d.y_feed;
      inlet->temperature_ = kT0;
      auto inlet_bc = std::make_unique<FlowBoundaryCondition>();
      inlet_bc->direction_ = BoundaryDirection::INLET;
      inlet_bc->location_ = BoundaryLocation::BOTTOM;
      inlet_bc->streams_.push_back(std::move(inlet));
      step.bottom_boundary_condition_ = std::move(inlet_bc);

      auto outlet_bc = std::make_unique<PressureBoundaryCondition>();
      outlet_bc->direction_ = BoundaryDirection::OUTLET;
      outlet_bc->location_ = BoundaryLocation::TOP;
      auto p_out = std::make_unique<ConstantBoundaryValue>();
      p_out->constant_value_ = kP0;
      outlet_bc->pressure_ = std::move(p_out);
      outlet_bc->mole_fraction_[kN2] = 1.0 - d.y_feed;
      outlet_bc->mole_fraction_[kO2] = d.y_feed;
      step.top_boundary_condition_ = std::move(outlet_bc);

      data.cycle.number_of_cycles_ = 1;
      data.settings.reltol = d.reltol;
      data.settings.max_step = d.max_step;
    }

    // Blend properties re-derived from the definitions in the file comment.
    struct Blend {
      double rho_mix = 0.0;  // particle density of the blend [kg/m3]
      std::vector<double> phi;  // solid-volume fractions [-]
      double ep = 0.0, cp = 0.0, k = 0.0;
    };

    Blend BlendOf(const std::vector<Solid>& solids) {
      Blend b;
      double v = 0.0;
      for (const auto& s : solids) v += s.w / s.rho;
      b.rho_mix = 1.0 / v;
      for (const auto& s : solids) {
        const double phi = (s.w / s.rho) * b.rho_mix;
        b.phi.push_back(phi);
        b.ep += phi * s.ep;
        b.k += phi * s.k;
        b.cp += s.w * s.cp;
      }
      return b;
    }

    // ------------------------------------------------------------------------
    // Field comparison (the method of test_metamorphic_relations.cpp)
    // ------------------------------------------------------------------------

    // End-of-step interior-cell snapshots of every field: "P", "T", "Tw",
    // "Ts" (LTNE), "y:<c>", "q:<c>". q is multiplied by `q_scale`.
    using Fields = std::map<std::string, std::vector<double>>;

    Fields CollectFields(const EquationLog& log, bool ltne, double q_scale) {
      Fields f;
      auto interior = [](const std::vector<double>& s) {
        return std::vector<double>(s.begin() + 1, s.end() - 1);
      };
      f["P"] = interior(log.GetVariable(0, 0, VariableType::PRESSURE).values.back());
      f["T"] = interior(log.GetVariable(0, 0, VariableType::TEMPERATURE).values.back());
      f["Tw"] = interior(log.GetVariable(0, 0, VariableType::WALL_TEMPERATURE).values.back());
      if (ltne) f["Ts"] = interior(log.GetVariable(0, 0, VariableType::SOLID_TEMPERATURE).values.back());
      for (const auto& c : {kN2, kO2}) {
        f["y:" + c] = interior(log.GetComponentVariable(0, 0, ComponentVariableType::MOLE_FRACTION, c).values.back());
        auto q = interior(log.GetComponentVariable(0, 0, ComponentVariableType::ADSORBED_CONCENTRATION, c).values.back());
        for (double& v : q) v *= q_scale;
        f["q:" + c] = q;
      }
      return f;
    }

    // Absolute tolerances of the unknowns (docs/numerics.md).
    double AbsTol(const std::string& key) {
      switch (key[0]) {
        case 'P': return 1e-6;  // [Pa]
        case 'T': return 1e-4;  // [K]
        default: return 1e-8;   // mole fraction [-], loading [mol/kg]
      }
    }

    struct FieldDiff {
      double d = 0.0;       // largest |a - b| in the field's scale
      double weight = 0.0;  // one error-test weight, reltol |x| + atol, same scale
    };

    // Largest |a - b| per field: pressure and loading relative to the field's
    // largest magnitude in `a`, temperatures in K, mole fractions as they are.
    std::map<std::string, FieldDiff> FieldDistance(const Fields& a, const Fields& b, double reltol) {
      std::map<std::string, FieldDiff> out;
      for (const auto& [key, va] : a) {
        const auto& vb = b.at(key);
        EXPECT_EQ(va.size(), vb.size()) << key;
        double magnitude = 0.0;
        for (double v : va) magnitude = std::max(magnitude, std::abs(v));
        const bool absolute = key[0] == 'T' || key[0] == 'y';
        const double scale = absolute ? 1.0 : magnitude;
        FieldDiff f;
        for (size_t i = 0; i < va.size() && i < vb.size(); ++i)
          f.d = std::max(f.d, std::abs(va[i] - vb[i]) / scale);
        f.weight = (reltol * magnitude + AbsTol(key)) / scale;
        out[key] = f;
      }
      return out;
    }

    Fields RunDeck(const MixedDeck& deck, double q_scale, bool& ok) {
      UserData data;
      BuildMixedDeck(data, deck);
      std::string error;
      ok = data.EvaluateAllExpressions(error) && data.PerformChecks(error);
      EXPECT_TRUE(ok) << error;
      if (!ok) return {};
      CoupledSolver solver{data};
      ok = solver.IntegrateCycle(data.cycle, false) >= 0;
      EXPECT_TRUE(ok) << "Solver did not complete";
      if (!ok) return {};
      return CollectFields(solver.GetLog(), deck.ltne, q_scale);
    }

    // TOLERANCE. Two decks that pose the same semi-discrete system are two
    // adaptive integrations of it; once a floating-point detail differs (a
    // density formed as (1 - eb) rho_mix w instead of (1 - eb) rho, a sum over
    // two solids instead of one) they take different steps and each sits
    // within its global error E of the exact semi-discrete solution. So they
    // may differ by 2E and no more. E is measured per field by re-running the
    // base deck at reltol / 100; one error-test weight reltol |x| + atol is
    // added for fields whose E estimate is at round-off. A wrong mixing rule
    // is a modelling difference: it does not shrink with the tolerance and
    // shows up orders of magnitude above this bound.
    void ExpectSameSolution(const char* relation, const MixedDeck& base, double base_q_scale,
                            const MixedDeck& transformed, double transformed_q_scale) {
      bool ok_a = false, ok_fine = false, ok_b = false;
      MixedDeck fine_deck = base;
      fine_deck.reltol = base.reltol * 1e-2;
      const Fields a = RunDeck(base, base_q_scale, ok_a);
      const Fields fine = RunDeck(fine_deck, base_q_scale, ok_fine);
      const Fields b = RunDeck(transformed, transformed_q_scale, ok_b);
      ASSERT_TRUE(ok_a && ok_fine && ok_b) << relation << ": a run failed";

      // Guards: the O2 front is inside the bed, there is a thermal front, and
      // both adsorbates are loaded, so every field the relation covers moves.
      const auto& y = a.at("y:" + kO2);
      ASSERT_GT(*std::max_element(y.begin(), y.end()), 0.9 * base.y_feed) << relation << ": front has not entered";
      ASSERT_LT(*std::min_element(y.begin(), y.end()), 0.1 * base.y_feed) << relation << ": front has left the bed";
      const auto& T = a.at("T");
      ASSERT_GT(*std::max_element(T.begin(), T.end()) - kT0, 1.0) << relation << ": no thermal front";
      for (const auto& c : {kN2, kO2}) {
        const auto& q = a.at("q:" + c);
        ASSERT_GT(*std::max_element(q.begin(), q.end()), 1e-3) << relation << ": " << c << " not loaded";
      }

      const auto error = FieldDistance(a, fine, base.reltol);
      const auto discrepancy = FieldDistance(a, b, base.reltol);
      for (const auto& [key, e] : error) {
        const double d = discrepancy.at(key).d;
        const double tolerance = 2.0 * e.d + e.weight;
        std::cout << "[" << relation << "] " << key << ": discrepancy " << d
                  << ", integration error " << e.d << ", tolerance " << tolerance << "\n";
        EXPECT_LE(d, tolerance) << relation << ": field " << key << " differs by " << d
                                << ", more than the integrator tolerance " << tolerance;
      }
    }

  } // namespace

  // ===========================================================================
  // 1. A blend of A and an identical copy A' = pure A
  // ===========================================================================
  //
  // With rho_A' = rho_A the blend's particle density, voidage, conductivity
  // and heat capacity are A's, and rho_b,A + rho_b,A' = rho_b,A(pure). Both
  // copies see the same gas and temperature and start from the same loading,
  // so q_A = q_A' everywhere and the source terms sum to the pure layer's:
  // the two decks are the same semi-discrete system. Checked at 50/50 and at
  // 25/75 (an unequal split catches a rule that averages instead of weighting)
  // in LTE and LTNE (the solid energy equation sums the copies' sources and
  // heat capacities separately in each mode).

  TEST(MixedLayers, IdenticalCopyBlendEqualsPureLayer) {
    for (bool ltne : {false, true}) {
      MixedDeck pure;
      pure.ltne = ltne;
      for (double w : {0.5, 0.25}) {
        MixedDeck blend = pure;
        blend.solids = {Solid{}, Solid{}};
        blend.solids[0].w = w;
        blend.solids[1].w = 1.0 - w;
        const std::string name = std::string(ltne ? "LTNE" : "LTE") + " copy " +
                                 std::to_string(static_cast<int>(100 * w)) + "/" +
                                 std::to_string(static_cast<int>(100 * (1 - w)));
        ExpectSameSolution(name.c_str(), pure, 1.0, blend, 1.0);
      }
    }
  }

  // ===========================================================================
  // 2. A blended with an inert = pure A at the blend's properties
  // ===========================================================================
  //
  // Blend A (w_A = 0.6) with an inert solid I of another density, particle
  // voidage, heat capacity and conductivity. The equivalent pure-A layer has
  //
  //   rho'  = w_A rho_mix        so (1 - eb) rho' = rho_b,A of the blend
  //   ep'   = phi_A ep_A + phi_I ep_I       (the blend's total voidage)
  //   k'    = phi_A k_A  + phi_I k_I        (the blend's solid conductivity)
  //   cp'   = (w_A cp_A + w_I cp_I) / w_A   (the inert's heat capacity added,
  //                                          per kg of A)
  //
  // so rho_b' cp' = rho_b (w_A cp_A + w_I cp_I), the blend's solid heat
  // capacity per bed volume. The inert holds no loading and releases no heat,
  // and A's isotherm and kinetics are per kg of A, so q_A is the same in both
  // decks. The logged loading is per kg of all solid, w_A q_A in the blend, so
  // the blend's is divided by w_A before the comparison. A wrong mixing rule
  // (arithmetic density, mass-weighted voidage or conductivity) changes the
  // bed capacity, voidage or heat dispersion: with these properties each is a
  // change of several per cent.

  TEST(MixedLayers, InertDiluentEqualsPureLayerAtBlendProperties) {
    const double w_A = 0.6;
    Solid A;
    A.w = w_A;
    A.dH_O2 = 25000.0;  // a clear thermal front with less adsorbent per bed volume
    Solid I;
    I.w = 1.0 - w_A;
    I.adsorbs = false;
    I.rho = 2600.0;
    I.cp = 750.0;
    I.k = 2.0;
    I.ep = 0.05;
    const Blend blend = BlendOf({A, I});

    // Guard: the blend is far from both a mass-weighted and an
    // arithmetic-density blend, so the relation discriminates the rules.
    const double rho_arith = w_A * A.rho + (1 - w_A) * I.rho;
    ASSERT_GT(std::abs(rho_arith - blend.rho_mix) / blend.rho_mix, 0.15);
    ASSERT_GT(std::abs((w_A * A.ep + (1 - w_A) * I.ep) - blend.ep), 0.03);

    for (bool ltne : {false, true}) {
      MixedDeck mixed;
      mixed.ltne = ltne;
      mixed.run_time = 25.0;  // less adsorbent per bed volume: the front is faster
      // At reltol 1e-6 the y_O2 discrepancy is 1.9e-7, 0.75 of its bound; it
      // falls to 4.6e-8 at 1e-7 and 2e-8 to 4e-8 at 1e-8 (the absolute-
      // tolerance floor), so it is the integrator's, and 1e-7 leaves every
      // field a factor of about 3 below its bound. A wrong rule gives 1e-3.
      mixed.reltol = 1e-7;
      mixed.solids = {A, I};

      MixedDeck pure = mixed;
      Solid Ap = A;
      Ap.w = 1.0;
      Ap.rho = w_A * blend.rho_mix;
      Ap.ep = blend.ep;
      Ap.k = blend.k;
      Ap.cp = blend.cp / w_A;
      pure.solids = {Ap};

      ExpectSameSolution(ltne ? "LTNE inert diluent" : "LTE inert diluent", pure, 1.0, mixed, 1.0 / w_A);
    }
  }

  // ===========================================================================
  // 3. Stoichiometric breakthrough time of a two-adsorbent layer
  // ===========================================================================
  //
  // A dilute, isothermal O2 tracer on Henry isotherms q_k = K_k p (K_k =
  // k_k e^{Q_k/RT}) in a blend of two adsorbents. Integrating the tracer
  // balance over the breakthrough, the flux-weighted first moment of the
  // outlet curve is the column's tracer capacity over the feed rate:
  //
  //   t_st = V P [ et/(R T) + (1 - eb) rho_mix sum_k w_k K_k ] / n_dot,
  //   et = eb + (1 - eb) sum_k phi_k ep_k
  //
  // exactly, whatever the dispersion, the LDF rates or the grid. The two
  // solids differ in density (1000 vs 2500 kg/m3), voidage, capacity and
  // kinetics; the guards check that the answer differs from either pure
  // solid's and from an arithmetic-density blend's by far more than the
  // tolerance.
  //
  // TOLERANCE 0.2 %: the closed form takes the bed isobaric at the outlet
  // pressure, and the capacity is proportional to p; the pressure drop is
  // asserted below 0.05 %. The rest is the trapezoidal quadrature of the
  // outlet flux at max_step = 1 s and the truncated tail (y_out within 1e-3 of
  // the feed). Measured: +1.6e-4. The guards keep every wrong rule several
  // tolerances away: an arithmetic-density blend gives +14 %, a mass-weighted
  // particle voidage -1.9 %, either pure solid more than 20 %.

  TEST(MixedLayers, LinearIsothermBreakthroughTimeIsMassWeightedCapacity) {
    // Two Henry adsorbents for O2 only; N2 does not adsorb.
    struct Henry { double w, rho, ep, k, Q, kldf; };
    const std::array<Henry, 2> solids{{
      {0.3, 1000.0, 0.40, 4.0e-10, 20000.0, 0.5},
      {0.7, 2500.0, 0.10, 2.0e-10, 18000.0, 0.2},
    }};
    const double eb = 0.4, L = 1.0, D = 0.1, n_dot = 0.005, y_feed = 0.02, run_time = 700.0;

    UserData data;
    data.bed.wall_.thickness_ = 0.005;
    data.bed.ambient_temperature_ = kT0;
    std::get<CylinderBed>(data.bed.shape_).diameter = D;
    data.fluid.AddComponent(kN2);
    data.fluid.AddComponent(kO2);
    const auto lid = data.bed.layers_.Add("Blend");
    auto& layer = data.bed.layers_.Get(lid);
    layer.momentum_model_ = MomentumModel::KarmenCozeny;
    layer.domain_.length = L;
    layer.domain_.n_cells = 40;
    layer.eb_ = eb;
    layer.particle_diameter_ = 0.003;
    for (size_t n = 0; n < solids.size(); ++n) {
      const auto& s = solids[n];
      const auto id = data.adsorbent_directory.Add("Henry " + std::to_string(n + 1));
      auto& ads = data.adsorbent_directory.Get(id);
      ads.solid_density_ = s.rho;
      ads.ep_ = s.ep;
      ads.heat_capacity_ = 900.0;
      ads.thermal_conductivity_ = 0.0;
      auto henry = std::make_unique<HenryIsotherm>();
      henry->params_[kO2].k = s.k;
      henry->params_[kO2].Q = s.Q;
      ads.adsorption_source_.AddIsotherm({kO2}, std::move(henry));
      auto ldf = std::make_unique<LinearDrivingForce>();
      auto mtc = std::make_unique<MassTransferCoefficientConstant>();
      mtc->constant_mtc_ = s.kldf;
      ldf->mtc_model_ = std::move(mtc);
      ads.adsorption_source_.mass_transfer_models_[kO2] = std::move(ldf);
      ads.adsorption_source_.heat_of_adsorption_[kO2] = 0.0;  // isothermal
      layer.composition_[id] = s.w;
    }
    // Isothermal and dispersion-free: the feed is at the bed temperature and
    // nothing releases heat.
    auto ihtc = std::make_unique<ConstantInternalHeatTransferModel>();
    ihtc->heat_transfer_coefficient_ = 0.0;
    layer.equations_.internal_heat_transfer_model_ = std::move(ihtc);
    auto ehtc = std::make_unique<ConstantExternalHeatTransferModel>();
    ehtc->heat_transfer_coefficient_ = 0.0;
    layer.equations_.external_heat_transfer_model_ = std::move(ehtc);
    auto disp = std::make_unique<DispersionModelConstant>();
    disp->dispersion_coefficient_ = 0.0;
    layer.equations_.dispersion_model_ = std::move(disp);

    data.ic.pressure = kP0;
    data.ic.temperature = kT0;
    data.ic.wall_temperature = kT0;
    data.ic.velocity = 0.0;
    data.ic.mole_fraction[kN2] = 1.0;
    data.ic.mole_fraction[kO2] = 0.0;
    data.cycle.AddStep("Feed");
    Step& step = data.cycle.GetStep("Feed");
    step.step_time_ = run_time;
    auto inlet = std::make_unique<BoundaryStream>();
    auto flow = std::make_unique<MolarFlowBoundaryValue>();
    static_cast<ConstantBoundaryValue*>(flow->profile_.get())->constant_value_ = n_dot;
    inlet->boundary_flow_ = std::move(flow);
    inlet->mole_fraction[kN2] = 1.0 - y_feed;
    inlet->mole_fraction[kO2] = y_feed;
    inlet->temperature_ = kT0;
    auto inlet_bc = std::make_unique<FlowBoundaryCondition>();
    inlet_bc->direction_ = BoundaryDirection::INLET;
    inlet_bc->location_ = BoundaryLocation::BOTTOM;
    inlet_bc->streams_.push_back(std::move(inlet));
    step.bottom_boundary_condition_ = std::move(inlet_bc);
    auto outlet_bc = std::make_unique<PressureBoundaryCondition>();
    outlet_bc->direction_ = BoundaryDirection::OUTLET;
    outlet_bc->location_ = BoundaryLocation::TOP;
    auto p_out = std::make_unique<ConstantBoundaryValue>();
    p_out->constant_value_ = kP0;
    outlet_bc->pressure_ = std::move(p_out);
    outlet_bc->mole_fraction_[kN2] = 1.0 - y_feed;
    outlet_bc->mole_fraction_[kO2] = y_feed;
    step.top_boundary_condition_ = std::move(outlet_bc);
    data.cycle.number_of_cycles_ = 1;
    data.settings.reltol = 1e-8;
    data.settings.max_step = 1.0;

    std::string error;
    ASSERT_TRUE(data.EvaluateAllExpressions(error) && data.PerformChecks(error)) << error;
    CoupledSolver solver{data};
    ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";

    // --- analytic stoichiometric time ----------------------------------------
    double v = 0.0;
    for (const auto& s : solids) v += s.w / s.rho;
    const double rho_mix = 1.0 / v;
    double ep = 0.0, wK = 0.0, wK_arith = 0.0, rho_arith = 0.0;
    std::array<double, 2> K{};
    for (size_t n = 0; n < solids.size(); ++n) {
      const auto& s = solids[n];
      K[n] = s.k * std::exp(s.Q / (kR * kT0));
      ep += (s.w / s.rho) * rho_mix * s.ep;
      wK += s.w * K[n];
      rho_arith += s.w * s.rho;
    }
    const double V = 0.25 * M_PI * D * D * L;
    const double et = eb + (1 - eb) * ep;
    auto t_st = [&](double et_, double solid_capacity) {
      return V * kP0 * (et_ / (kR * kT0) + solid_capacity) / n_dot;
    };
    const double solid_capacity = (1 - eb) * rho_mix * wK;  // [mol/(m3.Pa)]
    const double t_expected = t_st(et, solid_capacity);

    // --- guards ---------------------------------------------------------------
    const double tol = 0.002 * t_expected;
    ASSERT_GT(solid_capacity, 2.0 * et / (kR * kT0)) << "Adsorbed capacity too small to test the blend";
    ASSERT_LT(t_expected, 0.5 * run_time) << "Run too short for full breakthrough";
    const double t_arith = t_st(et, (1 - eb) * rho_arith * wK);
    EXPECT_GT(std::abs(t_arith - t_expected), 10 * tol) << "Arithmetic-density blend not separated";
    double ep_mass = 0.0;
    for (const auto& s : solids) ep_mass += s.w * s.ep;
    const double t_ep_mass = t_st(eb + (1 - eb) * ep_mass, solid_capacity);
    EXPECT_GT(std::abs(t_ep_mass - t_expected), 5 * tol) << "Mass-weighted voidage not separated";
    for (size_t n = 0; n < solids.size(); ++n) {
      const double t_pure = t_st(eb + (1 - eb) * solids[n].ep, (1 - eb) * solids[n].rho * K[n]);
      EXPECT_GT(std::abs(t_pure - t_expected), 10 * tol) << "Pure solid " << n + 1 << " not separated";
    }

    // --- measured flux-weighted first moment ---------------------------------
    const auto& log = solver.GetLog();
    const auto P = log.GetVariable(0, 0, VariableType::PRESSURE);
    ASSERT_LT((P.values.back()[1] - kP0) / kP0, 5e-4) << "Pressure drop too large for the isobaric closed form";
    auto [times, y_out] = TopFaceMoleFraction(log, 0, 0, kO2);
    const auto n_out = log.GetFaceVariable(0, 0, VariableType::MOLAR_FLOW_RATE);
    ASSERT_FALSE(times.empty());
    ASSERT_NEAR(y_out.back(), y_feed, 1e-3 * y_feed) << "Tracer has not broken through";
    std::vector<double> unbroken(times.size(), 0.0);
    for (size_t n = 0; n < times.size(); ++n)
      unbroken[n] = 1.0 - n_out.values[n].back() * y_out[n] / (n_dot * y_feed);
    const double t_measured = TimeIntegral(times, unbroken);

    std::cout << "[stoichiometric time] measured " << t_measured << " s, analytic " << t_expected
              << " s, relative error " << (t_measured - t_expected) / t_expected
              << " (arithmetic-density blend " << t_arith << " s)\n";
    EXPECT_NEAR(t_measured, t_expected, tol)
      << "Stoichiometric time of the blend disagrees with its mass-weighted capacity";
  }

  // ===========================================================================
  // 4a. Mole balance on a mixed bed
  // ===========================================================================
  //
  // The non-isothermal breakthrough of tests 1-2 on a blend of two different
  // adsorbents (different densities, voidages, capacities, kinetics and
  // heats) and an inert. Per component, the moles the boundary streams log
  // entering the bed must equal the change of the bed's inventory,
  //
  //   gas:      sum_i V_i et P_i y_i / (R T_i)
  //   adsorbed: sum_i V_i (1 - eb) rho_mix q_i,   q_i = sum_k w_k q_k,i (logged)
  //
  // with et and rho_mix re-derived here. The adsorbed inventory is the sum of
  // each solid's rho_b,k q_k because rho_b,k = (1 - eb) rho_mix w_k.
  //
  // TOLERANCE 1e-5 of the largest of the inflow, outflow and inventory
  // change: the stream flows are integrated in time with the trapezoidal rule
  // over the solver's steps (max_step 0.1 s, reltol 1e-8). Measured: 2.6e-7
  // (N2) and 6e-8 (O2), in LTE and LTNE. A wrong per-solid bulk density
  // shifts the adsorbed inventory by tens of per cent of the O2 retained.

  namespace {
    double StreamMolesIn(const BoundaryLog& blog, const std::string& name, const std::string& component,
                         double sign) {
      const auto d = blog.GetComponentVariableData(0, 0, name, BoundaryComponentVariableType::MOLAR_FLOW, component);
      return sign * IntegrateVariableInTime(d.variable_timeseries, d.time);
    }

    MixedDeck TwoAdsorbentsAndInert() {
      MixedDeck d;
      Solid A;
      A.w = 0.45;
      Solid B;
      B.w = 0.35;
      B.rho = 1800.0;
      B.cp = 1000.0;
      B.k = 0.2;
      B.ep = 0.2;
      B.m_O2 = 2.0;
      B.b_O2 = 4.0e-9;
      B.Q_O2 = 17000.0;
      B.dH_O2 = 22000.0;
      B.kldf_O2 = 0.1;
      B.k_N2 = 3.0e-10;
      B.dH_N2 = 14000.0;
      Solid I;
      I.w = 0.20;
      I.adsorbs = false;
      I.rho = 2600.0;
      I.cp = 750.0;
      I.k = 2.0;
      I.ep = 0.05;
      d.solids = {A, B, I};
      d.reltol = 1e-8;
      d.max_step = 0.1;
      return d;
    }
  } // namespace

  TEST(MixedLayers, MoleBalanceClosesOnMixedBed) {
    for (bool ltne : {false, true}) {
      MixedDeck deck = TwoAdsorbentsAndInert();
      deck.ltne = ltne;
      UserData data;
      BuildMixedDeck(data, deck);
      std::string error;
      ASSERT_TRUE(data.EvaluateAllExpressions(error) && data.PerformChecks(error)) << error;
      CoupledSolver solver{data};
      ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";

      const auto& log = solver.GetLog();
      const auto& blog = solver.GetBoundaryLog();
      const auto& dom = log.GetDomain();
      const Blend blend = BlendOf(deck.solids);
      const double et = deck.eb + (1 - deck.eb) * blend.ep;
      const double rho_b = (1 - deck.eb) * blend.rho_mix;

      const auto P = log.GetVariable(0, 0, VariableType::PRESSURE);
      const auto T = log.GetVariable(0, 0, VariableType::TEMPERATURE);
      ASSERT_GT(*std::max_element(T.values.back().begin(), T.values.back().end()) - kT0, 1.0)
        << "Not a non-isothermal run";

      for (const auto& c : {kN2, kO2}) {
        const auto y = log.GetComponentVariable(0, 0, ComponentVariableType::MOLE_FRACTION, c);
        const auto q = log.GetComponentVariable(0, 0, ComponentVariableType::ADSORBED_CONCENTRATION, c);
        auto inventory = [&](size_t n, double& gas, double& ads) {
          gas = ads = 0.0;
          for (size_t i = 1; i + 1 < P.values[n].size(); ++i) {
            gas += dom.volume[i] * et * P.values[n][i] * y.values[n][i] / (kR * T.values[n][i]);
            ads += dom.volume[i] * rho_b * q.values[n][i];
          }
        };
        double gas0, ads0, gas1, ads1;
        inventory(0, gas0, ads0);
        inventory(P.values.size() - 1, gas1, ads1);
        const double accumulated = (gas1 - gas0) + (ads1 - ads0);

        const double in = StreamMolesIn(blog, "Feed inlet", c, +1.0);
        const double out = -StreamMolesIn(blog, "Feed Pressure Boundary " +
                                          std::string(BoundaryLocationList[(int)BoundaryLocation::TOP]), c, -1.0);
        const double scale = std::max({std::abs(in), std::abs(out), std::abs(accumulated)});
        const double residual = (in - out - accumulated) / scale;
        std::cout << "[mole balance " << (ltne ? "LTNE" : "LTE") << "] " << c << ": in " << in << ", out " << out
                  << ", gas " << gas1 - gas0 << ", adsorbed " << ads1 - ads0 << ", relative residual "
                  << residual << "\n";
        EXPECT_LT(std::abs(residual), 1e-5) << c << ": streams and inventory disagree";
        if (c == kO2) {
          // Guard: the adsorbed phase holds most of the O2 retained.
          ASSERT_GT(ads1 - ads0, 2.0 * (gas1 - gas0)) << "Adsorbed O2 too small to test the per-solid densities";
        }
      }
    }
  }

  // ===========================================================================
  // 4b. Energy balance on a mixed bed: heat-driven desorption from two solids
  // ===========================================================================
  //
  // A closed, adiabatic bed of two Henry adsorbents for N2 (different
  // capacities, heats, LDF rates, densities, heat capacities and voidages),
  // heated electrically. Both ends are shut, so every cell obeys the same
  // ODEs (the derivation of ClosedBedEnergyTest.HeatDrivenDesorption, with a
  // loading per solid):
  //
  //   dq_k/dt = kldf_k ( K_k(T) P - q_k )
  //   et dC/dt = - sum_k rho_b,k dq_k/dt
  //   [et C Cp + sum_k rho_b,k cp_k - et R C] dT/dt
  //        = et R T dC/dt + sum_k dH_k rho_b,k dq_k/dt + Q / V,     P = C R T
  //
  // integrated independently with RK4 at dt = 1e-3 s. The comparison closes the
  // energy balance (heater input = sensible heat of gas and both solids plus
  // the latent heat of each solid) and the mole balance between the gas and
  // each solid. Guards: each solid's latent heat is at least 10 % of the
  // heater input, so dropping or mis-weighting either source shows.
  //
  // TOLERANCE 1e-6 of each quantity's change over the run: the residual is
  // the integrator's (reltol 1e-10) and the RK4 reference's (dt^4, negligible).
  // Measured: 2e-10 (T), 4e-10 (P), 6e-10 (q). Mis-weighting either solid's
  // heat capacity, density or latent heat moves T by per cent.

  TEST(MixedLayers, EnergyBalanceClosesForHeatDrivenDesorptionFromTwoSolids) {
    struct Henry { double w, rho, ep, cp, k, Q, kldf; };
    const std::array<Henry, 2> solids{{
      {0.5, 1000.0, 0.30, 30.0, 4.5e-11, 25000.0, 0.5},
      {0.5, 2000.0, 0.10, 60.0, 1.5e-9, 20000.0, 0.1},
    }};
    const double eb = 0.4, L = 1.0, D = 0.1, heater_power = 70.0, heating_time = 200.0;

    UserData data;
    data.bed.wall_.thickness_ = 0.005;
    std::get<CylinderBed>(data.bed.shape_).diameter = D;
    data.fluid.AddComponent(kN2);
    data.fluid.AddComponent(kO2);
    const auto lid = data.bed.layers_.Add("Blend");
    auto& layer = data.bed.layers_.Get(lid);
    layer.momentum_model_ = MomentumModel::KarmenCozeny;
    layer.domain_.length = L;
    layer.domain_.n_cells = 20;
    layer.eb_ = eb;
    layer.particle_diameter_ = 0.002;
    for (size_t n = 0; n < solids.size(); ++n) {
      const auto& s = solids[n];
      const auto id = data.adsorbent_directory.Add("Henry " + std::to_string(n + 1));
      auto& ads = data.adsorbent_directory.Get(id);
      ads.solid_density_ = s.rho;
      ads.ep_ = s.ep;
      ads.heat_capacity_ = s.cp;
      ads.thermal_conductivity_ = 0.0;
      auto henry = std::make_unique<HenryIsotherm>();
      henry->params_[kN2].k = s.k;
      henry->params_[kN2].Q = s.Q;
      ads.adsorption_source_.AddIsotherm({kN2}, std::move(henry));
      auto ldf = std::make_unique<LinearDrivingForce>();
      auto mtc = std::make_unique<MassTransferCoefficientConstant>();
      mtc->constant_mtc_ = s.kldf;
      ldf->mtc_model_ = std::move(mtc);
      ads.adsorption_source_.mass_transfer_models_[kN2] = std::move(ldf);
      ads.adsorption_source_.heat_of_adsorption_[kN2] = s.Q;  // van 't Hoff consistent
      layer.composition_[id] = s.w;
    }
    auto ihtc = std::make_unique<ConstantInternalHeatTransferModel>();
    ihtc->heat_transfer_coefficient_ = 0.0;
    layer.equations_.internal_heat_transfer_model_ = std::move(ihtc);
    auto ehtc = std::make_unique<ConstantExternalHeatTransferModel>();
    ehtc->heat_transfer_coefficient_ = 0.0;
    layer.equations_.external_heat_transfer_model_ = std::move(ehtc);

    const double T0 = 300.0;
    data.ic.pressure = kP0;
    data.ic.temperature = T0;
    data.ic.wall_temperature = T0;
    data.ic.velocity = 0.0;
    data.ic.mole_fraction[kN2] = 1.0;
    data.ic.mole_fraction[kO2] = 0.0;
    data.cycle.AddStep("Regenerate");
    Step& step = data.cycle.GetStep("Regenerate");
    step.step_time_ = heating_time;
    auto bottom = std::make_unique<ClosedBoundaryCondition>();
    bottom->location_ = BoundaryLocation::BOTTOM;
    step.bottom_boundary_condition_ = std::move(bottom);
    auto top = std::make_unique<ClosedBoundaryCondition>();
    top->location_ = BoundaryLocation::TOP;
    step.top_boundary_condition_ = std::move(top);
    EnergySource heater;
    heater.enabled_ = true;
    heater.power_ = heater_power;
    step.energy_sources_[lid] = heater;
    data.cycle.number_of_cycles_ = 1;
    data.settings.reltol = 1e-10;
    data.settings.max_step = 1.0;

    std::string error;
    ASSERT_TRUE(data.EvaluateAllExpressions(error) && data.PerformChecks(error)) << error;
    CoupledSolver solver{data};
    ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";

    // --- independent lumped reference ----------------------------------------
    double v = 0.0;
    for (const auto& s : solids) v += s.w / s.rho;
    const double rho_mix = 1.0 / v;
    double ep = 0.0;
    std::array<double, 2> rho_b{};
    for (size_t n = 0; n < 2; ++n) {
      ep += (solids[n].w / solids[n].rho) * rho_mix * solids[n].ep;
      rho_b[n] = (1 - eb) * rho_mix * solids[n].w;
    }
    const double et = eb + (1 - eb) * ep;
    const double V = 0.25 * M_PI * D * D * L;
    const double Cp = MixtureMolarCp(data.fluid, {{kN2, 1.0}, {kO2, 0.0}});
    const double q_vol = heater_power / V;
    auto q_eq = [&](size_t n, double P, double T) {
      return solids[n].k * P * std::exp(solids[n].Q / (kR * T));
    };

    using State = std::array<double, 4>;  // T, C, q1, q2
    auto derivatives = [&](const State& y) {
      const double T = y[0], C = y[1];
      const double P = C * kR * T;
      std::array<double, 2> dq{};
      double sink = 0.0, latent = 0.0, solid_cp = 0.0;
      for (size_t n = 0; n < 2; ++n) {
        dq[n] = solids[n].kldf * (q_eq(n, P, T) - y[2 + n]);
        sink += rho_b[n] * dq[n];
        latent += solids[n].Q * rho_b[n] * dq[n];
        solid_cp += rho_b[n] * solids[n].cp;
      }
      const double dC = -sink / et;
      const double lhs = et * C * Cp + solid_cp - et * kR * C;
      const double rhs = et * kR * T * dC + latent + q_vol;
      return State{rhs / lhs, dC, dq[0], dq[1]};
    };
    State y{T0, kP0 / (kR * T0), q_eq(0, kP0, T0), q_eq(1, kP0, T0)};
    const State y0 = y;
    const double dt = 1e-3;
    for (int n = 0; n < static_cast<int>(heating_time / dt); ++n) {
      const auto k1 = derivatives(y);
      State y2{}, y3{}, y4{};
      for (int j = 0; j < 4; ++j) y2[j] = y[j] + 0.5 * dt * k1[j];
      const auto k2 = derivatives(y2);
      for (int j = 0; j < 4; ++j) y3[j] = y[j] + 0.5 * dt * k2[j];
      const auto k3 = derivatives(y3);
      for (int j = 0; j < 4; ++j) y4[j] = y[j] + dt * k3[j];
      const auto k4 = derivatives(y4);
      for (int j = 0; j < 4; ++j) y[j] += dt / 6.0 * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]);
    }
    const double T_ref = y[0];
    const double P_ref = y[1] * kR * y[0];
    // The logged loading is per kg of all solid: sum_k w_k q_k.
    const double qlog0 = solids[0].w * y0[2] + solids[1].w * y0[3];
    const double qlog_ref = solids[0].w * y[2] + solids[1].w * y[3];

    // --- guards ---------------------------------------------------------------
    const double E_in = heater_power * heating_time;
    for (size_t n = 0; n < 2; ++n) {
      const double latent = solids[n].Q * rho_b[n] * (y0[2 + n] - y[2 + n]) * V;
      std::cout << "[energy balance] solid " << n + 1 << " latent heat " << latent << " J of " << E_in << " J\n";
      ASSERT_GT(latent, 0.1 * E_in) << "Solid " << n + 1 << " latent heat too small to test its source";
    }

    // --- comparison -----------------------------------------------------------
    const auto& log = solver.GetLog();
    const auto& dom = log.GetDomain();
    const double T_sim = InteriorMean(log.GetVariable(0, 0, VariableType::TEMPERATURE).values.back(), dom.volume);
    const double P_sim = InteriorMean(log.GetVariable(0, 0, VariableType::PRESSURE).values.back(), dom.volume);
    const double q_sim = InteriorMean(
      log.GetComponentVariable(0, 0, ComponentVariableType::ADSORBED_CONCENTRATION, kN2).values.back(), dom.volume);
    std::cout << "[energy balance] T " << T_sim << " vs " << T_ref << " K (error / change "
              << (T_sim - T_ref) / (T_ref - T0) << "), P " << P_sim << " vs " << P_ref << " Pa ("
              << (P_sim - P_ref) / (P_ref - kP0) << "), q " << q_sim << " vs " << qlog_ref << " mol/kg ("
              << (q_sim - qlog_ref) / (qlog0 - qlog_ref) << ")\n";
    EXPECT_NEAR(T_sim, T_ref, 1e-6 * (T_ref - T0)) << "Temperature disagrees with the lumped energy balance";
    EXPECT_NEAR(P_sim, P_ref, 1e-6 * std::abs(P_ref - kP0)) << "Pressure disagrees with the lumped model";
    EXPECT_NEAR(q_sim, qlog_ref, 1e-6 * std::abs(qlog0 - qlog_ref)) << "Loading disagrees with the lumped model";
  }

} // namespace adsim
