/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Order-of-accuracy verification of the convective face-reconstruction schemes.
//
// test_numerical_verification measures the observed order of the UPWIND scheme
// on a DISCONTINUOUS exact solution (a travelling step), via the outlet-curve
// variance. This file complements that with a SMOOTH exact solution -- an inert
// tracer with physical dispersion switched on, whose outlet breakthrough curve
// is the analytic Ogata-Banks profile -- and measures order the way CFD codes
// verify a discretisation when no exact FINITE-domain solution exists:
//
//   * ACCURACY is checked against the analytic Ogata-Banks solution. A raw
//     pointwise comparison carries a grid-INDEPENDENT floor, because the finite
//     column with its outlet boundary is not exactly the semi-infinite
//     Ogata-Banks solution (an O(1/Pe) difference). That floor is a property of
//     the physics + boundary, identical for every scheme, so it CANCELS when the
//     schemes are compared to each other -- which is what the scheme-comparison
//     tests use it for.
//
//   * ORDER is measured by self-convergence (Richardson): for a scheme of order
//     p, the L2 difference between successive grids satisfies
//         || F_h - F_2h ||  ~  C (2^p - 1) h^p,
//     so log2( ||F_2h - F_4h|| / ||F_h - F_2h|| ) = p. This needs no exact
//     solution and is immune to the boundary floor above.
//
// Three checks:
//   1. Upwind converges at first order on the smooth solution.
//   2. WENO reduces numerical diffusion vs upwind (more accurate at every grid).
//   3. van Leer reduces numerical diffusion vs upwind (same).
// and a conservation check:
//   4. Every scheme conserves each species through the bed ends (the bed-end
//      faces must carry the boundary composition) on a coarse adsorbing
//      breakthrough fed from either end, and keeps the cell mole fractions
//      summing to one at a steep front.
//
// NOTE on observed order: the breakthrough is measured at the OUTLET FACE, and
// the face reconstruction at the domain boundary is first-order (the bed-end
// faces take the upwind value in every scheme: the ghost holds a face value,
// not a cell-centre value). So the OUTLET-measured self-convergence order is capped
// at ~1 for every scheme, including WENO/van Leer -- the interior high-order
// reconstruction cannot raise the order of a quantity read at a first-order
// boundary. What the higher-order schemes do deliver is a smaller error CONSTANT
// (less interior numerical diffusion), which is what tests 2 and 3 verify. Raising
// the outlet-measured order would require a high-order boundary reconstruction --
// a worthwhile future change, noted in the V&V dossier.
//
// The stoichiometric first moment (front position, fixed by capacity) must be
// preserved by every scheme -- a conservation cross-guard that a broken
// reconstruction would violate.

#include "breakthrough_case.h"

#include <limits>

namespace adsim {

  using namespace adsim::validation;

  namespace {

    // Fixed time grid the outlet curves are sampled onto so that runs at
    // different cell counts (and hence different internal time steps) can be
    // differenced pointwise.
    const std::vector<double>& CommonTimeGrid() {
      static const std::vector<double> grid = [] {
        std::vector<double> g;
        for (double t = 3.0; t <= 200.0; t += 0.5) g.push_back(t);
        return g;
      }();
      return grid;
    }

    struct SchemeRun {
      double u_int;               // interstitial velocity [m/s]
      double Pe;                  // u_int L / D_phys (inf when D_phys <= 0)
      double m1;                  // first moment (stoichiometric time) [s]
      double variance;            // outlet-curve variance [s^2] (numerical diffusion on a shock)
      double l2_vs_analytic;      // L2 of (F_sim - Ogata-Banks(D_phys)); 0 when D_phys <= 0
      std::vector<double> F;      // outlet C/C0 sampled on CommonTimeGrid()
    };

    // Run one inert-tracer breakthrough with the given face scheme. D_phys > 0
    // turns physical dispersion on (smooth Ogata-Banks solution, analytic L2
    // populated); D_phys <= 0 leaves dispersion off (the exact solution is a
    // travelling step, so the outlet variance is pure numerical diffusion).
    SchemeRun RunScheme(AdvectionScheme scheme, int n_cells, double D_phys) {
      BreakthroughCase c;
      c.n_cells = n_cells;
      c.adsorbing = false;
      c.length = 1.0;
      c.diameter = 0.1;
      c.molar_flow = 0.005;
      c.y_feed = 0.05;
      c.run_time = 300.0;
      c.max_step = 0.1;
      c.reltol = 1e-9;

      UserData data;
      Build(data, c);

      auto& layer = data.bed.layers_.layers_.front();
      auto disp = std::make_unique<DispersionModelConstant>();
      disp->dispersion_coefficient_ = D_phys > 0.0 ? D_phys : 0.0;
      layer.equations_.dispersion_model_ = std::move(disp);

      data.settings.advection_scheme_ = scheme;

      CoupledSolver solver{data};
      std::string error;
      EXPECT_TRUE(data.PerformChecks(error)) << error;
      EXPECT_GE(solver.IntegrateCycle(data.cycle, false), 0)
        << "Solver did not complete (scheme=" << static_cast<int>(scheme)
        << ", n_cells=" << n_cells << ")";

      const auto& log = solver.GetLog();
      const auto& dom = log.GetDomain();

      const auto P = log.GetVariable(0, 0, VariableType::PRESSURE);
      const double P_in = P.values.back()[1];
      const double P_out = P.values.back()[P.values.back().size() - 2];
      EXPECT_LT(std::abs(P_in - P_out) / P_out, 0.01)
        << "Pressure drop too large for a constant-velocity Ogata-Banks fit";

      const auto vel = log.GetVariable(0, 0, VariableType::VELOCITY);
      const auto& u_face = vel.values.back();
      double u_sum = 0.0; int u_cnt = 0;
      for (size_t i = 1; i + 1 < u_face.size(); ++i) { u_sum += u_face[i]; ++u_cnt; }
      const double et = dom.et[1];
      const double u_int = (u_sum / u_cnt) / et;

      auto [times, y_out] = OutletSeries(log, c.feed_from_top);
      EXPECT_NEAR(y_out.back(), c.y_feed, 1e-2 * c.y_feed)
        << "Breakthrough incomplete (scheme=" << static_cast<int>(scheme)
        << ", n_cells=" << n_cells << ")";

      // Sample F = C/C0 onto the common grid, and accumulate the L2 error vs the
      // analytic Ogata-Banks profile at the physical D.
      std::vector<double> Fvals(y_out.size());
      for (size_t i = 0; i < y_out.size(); ++i) Fvals[i] = y_out[i] / c.y_feed;

      const auto& tg = CommonTimeGrid();
      std::vector<double> F(tg.size());
      double sum_sq = 0.0;
      for (size_t i = 0; i < tg.size(); ++i) {
        F[i] = interp_linear(times, Fvals, tg[i]);
        if (D_phys > 0.0) {
          const double ref = OgataBanks(c.length, u_int, D_phys, tg[i]);
          sum_sq += (F[i] - ref) * (F[i] - ref);
        }
      }
      const double l2 = D_phys > 0.0 ? std::sqrt(sum_sq / tg.size()) : 0.0;

      // Moments of the outlet RTD from integrals of (1 - F): m1 = int (1-F) dt,
      // m2 = 2 int t (1-F) dt, variance = m2 - m1^2.
      std::vector<double> unbroken(times.size()), t_unbroken(times.size());
      for (size_t i = 0; i < times.size(); ++i) {
        unbroken[i] = 1.0 - Fvals[i];
        t_unbroken[i] = times[i] * unbroken[i];
      }
      const double m1 = TimeIntegral(times, unbroken);
      const double variance = 2.0 * TimeIntegral(times, t_unbroken) - m1 * m1;

      const double Pe = D_phys > 0.0 ? u_int * c.length / D_phys
                                     : std::numeric_limits<double>::infinity();
      return {u_int, Pe, m1, variance, l2, std::move(F)};
    }

    double L2Diff(const std::vector<double>& a, const std::vector<double>& b) {
      double s = 0.0;
      for (size_t i = 0; i < a.size(); ++i) s += (a[i] - b[i]) * (a[i] - b[i]);
      return std::sqrt(s / a.size());
    }

    void CheckSmoothFrontGuards(const SchemeRun& r) {
      ASSERT_GT(r.Pe, 5.0)   << "Peclet " << r.Pe << " too low: column fully mixed";
      ASSERT_LT(r.Pe, 200.0) << "Peclet " << r.Pe << " too high: front is a shock";
    }

    constexpr double kDphys = 2.0e-3;  // [m2/s]

  } // namespace

  // ===========================================================================
  // 1. Upwind: first order on the smooth solution (self-convergence)
  // ===========================================================================

  TEST(AdvectionSchemeOrder, UpwindIsFirstOrderOnSmoothSolution) {
    std::vector<SchemeRun> r;
    for (int n : {40, 80, 160, 320})
      r.push_back(RunScheme(AdvectionScheme::Upwind, n, kDphys));
    for (const auto& run : r) CheckSmoothFrontGuards(run);

    // First moment (stoichiometric time) is grid-invariant.
    for (size_t i = 1; i < r.size(); ++i)
      EXPECT_NEAR(r[i].m1, r[0].m1, 0.01 * r[0].m1)
        << "Stoichiometric time drifted with grid refinement (upwind)";

    // Correctness anchor: the finest solution matches the analytic Ogata-Banks
    // profile to within the finite-column O(1/Pe) floor.
    EXPECT_LT(r.back().l2_vs_analytic, 0.05)
      << "Finest upwind solution does not match analytic Ogata-Banks (L2 "
      << r.back().l2_vs_analytic << ")";

    // Observed order by self-convergence: log2 of the ratio of successive-grid
    // L2 differences. Two independent estimates from the four grids.
    const double d1 = L2Diff(r[0].F, r[1].F);  // 40 vs 80
    const double d2 = L2Diff(r[1].F, r[2].F);  // 80 vs 160
    const double d3 = L2Diff(r[2].F, r[3].F);  // 160 vs 320
    ASSERT_GT(d3, 0.0);
    const double order_a = std::log2(d1 / d2);
    const double order_b = std::log2(d2 / d3);
    EXPECT_NEAR(order_a, 1.0, 0.3)
      << "Upwind self-convergence order (40->80->160) is " << order_a;
    EXPECT_NEAR(order_b, 1.0, 0.3)
      << "Upwind self-convergence order (80->160->320) is " << order_b;
  }

  // On the SHOCK case (dispersion off) the exact solution is a travelling step
  // whose outlet variance is zero, so the entire simulated variance is numerical
  // diffusion. Comparing a scheme's outlet variance to upwind's is therefore a
  // direct, unambiguous measure of how much less numerical diffusion it carries
  // -- the decisive verification of a higher-order reconstruction. Both must also
  // preserve the stoichiometric first moment (front position) and converge toward
  // the zero-variance exact solution under refinement.
  void CheckLessDiffusiveThanUpwind(AdvectionScheme scheme, const char* name) {
    // Two grids: the variance ratio to upwind is far below the 0.75 bound on
    // every grid (measured, WENO / van Leer: 0.21 / 0.13 on 40 cells, 0.15 /
    // 0.09 on 80, 0.11 / 0.07 on 160), so a 160-cell level adds run time, not
    // strength.
    const std::vector<int> grids = {40, 80};
    std::vector<SchemeRun> up, sc;
    for (int n : grids) {
      up.push_back(RunScheme(AdvectionScheme::Upwind, n, /*D_phys=*/0.0));
      sc.push_back(RunScheme(scheme, n, /*D_phys=*/0.0));
    }

    for (size_t i = 0; i < grids.size(); ++i) {
      ASSERT_GT(sc[i].variance, 0.0);

      // Front arrives at its physical stoichiometric time L/u_int: the scheme
      // must sharpen the front without moving it. (A few percent of slack covers
      // the effect of the [0,1] face-value clamp on the sharp shock, which can
      // nudge the measured arrival time; total mass stays conserved regardless.)
      const double t_arrival = 1.0 / sc[i].u_int;  // L / u_int, L = 1 m
      EXPECT_NEAR(sc[i].m1, t_arrival, 0.03 * t_arrival)
        << name << " front arrival time " << sc[i].m1 << " s departs from L/u_int = "
        << t_arrival << " s at n_cells = " << grids[i];

      // Strictly less numerical diffusion than upwind at every grid.
      EXPECT_LT(sc[i].variance, up[i].variance)
        << name << " is not less diffusive than upwind at n_cells = " << grids[i]
        << " (" << name << " variance " << sc[i].variance << " vs upwind "
        << up[i].variance << ")";
    }

    // ... and materially so: at least a 25% reduction at the finest grid.
    EXPECT_LT(sc.back().variance, 0.75 * up.back().variance)
      << name << " numerical-diffusion reduction over upwind is marginal at the "
      << "finest grid (" << sc.back().variance << " vs upwind " << up.back().variance << ")";

    // Converges toward the zero-variance exact solution under refinement.
    EXPECT_LT(sc.back().variance, sc.front().variance)
      << name << " variance did not decrease under grid refinement ("
      << sc.front().variance << " -> " << sc.back().variance << ")";
  }

  // ===========================================================================
  // 2. WENO carries less numerical diffusion than upwind
  // ===========================================================================

  TEST(AdvectionSchemeOrder, WENOIsLessDiffusiveThanUpwind) {
    CheckLessDiffusiveThanUpwind(AdvectionScheme::WENO, "WENO");
  }

  // ===========================================================================
  // 3. van Leer carries less numerical diffusion than upwind
  // ===========================================================================

  TEST(AdvectionSchemeOrder, VanLeerIsLessDiffusiveThanUpwind) {
    CheckLessDiffusiveThanUpwind(AdvectionScheme::VanLeer, "van Leer");
  }

  // ===========================================================================
  // 4. Each species is conserved through the bed ends, and the cell mole
  //    fractions keep summing to one, whatever the scheme
  // ===========================================================================

  namespace {

    struct SchemeConsistency {
      double worst_balance_error = 0.0;  // max over components of |in - held| / fed
      std::string worst_component;
      double inlet_face_error = 0.0;     // |y_O2 at the inlet face - y_feed| (final state)
      double face_sum_error = 0.0;       // max over faces of |sum_i y_i,face - 1| (final state)
      double cell_sum_error = 0.0;       // max over snapshots and cells of |sum_i y_i - 1|
    };

    // A dispersion-free adsorbing breakthrough with the given scheme. The
    // balance of each component compares the moles that crossed the two bed
    // ends, as the engine's discrete boundary flux (GetTotalComponentMolesIn:
    // the upwind bed-end face value times the face flow, the flux the CSS
    // mass-balance check and the stream logs assume), with the change of the
    // gas and adsorbed inventory.
    // favourable: O2 on a strongly favourable Langmuir isotherm instead of the
    // case's linear Henry one, so the front self-sharpens into a shock.
    SchemeConsistency RunConsistency(AdvectionScheme scheme, bool feed_from_top, int n_cells,
                                     double y_feed, double molar_flow, double run_time,
                                     bool favourable = false, double reltol = 1e-8) {
      BreakthroughCase c;
      c.n_cells = n_cells;
      c.adsorbing = !favourable;
      c.feed_from_top = feed_from_top;
      c.y_feed = y_feed;
      c.molar_flow = molar_flow;
      c.run_time = run_time;
      c.max_step = 1.0;
      c.reltol = reltol;

      UserData data;
      Build(data, c);
      data.settings.advection_scheme_ = scheme;
      if (favourable) {
        auto& adsorbent = data.adsorbent_directory.Get(data.bed.layers_.layers_.front().composition_.begin()->first);
        auto isotherm = std::make_unique<LangmuirIsotherm>();
        isotherm->params_["Oxygen"].m1 = 3.0;      // [mol/kg]
        isotherm->params_["Oxygen"].b1 = 2.0e-7;   // [1/Pa]: b p ~ 10 at the feed partial pressure
        isotherm->params_["Oxygen"].q1 = 20000.0;  // [J/mol]
        adsorbent.adsorption_source_.AddIsotherm({"Oxygen"}, std::move(isotherm));
        auto ldf = std::make_unique<LinearDrivingForce>();
        auto mtc = std::make_unique<MassTransferCoefficientConstant>();
        mtc->constant_mtc_ = 0.5;
        ldf->mtc_model_ = ToResistances(std::move(mtc));  // generic model, as the example decks store it
        adsorbent.adsorption_source_.mass_transfer_models_["Oxygen"] = std::move(ldf);
        adsorbent.adsorption_source_.heat_of_adsorption_["Oxygen"] = 0.0;
      }

      CoupledSolver solver{data};
      std::string error;
      EXPECT_TRUE(data.PerformChecks(error)) << error;
      EXPECT_GE(solver.IntegrateCycle(data.cycle, false), 0)
        << "Solver did not complete (scheme=" << static_cast<int>(scheme) << ")";

      const auto& log = solver.GetLog();
      SchemeConsistency out;
      for (const std::string comp : {"Oxygen", "Nitrogen"}) {
        const double in = GetTotalComponentMolesIn(log, 0, 0, BoundaryLocation::BOTTOM, comp)
                        + GetTotalComponentMolesIn(log, 0, 0, BoundaryLocation::TOP, comp);
        const double held = GasMoles(log, 0, 0, comp, true) - GasMoles(log, 0, 0, comp, false)
                          + AdsorbedMoles(log, 0, 0, comp, true) - AdsorbedMoles(log, 0, 0, comp, false);
        const double fed = c.molar_flow * c.run_time * ((comp == "Oxygen") ? c.y_feed : 1.0 - c.y_feed);
        const double rel = std::abs(in - held) / fed;
        if (rel > out.worst_balance_error) {
          out.worst_balance_error = rel;
          out.worst_component = comp;
        }
      }

      const auto y_total = log.GetVariable(0, 0, VariableType::TOTAL_MOLE_FRACTION);
      for (const auto& snapshot : y_total.values)
        for (size_t i = 1; i + 1 < snapshot.size(); ++i)
          out.cell_sum_error = std::max(out.cell_sum_error, std::abs(snapshot[i] - 1.0));

      const auto& y_o2 = data.state.mole_fraction_face.at("Oxygen").var;
      const auto& y_n2 = data.state.mole_fraction_face.at("Nitrogen").var;
      for (size_t n = 0; n < y_o2.size(); ++n)
        out.face_sum_error = std::max(out.face_sum_error, std::abs(y_o2[n] + y_n2[n] - 1.0));
      const size_t inlet_face = feed_from_top ? y_o2.size() - 1 : 0;
      out.inlet_face_error = std::abs(y_o2[inlet_face] - c.y_feed);
      return out;
    }

    std::string Label(AdvectionScheme scheme, bool feed_from_top) {
      return "scheme " + std::to_string(static_cast<int>(scheme))
           + (feed_from_top ? ", fed from the top" : ", fed from the bottom");
    }

    constexpr AdvectionScheme kAllSchemes[] = {
      AdvectionScheme::Upwind, AdvectionScheme::WENO, AdvectionScheme::VanLeer };

  } // namespace

  // A coarse (10-cell) grid keeps the inlet cell far from the feed for much of
  // the run, so a bed-end face that does not carry the boundary composition
  // shows at the percent level: WENO and van Leer used to put the average of
  // the inlet ghost (the feed composition itself) and the first cell there.
  TEST(AdvectionSchemeConservation, EachSpeciesIsConservedThroughTheBedEnds) {
    for (auto scheme : kAllSchemes) {
      for (bool feed_from_top : { false, true }) {
        const auto r = RunConsistency(scheme, feed_from_top, 10, 0.2, 0.005, 300.0);
        const std::string label = Label(scheme, feed_from_top);
        EXPECT_LT(r.worst_balance_error, 1e-5)
          << label << ": " << r.worst_component << " moles in through the bed ends differ from the "
          << "inventory change by " << r.worst_balance_error << " of the moles fed";
        EXPECT_LT(r.inlet_face_error, 1e-6) << label << ": the inlet face does not carry the feed composition";
        EXPECT_LT(r.face_sum_error, 1e-6) << label << ": face mole fractions do not sum to one";
      }
    }
  }

  // Every species and the total have their own balance, so only consistent
  // fluxes keep the cell mole fractions summing to one. A limited
  // reconstruction linearised about a steep front weights the downwind cell;
  // built component by component it let a departure of the sum grow (van
  // Leer reached a sum of 90 in a CO2/N2 breakthrough on 13X; without the
  // face rescaling this case reaches 1.05 with WENO and van Leer). A
  // self-sharpening front: O2 on a strongly favourable Langmuir isotherm.
  // With consistent faces the sum stays at about 1e-12 for every scheme at
  // reltol 1e-8 (upwind 1.7e-12, WENO 5e-13, van Leer 1.2e-12). While the
  // high-order schemes clamped their bed-end faces to [0, 1] (WP20-1, see
  // EveryFaceSumsToItsUpwindCellSum below), the outlet cell's departure was
  // not carried out of the bed and reached 1.6e-6 (WENO) and 3.7e-7 (van
  // Leer) here; 1e-9 separates the two by more than two decades each way.
  TEST(AdvectionSchemeConservation, CellMoleFractionsKeepSummingToOneAtASteepFront) {
    for (auto scheme : kAllSchemes) {
      const auto r = RunConsistency(scheme, false, 50, 0.15, 0.05, 2500.0, /*favourable=*/true);
      const std::string label = Label(scheme, false);
      EXPECT_LT(r.cell_sum_error, 1e-9) << label << ": the cell mole fractions sum to 1 +- " << r.cell_sum_error;
      EXPECT_LT(r.face_sum_error, 1e-6) << label << ": face mole fractions do not sum to one";
      EXPECT_LT(r.worst_balance_error, 1e-5)
        << label << ": " << r.worst_component << " not conserved (" << r.worst_balance_error << " of the moles fed)";
    }
  }

  // The discrete identity behind the two tests above, checked face by face on
  // a frozen state. Write the species balance of cell n as
  //   eps_t dz (dy_i/dt + y_i (dP/dt / P - dT/dt / T)) = -(F_e y_i,e - F_w y_i,w) + D_i - S_i
  // and the total balance as
  //   eps_t dz (dP/dt / P - dT/dt / T)                 = -(F_e - F_w) - S,
  // with F the face molar flux over the cell concentration, y_i,e / y_i,w the
  // face mole fractions, D_i the dispersive term and S_i the source. Summing
  // the species balances and subtracting sum_n = sum_i y_i,n times the total:
  //   eps_t dz d(sum_n)/dt = -(sum_e - sum_n) F_e + (sum_w - sum_n) F_w
  //                          + sum_i D_i - (sum_i S_i - sum_n S),
  // where sum_e, sum_w are the sums of the face mole fractions. Nothing pulls
  // sum_n back to one; with upwind faces (sum_e = sum_n, sum_w = sum_n-1 for
  // u > 0) a departure is carried from cell to cell and out of the bed. So the
  // sum behaves as under upwinding only if EVERY face, the two bed-end faces
  // included, sums to its upwind cell's sum. A bed-end face whose components
  // are clamped to [0, 1] breaks this whenever the end cell holds a component
  // marginally outside the range (a pure-N2 outlet at 1 + 1e-9): the outlet
  // face then sums to one instead of sum_n, the -(sum_e - sum_n) F_e term
  // cancels the carry term, and the outlet cell's sum is no longer carried
  // out but integrates every solver-tolerance error for as long as the
  // outflow lasts (WP20-1: 2.7e-3 with WENO in a valve surge and blowdown).
  // The state below puts both end cells and ghosts outside [0, 1] and gives
  // every cell a sum other than one.
  TEST(AdvectionSchemeConservation, EveryFaceSumsToItsUpwindCellSum) {
    for (auto scheme : kAllSchemes) {
      for (double sign : { 1.0, -1.0 }) {
        BreakthroughCase c;
        c.n_cells = 12;
        UserData data;
        Build(data, c);
        data.settings.advection_scheme_ = scheme;
        data.Initialise();

        auto& state = data.state;
        auto& y_o2 = state.mole_fraction.at("Oxygen").var;
        auto& y_n2 = state.mole_fraction.at("Nitrogen").var;
        const int n_cells = static_cast<int>(y_o2.size());  // with the two ghosts
        ASSERT_EQ(n_cells, c.n_cells + 2);
        for (int k = 0; k < n_cells; ++k) {
          // A front in the middle of the bed, so the high-order faces differ
          // from upwind, and a sum off one by a few 1e-6 in every cell.
          const double x = (k - 0.5 * n_cells) / 2.0;
          y_o2[k] = 0.2 / (1.0 + std::exp(4.0 * x));
          y_n2[k] = 1.0 - y_o2[k] + 1e-6 * (1 + k % 3);
        }
        for (int k : { 0, 1, n_cells - 2, n_cells - 1 }) {
          y_o2[k] = -2e-7;     // marginally below 0
          y_n2[k] = 1 + 3e-6;  // marginally above 1
        }
        for (auto& v : state.velocity.var) v = sign * 0.1;  // |u| >> u_s: hard upwind

        for (const auto& comp : data.fluid.components_) UpdateMoleFractionFace(&data, comp);

        const auto& f_o2 = state.mole_fraction_face.at("Oxygen").var;
        const auto& f_n2 = state.mole_fraction_face.at("Nitrogen").var;
        const int n_faces = static_cast<int>(f_o2.size());
        ASSERT_EQ(n_faces, n_cells - 1);
        const std::string label = "scheme " + std::to_string(static_cast<int>(scheme))
                                + (sign > 0 ? ", flow up the bed" : ", flow down the bed");
        double worst = 0.0, worst_departure_from_upwind = 0.0;
        int worst_face = -1;
        for (int n = 0; n < n_faces; ++n) {
          const int up = sign > 0 ? n : n + 1;
          const double err = std::abs((f_o2[n] + f_n2[n]) - (y_o2[up] + y_n2[up]));
          if (err > worst) { worst = err; worst_face = n; }
          worst_departure_from_upwind = std::max(worst_departure_from_upwind, std::abs(f_o2[n] - y_o2[up]));
        }
        EXPECT_LT(worst, 1e-14) << label << ": face " << worst_face << " of " << n_faces
                                << " does not sum to its upwind cell's sum";
        // Guard: the high-order faces really are high order here.
        if (scheme != AdvectionScheme::Upwind)
          EXPECT_GT(worst_departure_from_upwind, 1e-3) << label << ": the reconstruction is not active";
      }
    }
  }

  // ===========================================================================
  // 5. A component absent from the feed and the bed costs no extra work
  // ===========================================================================
  //
  // O2 is absent from the feed and the initial bed (y = 0 exactly), and hot
  // N2 is fed into the bed, so the run has flow, a pressure start-up and a
  // thermal front, and the O2 cells hold nothing but round-off. Upwind carries
  // that round-off along and out of the bed. The high-order schemes used to
  // clamp every reconstructed face to [0, 1], which zeroed the face of a cell
  // marginally below zero: negative round-off was trapped (its face and the
  // face's derivative were zero) while positive round-off left, and the face
  // switching between clamped and free as Newton iterated stalled the
  // iteration. On 40 / 80 cells WENO took 737 / 4187 steps with 264 / 2026
  // convergence failures, and van Leer 400 / 503 steps with 93 / 128, against
  // upwind's 241 / 246 steps with none (WP20-4; 25 to 70 times the run time
  // in two tests with an absent component). With the face bounded to [0, 1]
  // widened to its two cells' values (BoundFace in adsorbent_equations.cpp)
  // both took exactly upwind's steps and no convergence failure. Since the
  // energy balance's face temperature follows the scheme too, the thermal
  // front is sharper under WENO and van Leer and takes a few more steps to
  // resolve (measured on 40 / 80 cells: WENO 245 / 286, van Leer 237 / 304,
  // upwind 241 / 246, still no convergence failure in any run). The bounds
  // below allow 30 % more steps than upwind (the old van Leer counts were
  // 1.3 and 1.6 times over it) and convergence failures on 2 % of its steps
  // (20 to 400 times below the old counts), which is the failure mode.
  namespace {
    struct AbsentComponentRun {
      IntegrationStats stats;
      double outlet_flow = 0.0;   // [mol/s] at the end of the run
      double max_abs_o2 = 0.0;    // max |y_O2| over the cells and the log
    };

    AbsentComponentRun RunAbsentComponent(AdvectionScheme scheme, int n_cells) {
      BreakthroughCase c;
      c.n_cells = n_cells;
      c.y_feed = 0.0;               // O2 absent from the feed (and the bed: Build)
      c.feed_temperature = 360.0;   // a thermal front through pure N2
      c.bed_temperature = 300.0;
      c.molar_flow = 0.01;
      c.run_time = 2000.0;
      c.max_step = 20.0;
      c.reltol = 1e-8;

      UserData data;
      Build(data, c);
      data.settings.advection_scheme_ = scheme;
      CoupledSolver solver{data};
      std::string error;
      EXPECT_TRUE(data.PerformChecks(error)) << error;
      EXPECT_GE(solver.IntegrateCycle(data.cycle, false), 0)
        << "Solver did not complete (scheme=" << static_cast<int>(scheme) << ")";

      AbsentComponentRun r;
      r.stats = solver.GetStatistics();
      const auto& log = solver.GetLog();
      const auto n_dot = log.GetFaceVariable(0, 0, VariableType::MOLAR_FLOW_RATE);
      r.outlet_flow = n_dot.values.back().back();
      const auto y = log.GetComponentVariable(0, 0, ComponentVariableType::MOLE_FRACTION, "Oxygen");
      for (const auto& snapshot : y.values)
        for (size_t i = 1; i + 1 < snapshot.size(); ++i)
          r.max_abs_o2 = std::max(r.max_abs_o2, std::abs(snapshot[i]));
      return r;
    }
  } // namespace

  TEST(AdvectionSchemeRobustness, AnAbsentComponentCostsNoExtraSteps) {
    for (int n_cells : {40, 80}) {
      const auto up = RunAbsentComponent(AdvectionScheme::Upwind, n_cells);
      // Guards: gas flows through the bed, and O2 really is absent. Absent O2
      // is an exact fixed point of its balance, but the nonlinear solve leaves
      // noise at the level of the mole fractions' absolute tolerance (1e-8):
      // measured max |y_O2| 1.9e-9 / 2.4e-9 (40 / 80 cells) on Linux and
      // 1.3e-8 on macOS, and 9e-23 / 2.1e-9 with a constant gas Cp. The guard
      // is ten times that tolerance, far below any real leak of O2.
      std::printf("[absent O2] %d cells: upwind max |y_O2| %.3e, %d steps\n", n_cells, up.max_abs_o2, (int)up.stats.num_steps);
      ASSERT_GT(up.outlet_flow, 0.5 * 0.01) << "no flow through the bed";
      ASSERT_LT(up.max_abs_o2, 1e-7) << "O2 is not absent from the bed";
      ASSERT_GT(up.stats.num_steps, 100);

      for (auto scheme : { AdvectionScheme::WENO, AdvectionScheme::VanLeer }) {
        const auto r = RunAbsentComponent(scheme, n_cells);
        const std::string label = "scheme " + std::to_string(static_cast<int>(scheme))
                                + ", " + std::to_string(n_cells) + " cells";
        EXPECT_LE(r.stats.num_steps, 1.3 * up.stats.num_steps)
          << label << ": " << r.stats.num_steps << " steps against upwind's " << up.stats.num_steps;
        EXPECT_LE(r.stats.num_conv_fails, up.stats.num_conv_fails + 0.02 * up.stats.num_steps)
          << label << ": " << r.stats.num_conv_fails << " Newton convergence failures against upwind's "
          << up.stats.num_conv_fails << " in " << up.stats.num_steps << " steps";
      }
    }
  }

} // namespace adsim
