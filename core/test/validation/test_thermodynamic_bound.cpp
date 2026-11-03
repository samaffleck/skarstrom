/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Thermodynamic bound: at cyclic steady state the energy a cycle spends on
// separation can be no less than the minimum (reversible) work of producing
// its product streams from its feed.
//
// THE BOUND. Take the plant (bed plus its pumps) at cyclic steady state,
// exchanging heat only with surroundings at T0 and handling ideal gases at
// T0. The second law for a steady-flow control volume gives
//
//     W  >=  sum_out n b  -  sum_in n b,     b = R T0 [ sum_i y_i ln y_i + ln(P / P0) ]
//
// (b is the molar flow exergy relative to the pure components at P0 and T0;
// the reference terms cancel because every component's moles balance at CSS).
// With all products and feeds at P0 this is the textbook minimum work of
// separation, R T0 [ sum_p n_p sum_i y_pi ln y_pi - n_f sum_i y_fi ln y_fi ].
// The tests use isothermal decks (EnergyBalanceMode::ISOTHERMAL at T0), so
// no stream carries thermal exergy and the bound is exactly the one above.
//
// The bound is checked as a chain of three inequalities, each a separate
// statement about the model:
//
//   W_electric  >=  W_rev  >=  dB_bed  >=  W_min.
//
//   * W_electric is the engine's reported energy: the time integral of every
//     endpoint's logged pump power (StreamData::electrical_power_), the
//     E_total of the energy KPI (docs/performance-metrics.md). The pump
//     model is polytropic, P = n/(n-1) P1 Q1 [(P2/P1)^((n-1)/n) - 1] / eta,
//     with eta = 0.7 and n = 1.4 (the defaults).
//   * W_rev is the reversible isothermal work of the same pumps on the same
//     logged flows, int n_dot R T0 |ln(P_bed / P_atm)| dt, computed here from
//     the logged pressure and flow. W_electric >= W_rev holds if the pump
//     model has eta <= 1 and a polytropic work no smaller than the isothermal
//     one (n >= 1).
//   * dB_bed is sum_out - sum_in of the exergy of the instantaneous streams,
//     each at the pressure it has at the plant boundary: the suction or
//     discharge pressure for a pumped endpoint, the bed-end pressure for an
//     unpumped one. W_rev >= dB_bed is the second law for the BED alone (the
//     bed does no work: dispersion, pressure drop, finite-rate adsorption
//     and mixing can only destroy exergy). A model whose transport or
//     adsorption created exergy would violate it.
//   * W_min is the minimum work for the cycle-averaged (mixed) products at
//     their delivery pressures: the light product at the feed step's outlet
//     pressure, everything else at P_atm (a vent or a pump discharge). It is
//     the "ideal-mixing" bound, and it is <= dB_bed because mixing a time-
//     varying stream and throttling it to P_atm only destroy exergy.
//
// Tolerance: the cycle is not exactly at CSS, so the component moles do not
// balance exactly. Each mole of imbalance moves the exergy sums by at most
// R T0 (max |ln y| + max |ln(P/P0)|); that, times the measured imbalance, is
// the slack allowed on the bed's second law. The guards: the imbalance is
// small, the cycle separates (W_min > 0 by a margin), and the pumps do work.

#include "psa_cycle_case.h"
#include "adsim/objects/compressor.h"

#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace adsim {

  using namespace adsim::validation;

  namespace {

    constexpr double kT0 = 298.0;      // [K] the isothermal temperature
    constexpr double kPatm = 1.0e5;    // [Pa] supply / vent / pump discharge pressure

    double MixingTerm(const std::map<std::string, double>& y) {
      double s = 0.0;
      for (const auto& [c, yi] : y) if (yi > 0.0) s += yi * std::log(yi);
      return s;
    }

    // One endpoint's log, with what the test needs to know about it.
    struct Endpoint {
      std::string name;
      const StreamData* log = nullptr;
      BoundaryLocation location = BoundaryLocation::TOP;
      bool pumped = false;
      double delivery_pressure = kPatm;  // where its mixed product is delivered (W_min)
    };

    struct Ledger {
      double W_electric = 0.0;
      double W_rev = 0.0;
      double dB_bed = 0.0;
      double W_min = 0.0;
      double W_min_separation = 0.0;   // the composition part of W_min only
      double imbalance = 0.0;          // largest |in - out| over the components [mol]
      double max_log_term = 0.0;       // max |ln y| + max |ln(P/P0)| over the streams
    };

    // A pressure-swing cycle on the shared PSA case, isothermal, with pumps:
    //   Pressurise: pressure inlet at the feed end ramping to P_high, from the
    //               atmosphere (a compressor when P_high > P_atm), top closed;
    //   Feed:       0.02 mol/s through a compressor from P_atm, light product
    //               at P_high from the top;
    //   Blowdown:   feed end to P_low (a vacuum pump to P_atm when
    //               P_low < P_atm), top closed.
    void BuildPumpedCycle(UserData& data, double P_high, double P_low) {
      BuildPsaCycle(data, 30);
      data.bed.equations_.energy_balance_mode_ = EnergyBalanceMode::ISOTHERMAL;
      data.ic.pressure = P_low;
      {
        auto bc = std::make_unique<PressureBoundaryCondition>();
        bc->direction_ = BoundaryDirection::INLET;
        bc->location_ = BoundaryLocation::BOTTOM;
        auto ramp = std::make_unique<LinearBoundaryValue>();  // from the bed pressure, over the step
        ramp->end_value_ = P_high;
        bc->pressure_ = std::move(ramp);
        bc->temperature_ = kT0;
        bc->mole_fraction_["Nitrogen"] = 0.79;
        bc->mole_fraction_["Oxygen"] = 0.21;
        auto compressor = std::make_unique<Compressor>();
        compressor->use_pump_ = P_high > kPatm;
        compressor->suction_pressure_ = kPatm;
        bc->pump_ = std::move(compressor);
        data.cycle.GetStep("Pressurise").bottom_boundary_condition_ = std::move(bc);
      }
      {
        Step& feed = data.cycle.GetStep("Feed");
        auto* inlet = dynamic_cast<FlowBoundaryCondition*>(feed.bottom_boundary_condition_.get());
        auto* stream = dynamic_cast<BoundaryStream*>(inlet->streams_.front().get());
        auto compressor = std::make_unique<Compressor>();
        compressor->use_pump_ = true;
        compressor->suction_pressure_ = kPatm;
        stream->pump_ = std::move(compressor);
        auto* top = dynamic_cast<PressureBoundaryCondition*>(feed.top_boundary_condition_.get());
        static_cast<ConstantBoundaryValue*>(top->pressure_.get())->constant_value_ = P_high;
      }
      {
        auto* bc = dynamic_cast<PressureBoundaryCondition*>(
          data.cycle.GetStep("Blowdown").bottom_boundary_condition_.get());
        auto ramp = std::make_unique<ExponentialBoundaryValue>();
        ramp->end_value_ = P_low;
        ramp->step_time_ = 120.0;
        bc->pressure_ = std::move(ramp);
        auto vacuum = std::make_unique<VacuumPump>();
        vacuum->use_pump_ = P_low < kPatm;
        vacuum->discharge_pressure_ = kPatm;
        bc->pump_ = std::move(vacuum);
      }
    }

    Ledger RunAndAccount(UserData& data, double P_high, int cycles) {
      Ledger L;
      data.cycle.number_of_cycles_ = cycles;
      data.settings.reltol = 1e-8;
      data.settings.max_step = 2.0;
      std::string error;
      EXPECT_TRUE(data.PerformChecks(error)) << error;
      CoupledSolver solver{data};
      EXPECT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";

      // Every endpoint of the last cycle (each step's log is its last run).
      std::vector<Endpoint> endpoints;
      for (auto& step : data.cycle.step_list_) {
        for (auto* bc : {step.top_boundary_condition_.get(), step.bottom_boundary_condition_.get()}) {
          if (auto* f = dynamic_cast<FlowBoundaryCondition*>(bc)) {
            for (auto& s : f->streams_) {
              auto* b = dynamic_cast<BoundaryStream*>(s.get());
              endpoints.push_back({step.step_name_ + " " + s->name_, &s->stream_data_, f->location_,
                                   b && b->pump_ && b->pump_->use_pump_, kPatm});
            }
          }
          if (auto* p = dynamic_cast<PressureBoundaryCondition*>(bc)) {
            const bool light_product = step.step_name_ == "Feed" && p->location_ == BoundaryLocation::TOP;
            endpoints.push_back({step.step_name_ + " pressure boundary", &p->stream_data_, p->location_,
                                 p->pump_ && p->pump_->use_pump_, light_product ? P_high : kPatm});
          }
        }
      }

      const double RT = kR * kT0;
      std::map<std::string, double> component_out;  // over the plant, out minus in [mol]
      double min_y = 1.0, max_lnP = 0.0;
      for (const auto& e : endpoints) {
        const StreamData& sd = *e.log;
        // The logged flow is axial (positive upwards); out of the bed is +axial
        // at the top and -axial at the bottom. (The location comes from the
        // boundary condition: a pressure boundary's log does not record it.)
        const double out_sign = e.location == BoundaryLocation::TOP ? 1.0 : -1.0;
        const auto& t = sd.time_series_;
        std::map<std::string, double> moles;
        double w_e = 0.0, w_rev = 0.0, b = 0.0;
        auto rates = [&](size_t k, double& w_rev_rate, double& b_rate) {
          const double n_out = out_sign * sd.molar_flow_[k];
          const double P_bed = sd.pressure_boundary_[k];
          std::map<std::string, double> y;
          for (const auto& [c, v] : sd.mole_fraction_boundary_) {
            y[c] = v[k];
            if (std::abs(n_out) > 1e-12) min_y = std::min(min_y, std::max(v[k], 1e-300));
          }
          max_lnP = std::max(max_lnP, std::abs(std::log(P_bed / kPatm)));
          w_rev_rate = e.pumped ? std::abs(n_out) * RT * std::abs(std::log(P_bed / kPatm)) : 0.0;
          const double P_plant = e.pumped ? kPatm : P_bed;
          b_rate = n_out * RT * (MixingTerm(y) + std::log(P_plant / kPatm));
        };
        for (size_t k = 1; k < t.size(); ++k) {
          const double dt = t[k] - t[k - 1];
          w_e += 0.5 * (sd.electrical_power_[k] + sd.electrical_power_[k - 1]) * dt;
          double r1, b1, r0, b0;
          rates(k, r1, b1);
          rates(k - 1, r0, b0);
          w_rev += 0.5 * (r1 + r0) * dt;
          b += 0.5 * (b1 + b0) * dt;
          for (const auto& [c, v] : sd.component_molar_flow_)
            moles[c] += out_sign * 0.5 * (v[k] + v[k - 1]) * dt;
        }
        L.W_electric += w_e;
        L.W_rev += w_rev;
        L.dB_bed += b;

        // Mixed product (or feed) at its delivery pressure.
        double n_total = 0.0;
        for (const auto& [c, n] : moles) { n_total += n; component_out[c] += n; }
        if (n_total != 0.0) {
          std::map<std::string, double> y_mixed;
          for (const auto& [c, n] : moles) y_mixed[c] = n / n_total;
          const double P_delivery = n_total > 0.0 ? e.delivery_pressure : kPatm;  // feeds come from P_atm
          L.W_min += n_total * RT * (MixingTerm(y_mixed) + std::log(P_delivery / kPatm));
          L.W_min_separation += n_total * RT * MixingTerm(y_mixed);
        }
        std::cout << "  " << e.name << (e.pumped ? " (pumped)" : "") << ": out " << n_total
                  << " mol, W_electric " << w_e << " J, W_rev " << w_rev << " J, B_out " << b << " J\n";
      }
      for (const auto& [c, n] : component_out) L.imbalance = std::max(L.imbalance, std::abs(n));
      L.max_log_term = std::abs(std::log(min_y)) + max_lnP;
      return L;
    }

    void CheckChain(const char* name, const Ledger& L) {
      const double slack = kR * kT0 * L.max_log_term * L.imbalance;
      std::cout << "[" << name << "] W_electric " << L.W_electric << " J >= W_rev " << L.W_rev
                << " J >= dB_bed " << L.dB_bed << " J >= W_min " << L.W_min << " J (separation part "
                << L.W_min_separation << " J); second-law efficiency W_min/W_electric "
                << L.W_min / L.W_electric << "; CSS imbalance " << L.imbalance << " mol, slack " << slack
                << " J\n";
      // Guards: at CSS, separating, and the pumps do work.
      ASSERT_LT(L.imbalance, 1e-3) << "The cycle is not at cyclic steady state";
      ASSERT_GT(L.W_rev, 0.0);
      ASSERT_GT(L.W_min, 10.0 * slack) << "The bound is not resolved by the CSS accuracy";
      EXPECT_GE(L.W_electric, L.W_rev) << "Pump power below the reversible isothermal work";
      EXPECT_GE(L.W_rev + slack, L.dB_bed) << "The bed creates exergy";
      EXPECT_GE(L.dB_bed + slack, L.W_min) << "Mixing / throttling the products created exergy";
      EXPECT_GE(L.W_electric, L.W_min) << "Separation energy below the minimum work";
    }

  } // namespace

  // Vacuum swing: feed at ~1 bar through a compressor (it only covers the
  // pressure drop), light product at 1 bar, blowdown to 0.2 bar through a
  // vacuum pump to 1 bar, repressurisation from the atmosphere. Every plant
  // stream is at 1 bar, so W_min is the pure separation work.
  TEST(ThermodynamicBound, VacuumSwingEnergyExceedsMinimumSeparationWork) {
    UserData data;
    BuildPumpedCycle(data, 1.0e5, 0.2e5);
    const Ledger L = RunAndAccount(data, 1.0e5, 20);
    CheckChain("VSA", L);
    EXPECT_NEAR(L.W_min, L.W_min_separation, 1e-12 * std::abs(L.W_min));  // no pressure term
    EXPECT_GT(L.W_min_separation, 0.0) << "The products are not purer than the feed";
  }

  // Pressure swing 1 -> 3 bar: the light product leaves at 3 bar, so W_min
  // includes its pressure exergy as well as the separation work.
  TEST(ThermodynamicBound, PressureSwingEnergyExceedsMinimumWork) {
    UserData data;
    BuildPumpedCycle(data, 3.0e5, 1.0e5);
    const Ledger L = RunAndAccount(data, 3.0e5, 20);
    CheckChain("PSA", L);
    EXPECT_GT(L.W_min_separation, 0.0) << "The products are not purer than the feed";
  }

} // namespace adsim
