/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// A connection replays what left the donor's bed into the receiver
// (StreamData::GetOutflow). For a donor that is a pressure boundary at the
// BOTTOM of the bed, the moles the receiver takes in must equal the moles the
// donor gave out, with the same sign.
//
// Derivation (conservation across a connection): the receiver must take in
// exactly what the donor's bed gave out. With the axial sign convention of
// StreamData::molar_flow_ (positive upwards), a bottom outlet that releases
// gas has axial flow < 0; its outflow is -axial > 0; a bottom inlet that
// receives it has axial flow +outflow > 0. So receiver_axial = -donor_axial,
// asserted to 0.1 % (the receiver replays the donor's log interpolated at its
// own time points, so the two trapezoid sums differ by the interpolation error
// only; measured 6e-6 relative for the smooth blowdown and 3e-5 for the step
// change after the fix).
//
// Deck: the shared PSA cycle (psa_cycle_case.h, isothermal) plus a Rinse step
// as long as the blowdown (120 s). The blowdown's bottom pressure boundary (an
// exponential ramp to 1 bar, so the replayed flow is smooth) is connected to a
// no-valve flow inlet at the bottom of the rinse; the rinse top is a 1 bar
// pressure outlet. Before the fix the pressure boundary's log kept the default
// TOP, the receiver drew -1.8516 mol out of its bed instead of taking in the
// +1.8517 mol the blowdown released, and with a step-change blowdown IDA
// failed at the start of the rinse.

#include "../validation/psa_cycle_case.h"
#include "adsim/boundary_conditions/connection.h"

#include <iomanip>

namespace adsim {

  using namespace adsim::validation;

  namespace {

  void BuildRinse(UserData& data, bool smooth_blowdown) {
    BuildPsaCycle(data, 20);
    data.bed.equations_.energy_balance_mode_ = EnergyBalanceMode::ISOTHERMAL;
    // Rinse: the blowdown's bottom outflow (donor, a pressure boundary at the
    // bottom) is fed back in at the bottom; the top is a 1 bar outlet.
    data.cycle.AddStep("Rinse");
    Step& rinse = data.cycle.GetStep("Rinse");
    rinse.step_time_ = 120.0;  // = blowdown
    auto stream = std::make_unique<BoundaryStream>();
    auto flow = std::make_unique<MolarFlowBoundaryValue>();
    static_cast<ConstantBoundaryValue*>(flow->profile_.get())->constant_value_ = 0.001;
    stream->boundary_flow_ = std::move(flow);
    stream->mole_fraction["Nitrogen"] = 0.79;
    stream->mole_fraction["Oxygen"] = 0.21;
    stream->temperature_ = 298.0;
    stream->name_ = "Rinse in";
    IBoundaryStream* rinse_in = stream.get();
    auto in = std::make_unique<FlowBoundaryCondition>();
    in->direction_ = BoundaryDirection::INLET;
    in->location_ = BoundaryLocation::BOTTOM;
    in->streams_.push_back(std::move(stream));
    rinse.bottom_boundary_condition_ = std::move(in);
    auto out = std::make_unique<PressureBoundaryCondition>();
    out->direction_ = BoundaryDirection::OUTLET;
    out->location_ = BoundaryLocation::TOP;
    out->pressure_ = std::make_unique<ConstantBoundaryValue>(1.0e5);
    out->mole_fraction_["Nitrogen"] = 0.79;
    out->mole_fraction_["Oxygen"] = 0.21;
    rinse.top_boundary_condition_ = std::move(out);
    data.cycle.EnsureIds();
    auto* donor = dynamic_cast<PressureBoundaryCondition*>(data.cycle.GetStep("Blowdown").bottom_boundary_condition_.get());
    ASSERT_NE(donor, nullptr);
    ASSERT_EQ(donor->location_, BoundaryLocation::BOTTOM);
    auto a = MakeEndpoint(data.cycle.GetStep("Blowdown"), *donor);
    auto b = MakeEndpoint(data.cycle.GetStep("Rinse"), *rinse_in);
    ASSERT_TRUE(MakeMutualConnection(a, b));
    if (smooth_blowdown) {
      auto e = std::make_unique<ExponentialBoundaryValue>();
      e->end_value_ = 1.0e5;
      e->step_time_ = 120.0;
      donor->pressure_ = std::move(e);
    }
  }

  } // namespace

  // One cycle: no replay has happened yet; the donor's log alone shows the
  // problem. Its outflow (what a connection replays) has the wrong sign.
  TEST(ConnectionSign, BottomPressureBoundaryLogRecordsItsBedEnd) {
    UserData data;
    BuildRinse(data, true);
    data.cycle.number_of_cycles_ = 1;
    data.settings.reltol = 1e-6;
    data.settings.max_step = 2.0;
    std::string err;
    ASSERT_TRUE(data.PerformChecks(err)) << err;
    CoupledSolver solver{data};
    ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0);
    auto* donor = dynamic_cast<PressureBoundaryCondition*>(data.cycle.GetStep("Blowdown").bottom_boundary_condition_.get());
    const auto& sd = donor->stream_data_;
    ASSERT_FALSE(sd.time_series_.empty());
    const double t_mid = 0.5 * (sd.time_series_.front() + sd.time_series_.back());
    std::cout << "axial flow " << sd.GetMolarFlow(t_mid) << " mol/s, GetOutflow " << sd.GetOutflow(t_mid) << " mol/s\n";
    ASSERT_LT(sd.GetMolarFlow(t_mid), 0.0) << "Blowdown gas leaves through the bottom (downwards)";
    EXPECT_EQ(sd.location_, BoundaryLocation::BOTTOM);
    EXPECT_GT(sd.GetOutflow(t_mid), 0.0) << "Gas leaving the bed must be a positive outflow";
  }

  // Runs two cycles (the second replays the first's blowdown into the rinse)
  // and checks the receiver took in what the donor released.
  void ExpectRinseTakesInTheBlowdown(bool smooth_blowdown) {
    UserData data;
    BuildRinse(data, smooth_blowdown);
    auto* donor = dynamic_cast<PressureBoundaryCondition*>(data.cycle.GetStep("Blowdown").bottom_boundary_condition_.get());
    IBoundaryStream* rinse_in = dynamic_cast<FlowBoundaryCondition*>(data.cycle.GetStep("Rinse").bottom_boundary_condition_.get())->streams_.front().get();
    data.cycle.number_of_cycles_ = 2;
    data.settings.reltol = 1e-6;
    data.settings.max_step = 2.0;
    std::string err;
    ASSERT_TRUE(data.PerformChecks(err)) << err;
    CoupledSolver solver{data};
    ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0);

    auto moles = [](const StreamData& sd) {
      double n = 0.0;
      for (size_t k = 1; k < sd.time_series_.size(); ++k)
        n += 0.5 * (sd.molar_flow_[k] + sd.molar_flow_[k - 1]) * (sd.time_series_[k] - sd.time_series_[k - 1]);
      return n;  // along the bed axis, positive upwards
    };
    // Donor at the bottom: gas leaving the bed flows downwards (negative).
    const double donor_axial = moles(donor->stream_data_);
    // Receiver at the bottom: gas entering the bed flows upwards (positive).
    const double receiver_axial = moles(rinse_in->stream_data_);
    std::cout << std::setprecision(8) << "donor axial moles " << donor_axial << ", receiver axial moles " << receiver_axial
              << ", donor logged location " << (donor->stream_data_.location_ == BoundaryLocation::TOP ? "TOP" : "BOTTOM") << "\n";
    ASSERT_LT(donor_axial, -0.1) << "The blowdown must release gas";
    EXPECT_EQ(donor->stream_data_.location_, BoundaryLocation::BOTTOM)
      << "The pressure boundary's log does not record the bed end it was taken at";
    EXPECT_NEAR(receiver_axial, -donor_axial, 1e-3 * std::abs(donor_axial))
      << "What left the donor's bed must enter the receiver's";
  }

  TEST(ConnectionSign, BottomPressureBoundaryDonorIsReplayedAsInflow) {
    ExpectRinseTakesInTheBlowdown(true);
  }

  // A step-change blowdown (constant 1 bar outlet) releases a large pulse at
  // the start. With the wrong sign the rinse inlet was asked to extract that
  // pulse from its bed and IDA failed at the start of the rinse.
  TEST(ConnectionSign, StepChangeBlowdownDonorIsReplayedAsInflow) {
    ExpectRinseTakesInTheBlowdown(false);
  }

} // namespace adsim
