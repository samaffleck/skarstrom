// STL includes
#include <filesystem>

// AdSim includes
#include "adsim/solver/cycle.h"
#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/boundary_conditions/boundary_value.h"
#include "adsim/boundary_conditions/flow_boundary_value.h"
#include "adsim/boundary_conditions/connection.h"

// GTest includes
#include <gtest/gtest.h>


namespace adsim {

  const std::filesystem::path cycles_path = std::filesystem::current_path() / ".." / "data" / "cycles";

  // Generates the cycle templates in data/cycles/<name>.json (front ends offer
  // them as starting points) and round-trips each one. Keep `name` shell-safe
  // (lowercase, [a-z0-9_-]); see cmake/check_data_names.cmake.
  void SaveCycle(const std::string& name, const Cycle& cycle) {
    std::filesystem::create_directories(cycles_path);
    std::filesystem::path json_path = cycles_path / (name + ".json");
    SaveCycleToJSON(json_path.string(), cycle);
  }

  void VerifyCycleRoundTrip(const std::string& name, const Cycle& original) {
    std::filesystem::path json_path = cycles_path / (name + ".json");

    Cycle loaded;
    LoadCycleFromJSON(json_path.string(), loaded);

    EXPECT_EQ(loaded.number_of_cycles_, original.number_of_cycles_);
    EXPECT_EQ(loaded.step_list_.size(), original.step_list_.size());

    for (size_t i = 0; i < loaded.step_list_.size(); ++i) {
      EXPECT_EQ(loaded.step_list_[i].step_name_, original.step_list_[i].step_name_);
      EXPECT_DOUBLE_EQ(loaded.step_list_[i].step_time_, original.step_list_[i].step_time_);
    }
  }

  // ---------------------------------------------------------------------------
  // 4-Step VPSA: Feed, Blowdown, Evacuation, Pressurisation
  // Based on a typical vacuum pressure swing adsorption cycle for CO2/N2
  // ---------------------------------------------------------------------------
  TEST(CycleDesignTest, FourStepVPSA) {
    Cycle cycle;
    cycle.number_of_cycles_ = 5;

    // Feed step: flow inlet at bottom, pressure outlet at top
    cycle.AddStep("Feed");
    Step& feed = cycle.GetStep("Feed");
    feed.step_time_ = 60.0;

    auto feed_flow = std::make_unique<FlowBoundaryCondition>();
    feed_flow->location_ = BoundaryLocation::BOTTOM;
    feed_flow->direction_ = BoundaryDirection::INLET;
    auto feed_stream = std::make_unique<BoundaryStream>();
    feed_stream->name_ = "Feed inlet";
    auto vel = std::make_unique<VelocityBoundaryValue>();
    static_cast<ConstantBoundaryValue*>(vel->profile_.get())->constant_value_ = 0.5;
    feed_stream->boundary_flow_ = std::move(vel);
    feed_stream->temperature_ = 298.0;
    feed_stream->mole_fraction["CarbonDioxide"] = 0.15;
    feed_stream->mole_fraction["Nitrogen"] = 0.85;
    feed_flow->streams_.push_back(std::move(feed_stream));
    feed.bottom_boundary_condition_ = std::move(feed_flow);

    auto feed_pout = std::make_unique<PressureBoundaryCondition>();
    feed_pout->location_ = BoundaryLocation::TOP;
    feed_pout->direction_ = BoundaryDirection::OUTLET;
    auto feed_p = std::make_unique<ConstantBoundaryValue>();
    feed_p->constant_value_ = 1.0e5;
    feed_pout->pressure_ = std::move(feed_p);
    feed.top_boundary_condition_ = std::move(feed_pout);

    // Blowdown step: pressure outlet at top, closed bottom
    cycle.AddStep("Blowdown");
    Step& blowdown = cycle.GetStep("Blowdown");
    blowdown.step_time_ = 60.0;

    auto bd_pout = std::make_unique<PressureBoundaryCondition>();
    bd_pout->location_ = BoundaryLocation::TOP;
    bd_pout->direction_ = BoundaryDirection::OUTLET;
    auto bd_p = std::make_unique<ExponentialBoundaryValue>();
    bd_p->end_value_ = 0.1e5;
    bd_p->step_time_ = 60.0;
    bd_pout->pressure_ = std::move(bd_p);
    blowdown.top_boundary_condition_ = std::move(bd_pout);

    auto bd_closed = std::make_unique<ClosedBoundaryCondition>();
    bd_closed->location_ = BoundaryLocation::BOTTOM;
    blowdown.bottom_boundary_condition_ = std::move(bd_closed);

    // Evacuation step: closed top, pressure outlet at bottom
    cycle.AddStep("Evacuation");
    Step& evac = cycle.GetStep("Evacuation");
    evac.step_time_ = 200.0;

    auto evac_closed = std::make_unique<ClosedBoundaryCondition>();
    evac_closed->location_ = BoundaryLocation::TOP;
    evac.top_boundary_condition_ = std::move(evac_closed);

    auto evac_pout = std::make_unique<PressureBoundaryCondition>();
    evac_pout->location_ = BoundaryLocation::BOTTOM;
    evac_pout->direction_ = BoundaryDirection::OUTLET;
    auto evac_p = std::make_unique<ExponentialBoundaryValue>();
    evac_p->end_value_ = 0.01e5;
    evac_p->step_time_ = 200.0;
    evac_pout->pressure_ = std::move(evac_p);
    evac.bottom_boundary_condition_ = std::move(evac_pout);

    // Pressurisation step: pressure inlet at bottom, closed top
    cycle.AddStep("Pressurisation");
    Step& press = cycle.GetStep("Pressurisation");
    press.step_time_ = 30.0;

    auto press_pin = std::make_unique<PressureBoundaryCondition>();
    press_pin->location_ = BoundaryLocation::BOTTOM;
    press_pin->direction_ = BoundaryDirection::INLET;
    auto press_p = std::make_unique<LinearBoundaryValue>();
    press_p->start_value_ = 0.01e5;
    press_p->end_value_ = 1.0e5;
    press_p->ramp_time_ = 30.0;
    press_pin->pressure_ = std::move(press_p);
    press_pin->temperature_ = 298.0;
    press_pin->mole_fraction_["CarbonDioxide"] = 0.15;
    press_pin->mole_fraction_["Nitrogen"] = 0.85;
    press.bottom_boundary_condition_ = std::move(press_pin);

    auto press_closed = std::make_unique<ClosedBoundaryCondition>();
    press_closed->location_ = BoundaryLocation::TOP;
    press.top_boundary_condition_ = std::move(press_closed);

    SaveCycle("4_step_vpsa", cycle);
    VerifyCycleRoundTrip("4_step_vpsa", cycle);
  }

  // ---------------------------------------------------------------------------
  // Skarstrom Cycle: Feed, Blowdown, Purge, Pressurisation
  // Classic 4-step PSA cycle with a light product purge step
  // ---------------------------------------------------------------------------
  TEST(CycleDesignTest, SkarstromCycle) {
    Cycle cycle;
    cycle.number_of_cycles_ = 10;

    // Feed step: flow inlet at bottom, pressure outlet at top
    cycle.AddStep("Feed");
    Step& feed = cycle.GetStep("Feed");
    feed.step_time_ = 120.0;

    auto feed_flow = std::make_unique<FlowBoundaryCondition>();
    feed_flow->location_ = BoundaryLocation::BOTTOM;
    feed_flow->direction_ = BoundaryDirection::INLET;
    auto feed_stream = std::make_unique<BoundaryStream>();
    feed_stream->name_ = "Feed inlet";
    auto vel = std::make_unique<VelocityBoundaryValue>();
    static_cast<ConstantBoundaryValue*>(vel->profile_.get())->constant_value_ = 0.3;
    feed_stream->boundary_flow_ = std::move(vel);
    feed_stream->temperature_ = 298.0;
    feed_stream->mole_fraction["CarbonDioxide"] = 0.15;
    feed_stream->mole_fraction["Nitrogen"] = 0.85;
    feed_flow->streams_.push_back(std::move(feed_stream));
    feed.bottom_boundary_condition_ = std::move(feed_flow);

    auto feed_pout = std::make_unique<PressureBoundaryCondition>();
    feed_pout->location_ = BoundaryLocation::TOP;
    feed_pout->direction_ = BoundaryDirection::OUTLET;
    auto feed_p = std::make_unique<ConstantBoundaryValue>();
    feed_p->constant_value_ = 3.0e5;
    feed_pout->pressure_ = std::move(feed_p);
    feed.top_boundary_condition_ = std::move(feed_pout);

    // Blowdown step: pressure outlet at bottom, closed top
    cycle.AddStep("Blowdown");
    Step& blowdown = cycle.GetStep("Blowdown");
    blowdown.step_time_ = 30.0;

    auto bd_closed = std::make_unique<ClosedBoundaryCondition>();
    bd_closed->location_ = BoundaryLocation::TOP;
    blowdown.top_boundary_condition_ = std::move(bd_closed);

    auto bd_pout = std::make_unique<PressureBoundaryCondition>();
    bd_pout->location_ = BoundaryLocation::BOTTOM;
    bd_pout->direction_ = BoundaryDirection::OUTLET;
    auto bd_p = std::make_unique<ExponentialBoundaryValue>();
    bd_p->end_value_ = 1.0e5;
    bd_p->step_time_ = 30.0;
    bd_pout->pressure_ = std::move(bd_p);
    blowdown.bottom_boundary_condition_ = std::move(bd_pout);

    // Purge step: flow inlet at top (light product purge), pressure outlet at bottom
    cycle.AddStep("Purge");
    Step& purge = cycle.GetStep("Purge");
    purge.step_time_ = 120.0;

    auto purge_flow = std::make_unique<FlowBoundaryCondition>();
    purge_flow->location_ = BoundaryLocation::TOP;
    purge_flow->direction_ = BoundaryDirection::INLET;
    auto purge_stream = std::make_unique<BoundaryStream>();
    purge_stream->name_ = "Purge inlet";
    auto purge_vel = std::make_unique<VelocityBoundaryValue>();
    static_cast<ConstantBoundaryValue*>(purge_vel->profile_.get())->constant_value_ = 0.2;
    purge_stream->boundary_flow_ = std::move(purge_vel);
    purge_stream->temperature_ = 298.0;
    purge_stream->mole_fraction["CarbonDioxide"] = 0.0;
    purge_stream->mole_fraction["Nitrogen"] = 1.0;
    purge_flow->streams_.push_back(std::move(purge_stream));
    purge.top_boundary_condition_ = std::move(purge_flow);

    auto purge_pout = std::make_unique<PressureBoundaryCondition>();
    purge_pout->location_ = BoundaryLocation::BOTTOM;
    purge_pout->direction_ = BoundaryDirection::OUTLET;
    auto purge_p = std::make_unique<ConstantBoundaryValue>();
    purge_p->constant_value_ = 1.0e5;
    purge_pout->pressure_ = std::move(purge_p);
    purge.bottom_boundary_condition_ = std::move(purge_pout);

    // Pressurisation step: pressure inlet at bottom, closed top
    cycle.AddStep("Pressurisation");
    Step& press = cycle.GetStep("Pressurisation");
    press.step_time_ = 30.0;

    auto press_pin = std::make_unique<PressureBoundaryCondition>();
    press_pin->location_ = BoundaryLocation::BOTTOM;
    press_pin->direction_ = BoundaryDirection::INLET;
    auto press_p = std::make_unique<LinearBoundaryValue>();
    press_p->start_value_ = 1.0e5;
    press_p->end_value_ = 3.0e5;
    press_p->ramp_time_ = 30.0;
    press_pin->pressure_ = std::move(press_p);
    press_pin->temperature_ = 298.0;
    press_pin->mole_fraction_["CarbonDioxide"] = 0.15;
    press_pin->mole_fraction_["Nitrogen"] = 0.85;
    press.bottom_boundary_condition_ = std::move(press_pin);

    auto press_closed = std::make_unique<ClosedBoundaryCondition>();
    press_closed->location_ = BoundaryLocation::TOP;
    press.top_boundary_condition_ = std::move(press_closed);

    SaveCycle("skarstrom_cycle", cycle);
    VerifyCycleRoundTrip("skarstrom_cycle", cycle);
  }

  // ---------------------------------------------------------------------------
  // 4-Step VPSA with LPP: Feed, Blowdown, Evacuation, Light Product Pressurisation
  // Pressurisation uses light product (N2) from the top instead of feed from bottom
  // ---------------------------------------------------------------------------
  TEST(CycleDesignTest, FourStepVPSA_LPP) {
    Cycle cycle;
    cycle.number_of_cycles_ = 5;

    // Feed step
    cycle.AddStep("Feed");
    Step& feed = cycle.GetStep("Feed");
    feed.step_time_ = 60.0;

    auto feed_flow = std::make_unique<FlowBoundaryCondition>();
    feed_flow->location_ = BoundaryLocation::BOTTOM;
    feed_flow->direction_ = BoundaryDirection::INLET;
    auto feed_stream = std::make_unique<BoundaryStream>();
    feed_stream->name_ = "Feed inlet";
    auto vel = std::make_unique<VelocityBoundaryValue>();
    static_cast<ConstantBoundaryValue*>(vel->profile_.get())->constant_value_ = 0.5;
    feed_stream->boundary_flow_ = std::move(vel);
    feed_stream->temperature_ = 298.0;
    feed_stream->mole_fraction["CarbonDioxide"] = 0.15;
    feed_stream->mole_fraction["Nitrogen"] = 0.85;
    feed_flow->streams_.push_back(std::move(feed_stream));
    feed.bottom_boundary_condition_ = std::move(feed_flow);

    auto feed_pout = std::make_unique<PressureBoundaryCondition>();
    feed_pout->location_ = BoundaryLocation::TOP;
    feed_pout->direction_ = BoundaryDirection::OUTLET;
    auto feed_p = std::make_unique<ConstantBoundaryValue>();
    feed_p->constant_value_ = 1.0e5;
    feed_pout->pressure_ = std::move(feed_p);
    feed.top_boundary_condition_ = std::move(feed_pout);

    // Blowdown step
    cycle.AddStep("Blowdown");
    Step& blowdown = cycle.GetStep("Blowdown");
    blowdown.step_time_ = 60.0;

    auto bd_pout = std::make_unique<PressureBoundaryCondition>();
    bd_pout->location_ = BoundaryLocation::TOP;
    bd_pout->direction_ = BoundaryDirection::OUTLET;
    auto bd_p = std::make_unique<ExponentialBoundaryValue>();
    bd_p->end_value_ = 0.1e5;
    bd_p->step_time_ = 60.0;
    bd_pout->pressure_ = std::move(bd_p);
    blowdown.top_boundary_condition_ = std::move(bd_pout);

    auto bd_closed = std::make_unique<ClosedBoundaryCondition>();
    bd_closed->location_ = BoundaryLocation::BOTTOM;
    blowdown.bottom_boundary_condition_ = std::move(bd_closed);

    // Evacuation step
    cycle.AddStep("Evacuation");
    Step& evac = cycle.GetStep("Evacuation");
    evac.step_time_ = 200.0;

    auto evac_closed = std::make_unique<ClosedBoundaryCondition>();
    evac_closed->location_ = BoundaryLocation::TOP;
    evac.top_boundary_condition_ = std::move(evac_closed);

    auto evac_pout = std::make_unique<PressureBoundaryCondition>();
    evac_pout->location_ = BoundaryLocation::BOTTOM;
    evac_pout->direction_ = BoundaryDirection::OUTLET;
    auto evac_p = std::make_unique<ExponentialBoundaryValue>();
    evac_p->end_value_ = 0.01e5;
    evac_p->step_time_ = 200.0;
    evac_pout->pressure_ = std::move(evac_p);
    evac.bottom_boundary_condition_ = std::move(evac_pout);

    // Light Product Pressurisation: light product (N2) pushed in from top
    cycle.AddStep("LPP");
    Step& lpp = cycle.GetStep("LPP");
    lpp.step_time_ = 30.0;

    auto lpp_flow = std::make_unique<FlowBoundaryCondition>();
    lpp_flow->location_ = BoundaryLocation::TOP;
    lpp_flow->direction_ = BoundaryDirection::INLET;
    auto lpp_stream = std::make_unique<BoundaryStream>();
    lpp_stream->name_ = "Light product";
    auto lpp_vel = std::make_unique<VelocityBoundaryValue>();
    static_cast<ConstantBoundaryValue*>(lpp_vel->profile_.get())->constant_value_ = 0.3;
    lpp_stream->boundary_flow_ = std::move(lpp_vel);
    lpp_stream->temperature_ = 298.0;
    lpp_stream->mole_fraction["CarbonDioxide"] = 0.0;
    lpp_stream->mole_fraction["Nitrogen"] = 1.0;
    lpp_flow->streams_.push_back(std::move(lpp_stream));
    lpp.top_boundary_condition_ = std::move(lpp_flow);

    auto lpp_closed = std::make_unique<ClosedBoundaryCondition>();
    lpp_closed->location_ = BoundaryLocation::BOTTOM;
    lpp.bottom_boundary_condition_ = std::move(lpp_closed);

    SaveCycle("4_step_vpsa_lpp", cycle);
    VerifyCycleRoundTrip("4_step_vpsa_lpp", cycle);
  }

  // ---------------------------------------------------------------------------
  // 6-Step VPSA with Rinse: Feed, Rinse, Blowdown, Evacuation, LPP, Pressurisation
  // Includes a heavy reflux (rinse) step for improved CO2 purity
  // ---------------------------------------------------------------------------
  TEST(CycleDesignTest, SixStepVPSA_Rinse) {
    Cycle cycle;
    cycle.number_of_cycles_ = 5;

    // Feed step
    cycle.AddStep("Feed");
    Step& feed = cycle.GetStep("Feed");
    feed.step_time_ = 60.0;

    auto feed_flow = std::make_unique<FlowBoundaryCondition>();
    feed_flow->location_ = BoundaryLocation::BOTTOM;
    feed_flow->direction_ = BoundaryDirection::INLET;
    auto feed_stream = std::make_unique<BoundaryStream>();
    feed_stream->name_ = "Feed inlet";
    auto vel = std::make_unique<VelocityBoundaryValue>();
    static_cast<ConstantBoundaryValue*>(vel->profile_.get())->constant_value_ = 0.5;
    feed_stream->boundary_flow_ = std::move(vel);
    feed_stream->temperature_ = 298.0;
    feed_stream->mole_fraction["CarbonDioxide"] = 0.15;
    feed_stream->mole_fraction["Nitrogen"] = 0.85;
    feed_flow->streams_.push_back(std::move(feed_stream));
    feed.bottom_boundary_condition_ = std::move(feed_flow);

    auto feed_pout = std::make_unique<PressureBoundaryCondition>();
    feed_pout->location_ = BoundaryLocation::TOP;
    feed_pout->direction_ = BoundaryDirection::OUTLET;
    auto feed_p = std::make_unique<ConstantBoundaryValue>();
    feed_p->constant_value_ = 1.5e5;
    feed_pout->pressure_ = std::move(feed_p);
    feed.top_boundary_condition_ = std::move(feed_pout);

    // Rinse step: CO2-rich stream pushed in from bottom, outlet at top
    cycle.AddStep("Rinse");
    Step& rinse = cycle.GetStep("Rinse");
    rinse.step_time_ = 40.0;

    auto rinse_flow = std::make_unique<FlowBoundaryCondition>();
    rinse_flow->location_ = BoundaryLocation::BOTTOM;
    rinse_flow->direction_ = BoundaryDirection::INLET;
    auto rinse_stream = std::make_unique<BoundaryStream>();
    rinse_stream->name_ = "Rinse inlet";
    auto rinse_vel = std::make_unique<VelocityBoundaryValue>();
    static_cast<ConstantBoundaryValue*>(rinse_vel->profile_.get())->constant_value_ = 0.3;
    rinse_stream->boundary_flow_ = std::move(rinse_vel);
    rinse_stream->temperature_ = 298.0;
    rinse_stream->mole_fraction["CarbonDioxide"] = 0.95;
    rinse_stream->mole_fraction["Nitrogen"] = 0.05;
    rinse_flow->streams_.push_back(std::move(rinse_stream));
    rinse.bottom_boundary_condition_ = std::move(rinse_flow);

    auto rinse_pout = std::make_unique<PressureBoundaryCondition>();
    rinse_pout->location_ = BoundaryLocation::TOP;
    rinse_pout->direction_ = BoundaryDirection::OUTLET;
    auto rinse_p = std::make_unique<ConstantBoundaryValue>();
    rinse_p->constant_value_ = 1.5e5;
    rinse_pout->pressure_ = std::move(rinse_p);
    rinse.top_boundary_condition_ = std::move(rinse_pout);

    // Blowdown step
    cycle.AddStep("Blowdown");
    Step& blowdown = cycle.GetStep("Blowdown");
    blowdown.step_time_ = 60.0;

    auto bd_pout = std::make_unique<PressureBoundaryCondition>();
    bd_pout->location_ = BoundaryLocation::TOP;
    bd_pout->direction_ = BoundaryDirection::OUTLET;
    auto bd_p = std::make_unique<ExponentialBoundaryValue>();
    bd_p->end_value_ = 0.1e5;
    bd_p->step_time_ = 60.0;
    bd_pout->pressure_ = std::move(bd_p);
    blowdown.top_boundary_condition_ = std::move(bd_pout);

    auto bd_closed = std::make_unique<ClosedBoundaryCondition>();
    bd_closed->location_ = BoundaryLocation::BOTTOM;
    blowdown.bottom_boundary_condition_ = std::move(bd_closed);

    // Evacuation step
    cycle.AddStep("Evacuation");
    Step& evac = cycle.GetStep("Evacuation");
    evac.step_time_ = 200.0;

    auto evac_closed = std::make_unique<ClosedBoundaryCondition>();
    evac_closed->location_ = BoundaryLocation::TOP;
    evac.top_boundary_condition_ = std::move(evac_closed);

    auto evac_pout = std::make_unique<PressureBoundaryCondition>();
    evac_pout->location_ = BoundaryLocation::BOTTOM;
    evac_pout->direction_ = BoundaryDirection::OUTLET;
    auto evac_p = std::make_unique<ExponentialBoundaryValue>();
    evac_p->end_value_ = 0.01e5;
    evac_p->step_time_ = 200.0;
    evac_pout->pressure_ = std::move(evac_p);
    evac.bottom_boundary_condition_ = std::move(evac_pout);

    // LPP step
    cycle.AddStep("LPP");
    Step& lpp = cycle.GetStep("LPP");
    lpp.step_time_ = 20.0;

    auto lpp_flow = std::make_unique<FlowBoundaryCondition>();
    lpp_flow->location_ = BoundaryLocation::TOP;
    lpp_flow->direction_ = BoundaryDirection::INLET;
    auto lpp_stream = std::make_unique<BoundaryStream>();
    lpp_stream->name_ = "Light product";
    auto lpp_vel = std::make_unique<VelocityBoundaryValue>();
    static_cast<ConstantBoundaryValue*>(lpp_vel->profile_.get())->constant_value_ = 0.2;
    lpp_stream->boundary_flow_ = std::move(lpp_vel);
    lpp_stream->temperature_ = 298.0;
    lpp_stream->mole_fraction["CarbonDioxide"] = 0.0;
    lpp_stream->mole_fraction["Nitrogen"] = 1.0;
    lpp_flow->streams_.push_back(std::move(lpp_stream));
    lpp.top_boundary_condition_ = std::move(lpp_flow);

    auto lpp_closed = std::make_unique<ClosedBoundaryCondition>();
    lpp_closed->location_ = BoundaryLocation::BOTTOM;
    lpp.bottom_boundary_condition_ = std::move(lpp_closed);

    // Pressurisation step
    cycle.AddStep("Pressurisation");
    Step& press = cycle.GetStep("Pressurisation");
    press.step_time_ = 30.0;

    auto press_pin = std::make_unique<PressureBoundaryCondition>();
    press_pin->location_ = BoundaryLocation::BOTTOM;
    press_pin->direction_ = BoundaryDirection::INLET;
    auto press_p = std::make_unique<LinearBoundaryValue>();
    press_p->start_value_ = 0.01e5;
    press_p->end_value_ = 1.5e5;
    press_p->ramp_time_ = 30.0;
    press_pin->pressure_ = std::move(press_p);
    press_pin->temperature_ = 298.0;
    press_pin->mole_fraction_["CarbonDioxide"] = 0.15;
    press_pin->mole_fraction_["Nitrogen"] = 0.85;
    press.bottom_boundary_condition_ = std::move(press_pin);

    auto press_closed = std::make_unique<ClosedBoundaryCondition>();
    press_closed->location_ = BoundaryLocation::TOP;
    press.top_boundary_condition_ = std::move(press_closed);

    SaveCycle("6_step_vpsa_rinse", cycle);
    VerifyCycleRoundTrip("6_step_vpsa_rinse", cycle);
  }

  // ---------------------------------------------------------------------------
  // Simple Breakthrough: Single feed step (useful for breakthrough curve studies)
  // ---------------------------------------------------------------------------
  TEST(CycleDesignTest, Breakthrough) {
    Cycle cycle;
    cycle.number_of_cycles_ = 1;

    cycle.AddStep("Feed");
    Step& feed = cycle.GetStep("Feed");
    feed.step_time_ = 600.0;

    auto feed_flow = std::make_unique<FlowBoundaryCondition>();
    feed_flow->location_ = BoundaryLocation::BOTTOM;
    feed_flow->direction_ = BoundaryDirection::INLET;
    auto feed_stream = std::make_unique<BoundaryStream>();
    feed_stream->name_ = "Feed inlet";
    auto vel = std::make_unique<VelocityBoundaryValue>();
    static_cast<ConstantBoundaryValue*>(vel->profile_.get())->constant_value_ = 0.3;
    feed_stream->boundary_flow_ = std::move(vel);
    feed_stream->temperature_ = 298.0;
    feed_stream->mole_fraction["CarbonDioxide"] = 0.15;
    feed_stream->mole_fraction["Nitrogen"] = 0.85;
    feed_flow->streams_.push_back(std::move(feed_stream));
    feed.bottom_boundary_condition_ = std::move(feed_flow);

    auto feed_pout = std::make_unique<PressureBoundaryCondition>();
    feed_pout->location_ = BoundaryLocation::TOP;
    feed_pout->direction_ = BoundaryDirection::OUTLET;
    auto feed_p = std::make_unique<ConstantBoundaryValue>();
    feed_p->constant_value_ = 1.0e5;
    feed_pout->pressure_ = std::move(feed_p);
    feed.top_boundary_condition_ = std::move(feed_pout);

    SaveCycle("breakthrough", cycle);
    VerifyCycleRoundTrip("breakthrough", cycle);
  }

  // ---------------------------------------------------------------------------
  // DuplicateStep
  // ---------------------------------------------------------------------------

  // Build a step with a pressure BC at the top and a flow BC (one stream) at
  // the bottom — exercises both polymorphic BC branches in DuplicateStep.
  static void BuildSampleStep(Cycle& cycle, const std::string& name, double step_time,
                              double pressure_pa, double velocity) {
    cycle.AddStep(name);
    Step& s = cycle.GetStep(name);
    s.step_time_ = step_time;

    auto p_bc = std::make_unique<PressureBoundaryCondition>();
    p_bc->location_ = BoundaryLocation::TOP;
    p_bc->direction_ = BoundaryDirection::OUTLET;
    auto p_val = std::make_unique<ConstantBoundaryValue>();
    p_val->constant_value_ = pressure_pa;
    p_bc->pressure_ = std::move(p_val);
    s.top_boundary_condition_ = std::move(p_bc);

    auto flow = std::make_unique<FlowBoundaryCondition>();
    flow->location_ = BoundaryLocation::BOTTOM;
    flow->direction_ = BoundaryDirection::INLET;
    auto stream = std::make_unique<BoundaryStream>();
    stream->name_ = "Inlet";
    auto vel = std::make_unique<VelocityBoundaryValue>();
    static_cast<ConstantBoundaryValue*>(vel->profile_.get())->constant_value_ = velocity;
    stream->boundary_flow_ = std::move(vel);
    stream->temperature_ = 298.0;
    flow->streams_.push_back(std::move(stream));
    s.bottom_boundary_condition_ = std::move(flow);
  }

  TEST(CycleDuplicateStep, InsertsCopyImmediatelyAfterOriginal) {
    Cycle cycle;
    BuildSampleStep(cycle, "Feed",   60.0, 1.0e5, 0.5);
    BuildSampleStep(cycle, "Purge", 30.0, 1.0e4, 0.2);
    ASSERT_EQ(cycle.step_list_.size(), 2u);

    size_t new_idx = cycle.DuplicateStep(0);
    EXPECT_EQ(new_idx, 1u);
    ASSERT_EQ(cycle.step_list_.size(), 3u);
    EXPECT_EQ(cycle.step_list_[0].step_name_, "Feed");
    EXPECT_EQ(cycle.step_list_[1].step_name_, "Feed (copy)");
    EXPECT_EQ(cycle.step_list_[2].step_name_, "Purge");
  }

  TEST(CycleDuplicateStep, OutOfRangeReturnsSentinel) {
    Cycle cycle;
    BuildSampleStep(cycle, "Feed", 60.0, 1.0e5, 0.5);

    size_t res = cycle.DuplicateStep(99);
    EXPECT_EQ(res, static_cast<size_t>(-1));
    EXPECT_EQ(cycle.step_list_.size(), 1u);
  }

  TEST(CycleDuplicateStep, MakesUniqueNameWhenCopyAlreadyExists) {
    Cycle cycle;
    BuildSampleStep(cycle, "Feed", 60.0, 1.0e5, 0.5);

    cycle.DuplicateStep(0);
    cycle.DuplicateStep(0);
    cycle.DuplicateStep(0);

    ASSERT_EQ(cycle.step_list_.size(), 4u);
    // Order: Feed, Feed (copy 3), Feed (copy 2), Feed (copy) — newer copies
    // are inserted right after the original, pushing older copies down.
    std::set<std::string> names;
    for (const auto& s : cycle.step_list_) names.insert(s.step_name_);
    EXPECT_EQ(names.size(), 4u);  // all four names are unique
    EXPECT_TRUE(names.count("Feed"));
    EXPECT_TRUE(names.count("Feed (copy)"));
    EXPECT_TRUE(names.count("Feed (copy 2)"));
    EXPECT_TRUE(names.count("Feed (copy 3)"));
  }

  TEST(CycleDuplicateStep, DeepCopiesBoundaryConditions) {
    Cycle cycle;
    BuildSampleStep(cycle, "Feed", 60.0, 1.0e5, 0.5);

    cycle.DuplicateStep(0);
    ASSERT_EQ(cycle.step_list_.size(), 2u);

    auto* orig_p = dynamic_cast<PressureBoundaryCondition*>(
        cycle.step_list_[0].top_boundary_condition_.get());
    auto* copy_p = dynamic_cast<PressureBoundaryCondition*>(
        cycle.step_list_[1].top_boundary_condition_.get());
    ASSERT_NE(orig_p, nullptr);
    ASSERT_NE(copy_p, nullptr);

    // Values match.
    auto* orig_val = dynamic_cast<ConstantBoundaryValue*>(orig_p->pressure_.get());
    auto* copy_val = dynamic_cast<ConstantBoundaryValue*>(copy_p->pressure_.get());
    ASSERT_NE(orig_val, nullptr);
    ASSERT_NE(copy_val, nullptr);
    EXPECT_DOUBLE_EQ(orig_val->constant_value_, copy_val->constant_value_);

    // But pointers are distinct — mutating the copy must not affect the original.
    EXPECT_NE(orig_p, copy_p);
    EXPECT_NE(orig_p->pressure_.get(), copy_p->pressure_.get());

    copy_val->constant_value_ = 99.0;
    EXPECT_DOUBLE_EQ(orig_val->constant_value_, 1.0e5);

    // Same independence on the flow side.
    auto* orig_f = dynamic_cast<FlowBoundaryCondition*>(
        cycle.step_list_[0].bottom_boundary_condition_.get());
    auto* copy_f = dynamic_cast<FlowBoundaryCondition*>(
        cycle.step_list_[1].bottom_boundary_condition_.get());
    ASSERT_NE(orig_f, nullptr);
    ASSERT_NE(copy_f, nullptr);
    EXPECT_NE(orig_f, copy_f);
    ASSERT_EQ(copy_f->streams_.size(), 1u);
    EXPECT_NE(orig_f->streams_[0].get(), copy_f->streams_[0].get());
  }

  TEST(CycleDuplicateStep, ConnectionsOnDuplicateAreCleared) {
    Cycle cycle;
    BuildSampleStep(cycle, "Feed",  60.0, 1.0e5, 0.5);
    BuildSampleStep(cycle, "Purge", 30.0, 1.0e4, 0.2);

    // Connect Feed's bottom flow stream to Purge's top pressure boundary —
    // a stream↔pressure pairing is the only kind MakeMutualConnection supports.
    auto* feed_flow = dynamic_cast<FlowBoundaryCondition*>(
        cycle.step_list_[0].bottom_boundary_condition_.get());
    auto* purge_p = dynamic_cast<PressureBoundaryCondition*>(
        cycle.step_list_[1].top_boundary_condition_.get());
    ASSERT_NE(feed_flow, nullptr);
    ASSERT_NE(purge_p, nullptr);
    ASSERT_EQ(feed_flow->streams_.size(), 1u);
    auto& feed_stream = *feed_flow->streams_[0];

    ConnectionEndpoint a{ &feed_stream.connection_, &cycle.step_list_[0], &feed_stream, nullptr };
    ConnectionEndpoint b{ &purge_p->connection_,    &cycle.step_list_[1], nullptr, purge_p };
    ASSERT_TRUE(MakeMutualConnection(a, b));
    ASSERT_TRUE(feed_stream.connection_.IsConnected());

    // Duplicate Feed. The duplicate lives at index 1, pushing Purge to index 2.
    cycle.DuplicateStep(0);
    ASSERT_EQ(cycle.step_list_.size(), 3u);

    auto* copy_flow = dynamic_cast<FlowBoundaryCondition*>(
        cycle.step_list_[1].bottom_boundary_condition_.get());
    ASSERT_NE(copy_flow, nullptr);
    ASSERT_EQ(copy_flow->streams_.size(), 1u);
    auto& copy_stream = *copy_flow->streams_[0];

    EXPECT_FALSE(copy_stream.connection_.IsConnected());
    EXPECT_EQ(copy_stream.connection_.GetID(), 0u);

    // Original connection unchanged.
    EXPECT_TRUE(feed_stream.connection_.IsConnected());
  }

} // namespace adsim
