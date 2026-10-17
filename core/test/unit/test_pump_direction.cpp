/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Tests for EnsurePumpMatchesDirection, the helper the GUI boundary-condition
// menus call every frame to keep a BC's pump type consistent with its flow
// direction (INLET -> compressor, OUTLET -> vacuum pump). It must only
// recreate the pump on a genuine type mismatch; a matching pump keeps its
// user settings (regression for G-02, where a shared static "last direction"
// recreated pumps whenever BCs of different directions were shown in turn).

// STL includes
#include <memory>

// adsim includes
#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/objects/compressor.h"

// Gtest includes
#include <gtest/gtest.h>


namespace adsim {

  TEST(PumpDirectionTest, PumpTypeForDirection) {
    EXPECT_EQ(PumpTypeForDirection(BoundaryDirection::INLET), PumpType::COMPRESSOR);
    EXPECT_EQ(PumpTypeForDirection(BoundaryDirection::OUTLET), PumpType::VACUUM_PUMP);
  }

  TEST(PumpDirectionTest, NullPumpIsCreated) {
    std::unique_ptr<IPump> pump;
    EXPECT_TRUE(EnsurePumpMatchesDirection(pump, BoundaryDirection::OUTLET));
    ASSERT_NE(pump, nullptr);
    EXPECT_EQ(pump->GetType(), PumpType::VACUUM_PUMP);
  }

  TEST(PumpDirectionTest, MismatchedPumpIsReplaced) {
    std::unique_ptr<IPump> pump = CreatePump(PumpType::COMPRESSOR);
    pump->use_pump_ = true;
    EXPECT_TRUE(EnsurePumpMatchesDirection(pump, BoundaryDirection::OUTLET));
    EXPECT_EQ(pump->GetType(), PumpType::VACUUM_PUMP);
    EXPECT_FALSE(pump->use_pump_); // fresh default pump

    EXPECT_TRUE(EnsurePumpMatchesDirection(pump, BoundaryDirection::INLET));
    EXPECT_EQ(pump->GetType(), PumpType::COMPRESSOR);
  }

  TEST(PumpDirectionTest, MatchingPumpKeepsSettings) {
    std::unique_ptr<IPump> pump = CreatePump(PumpType::COMPRESSOR);
    pump->use_pump_ = true;
    pump->poly_index_ = 1.25;
    const IPump* before = pump.get();

    EXPECT_FALSE(EnsurePumpMatchesDirection(pump, BoundaryDirection::INLET));
    EXPECT_EQ(pump.get(), before);
    EXPECT_TRUE(pump->use_pump_);
    EXPECT_DOUBLE_EQ(static_cast<double>(pump->poly_index_), 1.25);
  }

  // G-02 scenario: an INLET pressure BC (step A) and an OUTLET pressure BC
  // (step B), each with the power calculation enabled, shown alternately as
  // the user browses steps. Neither pump may be recreated.
  TEST(PumpDirectionTest, BrowsingBoundaryConditionsOfOppositeDirectionKeepsPumps) {
    PressureBoundaryCondition bc_a;
    bc_a.direction_ = BoundaryDirection::INLET;
    EnsurePumpMatchesDirection(bc_a.pump_, bc_a.direction_);
    bc_a.pump_->use_pump_ = true;

    PressureBoundaryCondition bc_b;
    bc_b.direction_ = BoundaryDirection::OUTLET;
    EnsurePumpMatchesDirection(bc_b.pump_, bc_b.direction_);
    bc_b.pump_->use_pump_ = true;

    const IPump* pump_a = bc_a.pump_.get();
    const IPump* pump_b = bc_b.pump_.get();

    for (int frame = 0; frame < 4; ++frame) {
      EXPECT_FALSE(EnsurePumpMatchesDirection(bc_a.pump_, bc_a.direction_));
      EXPECT_FALSE(EnsurePumpMatchesDirection(bc_b.pump_, bc_b.direction_));
    }
    EXPECT_EQ(bc_a.pump_.get(), pump_a);
    EXPECT_EQ(bc_b.pump_.get(), pump_b);
    EXPECT_TRUE(bc_a.pump_->use_pump_);
    EXPECT_TRUE(bc_b.pump_->use_pump_);

    // A genuine direction change still swaps the pump type.
    bc_a.direction_ = BoundaryDirection::OUTLET;
    EXPECT_TRUE(EnsurePumpMatchesDirection(bc_a.pump_, bc_a.direction_));
    EXPECT_EQ(bc_a.pump_->GetType(), PumpType::VACUUM_PUMP);
  }

} // namespace adsim
