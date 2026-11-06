/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
// SaveAdsorbentToJSON reports failures: a file it cannot write throws
// std::runtime_error naming the file and the reason, and leaves no temporary
// file behind. A successful save replaces the target and leaves only it.

#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

#include <gtest/gtest.h>

#include "adsim/objects/porous_media.h"

namespace adsim {
namespace {

namespace fs = std::filesystem;

class SaveAdsorbentErrors : public ::testing::Test {
 protected:
  void SetUp() override {
    dir_ = fs::temp_directory_path() /
           ("adsim_save_adsorbent_" + std::string(::testing::UnitTest::GetInstance()->current_test_info()->name()));
    fs::remove_all(dir_);
    fs::create_directories(dir_);
    adsorbent_.name_ = "Test sorbent";
    adsorbent_.solid_density_ = 1234.0;
  }
  void TearDown() override {
    std::error_code ec;
    fs::remove_all(dir_, ec);
  }
  int FilesIn(const fs::path& d) const {
    int n = 0;
    for ([[maybe_unused]] const auto& e : fs::directory_iterator(d)) ++n;
    return n;
  }

  fs::path dir_;
  Adsorbent adsorbent_;
};

TEST_F(SaveAdsorbentErrors, MissingDirectoryThrowsWithThePath) {
  const fs::path target = dir_ / "no_such_dir" / "sorbent.json";
  try {
    SaveAdsorbentToJSON(target.string(), adsorbent_);
    FAIL() << "expected std::runtime_error";
  } catch (const std::runtime_error& e) {
    const std::string what = e.what();
    EXPECT_NE(what.find(target.string()), std::string::npos) << what;
    EXPECT_NE(what.find("cannot create the temporary file"), std::string::npos) << what;
  }
  EXPECT_FALSE(fs::exists(target));
  EXPECT_EQ(FilesIn(dir_), 0);  // nothing left behind
}

TEST_F(SaveAdsorbentErrors, SuccessReplacesTheTargetAndLeavesNoTemporaryFile) {
  const fs::path target = dir_ / "sorbent.json";
  {
    std::ofstream old(target);
    old << "old contents";
  }
  ASSERT_NO_THROW(SaveAdsorbentToJSON(target.string(), adsorbent_));
  EXPECT_EQ(FilesIn(dir_), 1);
  Adsorbent back;
  LoadAdsorbentFromJSON(target.string(), back);
  EXPECT_EQ(back.name_, "Test sorbent");
  EXPECT_DOUBLE_EQ(back.solid_density_.Get(), 1234.0);
}

}  // namespace
}  // namespace adsim
