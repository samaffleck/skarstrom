// The generated adsim/version.h must agree with itself: the version string is
// built from project(AdSim VERSION ...) in the root CMakeLists.txt, so the
// string and the numeric components can never drift apart.
#include <gtest/gtest.h>

#include <string>

#include "adsim/version.h"

TEST(Version, StringMatchesComponents) {
  const std::string expected = std::to_string(ADSIM_VERSION_MAJOR) + "." +
                               std::to_string(ADSIM_VERSION_MINOR) + "." +
                               std::to_string(ADSIM_VERSION_PATCH) + ADSIM_VERSION_SUFFIX;
  EXPECT_EQ(std::string(ADSIM_VERSION_STRING), expected);
  EXPECT_EQ(std::string(adsim::kVersionString), expected);
}

TEST(Version, IsSemanticVersion) {
  // SemVer MAJOR.MINOR.PATCH (the open-source releases restarted at 0.1.0).
  // Guards against reverting to the old calendar versions such as 2026.3.0.
  EXPECT_GE(adsim::kVersionMajor, 0);
  EXPECT_LT(adsim::kVersionMajor, 1000) << "looks like a calendar year, not a SemVer major";
  EXPECT_GE(adsim::kVersionMinor, 0);
  EXPECT_GE(adsim::kVersionPatch, 0);
}
