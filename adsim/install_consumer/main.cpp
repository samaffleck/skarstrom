// Downstream smoke test for the installed adsim package (see CMakeLists.txt).
//
// Includes headers that reach every third-party dependency of the public API
// (cereal, Eigen, SUNDIALS, CoolProp, dlib), loads a saved deck, runs it with
// the IDA solver and prints the outcome. Exit code 0 means the deck ran.

#include "adsim/version.h"
#include "adsim/objects/user_data.h"        // cereal, Eigen, CoolProp
#include "adsim/solver/coupled_solver.h"    // SUNDIALS IDA
#include "adsim/solver/optimisation.h"      // dlib

#include <cstdio>
#include <exception>

int main(int argc, char** argv) {
  std::printf("adsim %s (package consumer)\n", adsim::kVersionString);
  if (argc < 2) {
    std::fprintf(stderr, "usage: %s deck.json\n", argv[0]);
    return 2;
  }
  try {
    adsim::UserData data;
    adsim::LoadFromJson(data, argv[1]);
    adsim::CoupledSolver solver{data};
    const int rc = solver.IntegrateCycle(data.cycle, false);
    const auto stats = solver.GetStatistics();
    std::printf("deck=%s rc=%d steps=%ld residual_evals=%ld\n", argv[1], rc,
                stats.num_steps, stats.num_res_evals);
    return rc < 0 ? 1 : 0;
  } catch (const std::exception& e) {
    std::fprintf(stderr, "error: %s\n", e.what());
    return 1;
  }
}
