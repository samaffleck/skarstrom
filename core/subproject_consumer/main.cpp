// Front-end smoke test for a subproject build of the engine (see CMakeLists.txt).
//
// Touches the public API a front end uses: the Simulation facade, deck I/O,
// the structural deck edits, the parameter registry and the isotherm fitter.
// Runs one cycle of the deck given on the command line. Exit code 0 = ok.

#include <cstdio>
#include <exception>
#include <fstream>
#include <sstream>
#include <string>

#include "adsim/equations/isotherm_kinds.h"
#include "adsim/objects/deck_edit.h"
#include "adsim/objects/deck_io.h"
#include "adsim/simulation.h"
#include "adsim/solver/isotherm_fit.h"
#include "adsim/version.h"

int main(int argc, char** argv) {
  std::printf("adsim %s (subproject consumer)\n", adsim::kVersionString);
  if (argc < 2) {
    std::fprintf(stderr, "usage: %s deck.json\n", argv[0]);
    return 2;
  }
  try {
    std::ifstream in(argv[1], std::ios::binary);
    std::stringstream text;
    text << in.rdbuf();
    adsim::UserData data;
    adsim::LoadDeckFlexible(data, text.str());

    std::string err;
    if (!adsim::AddComponentToDeck(data, "Argon", err)) throw std::runtime_error(err);
    adsim::RemoveComponentFromDeck(data, "Argon");

    adsim::Simulation sim;
    sim.SetData(std::move(data));
    sim.SetParameter("cycle.number_of_cycles", 1);
    const adsim::RunResult r = sim.Run();
    std::printf("run: %s (%s), %d cycle(s)\n", adsim::RunStatusName(r.status), r.message.c_str(), r.cycles_run);

    std::vector<adsim::isofit::Point> pts;
    for (double p : {1e4, 5e4, 1e5, 3e5}) {
      adsim::isofit::Point pt;
      pt.P = p;
      pt.q = 4.0 * 2e-5 * p / (1.0 + 2e-5 * p);
      pt.T = 298.15;
      pts.push_back(pt);
    }
    adsim::IsothermType sips{};
    adsim::KindFromString("sips", sips);
    std::vector<adsim::isofit::ParamSpec> specs;
    double t_ref = 0.0;
    if (!adsim::isofit::DefaultSpecs(sips, pts, specs, t_ref, err)) throw std::runtime_error(err);
    const auto fit = adsim::isofit::Fit(sips, pts, specs, t_ref);
    std::printf("isotherm fit: %s, R2 = %.6f\n", fit.ok ? "ok" : fit.error.c_str(), fit.r2);
    return (r.Ok() && fit.ok) ? 0 : 1;
  } catch (const std::exception& e) {
    std::fprintf(stderr, "error: %s\n", e.what());
    return 1;
  }
}
