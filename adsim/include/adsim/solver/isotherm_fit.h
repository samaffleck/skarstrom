/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// Isotherm fitting from (p, q[, T]) points. Used by the Python package
// (skarstrom.fit_isotherm) and by any front end that offers isotherm fitting;
// the implementation is in adsim/src/solver/isotherm_fit.cpp.
//
// Every model evaluation goes through the ENGINE's isotherm code: a real
// IIsotherm object is built with CreateIsotherm, its parameters
// are addressed through IsothermFields (the single place that knows each
// kind's parameter layout) and q is computed with EvalIsotherm (the
// IsothermFunction free functions the solver itself uses). There is no
// second copy of any isotherm formula here or in any front end.
//
// Fitting = minimise SSE = sum_i (q_model(P_i, T_i, y_i) - q_i)^2 over the free
// parameters:
//   1. dlib::find_min_global (MaxLIPO + trust region) over a transformed box
//      (log10 for positive scale parameters, linear otherwise), so parameters
//      that span many decades (b ~ 1e-10 1/Pa next to Q ~ 3e4 J/mol) are
//      searched evenly;
//   2. a bounded Levenberg-Marquardt polish (box-projected, finite-difference
//      Jacobian) started from the best few distinct points of step 1.
// An Arrhenius pre-exponential (b in b*exp(Q/RT)) is searched as its value at
// the reference temperature T_ref, b_ref = b*exp(Q/(R*T_ref)); this removes the
// near-total b/Q correlation that otherwise makes the search ill-conditioned.
// Bounds of such a "paired" parameter therefore apply to b_ref.
//
// Goodness of fit: SSE, RMSE = sqrt(SSE/n), R^2 = 1 - SSE/SST. Parameter
// uncertainty is the usual linearised estimate: Cov = s^2 (J^T J)^-1 with
// s^2 = SSE/(n - p) and J the Jacobian of q w.r.t. the free parameters at the
// optimum (central differences). Reported as one standard error; it is only
// meaningful for a well-identified fit with roughly independent, homoscedastic
// errors, and is null when J^T J is (numerically) singular.

#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "adsim/equations/isotherm.h"
#include "adsim/equations/isotherm_kinds.h"

namespace adsim {
namespace isofit {

// One measured equilibrium point, SI base units.
struct Point {
  double P = 0.0;   // total pressure [Pa]
  double q = 0.0;   // loading [mol/kg]
  double T = 298.15;  // temperature [K]
  double y = 1.0;   // mole fraction (1 for pure-component data)
};

enum class Scale { Linear, Log };

// One isotherm parameter as the fitter sees it. Values and bounds are SI base
// (the units the deck stores). For a paired parameter (pair_energy non-empty)
// `lower`/`upper` bound the reference-temperature value b*exp(E/(R*T_ref)).
struct ParamSpec {
  std::string name;
  double value = 0.0;   // initial guess, or the held value when fixed
  double lower = 0.0;
  double upper = 0.0;
  bool fixed = false;
  Scale scale = Scale::Linear;
  std::string pair_energy;  // name of the energy parameter E in b*exp(E/RT)
  std::string note;         // why a default was chosen / fixed (UI hint)
};

struct ParamResult {
  std::string name;
  double value = 0.0;
  double std_error = std::numeric_limits<double>::quiet_NaN();  // NaN => n/a
  bool fixed = false;
};

struct FitResult {
  bool ok = false;
  std::string error;
  std::string warning;
  std::vector<ParamResult> params;
  double sse = 0.0;
  double rmse = 0.0;
  double r2 = 0.0;
  int n = 0;            // points used
  int dof = 0;          // n - free parameters
  int evaluations = 0;  // SSE evaluations spent
  double t_ref = 0.0;   // reference temperature used for paired parameters
  std::vector<double> q_fit;  // model at each input point
};

struct FitOptions {
  int max_global_evals = 0;  // 0 => automatic (100 + 80 * free parameters)
  int max_local_evals = 4000;
};

// ---------------------------------------------------------------------------
// Model: a single-component engine isotherm addressed by parameter name.
// ---------------------------------------------------------------------------

class Model {
 public:
  explicit Model(IsothermType type) : type_(type) {
    iso_ = CreateIsotherm(type, kComp);
    if (iso_) fields_ = IsothermFields(iso_.get(), kComp);
  }

  bool valid() const { return iso_ != nullptr; }
  IsothermType type() const { return type_; }

  std::vector<std::string> Names() const {
    std::vector<std::string> out;
    for (const auto& f : fields_) out.push_back(f.first);
    return out;
  }

  ExpressionParam* Field(const std::string& name) const {
    for (const auto& f : fields_) if (f.first == name) return f.second;
    return nullptr;
  }

  bool Set(const std::string& name, double v) {
    auto* p = Field(name);
    if (!p) return false;
    p->Get() = v;  // value only; the fitter never needs the expression text
    return true;
  }

  double Get(const std::string& name) const {
    auto* p = Field(name);
    return p ? p->Get() : std::numeric_limits<double>::quiet_NaN();
  }

  double Eval(double T, double P, double y) const {
    return EvalIsotherm(iso_.get(), kComp, T, P, y);
  }

 private:
  static constexpr const char* kComp = "fit";
  IsothermType type_;
  std::unique_ptr<IIsotherm> iso_;
  std::vector<std::pair<std::string, ExpressionParam*>> fields_;
};

// The kinds the tool offers. Each has data-driven default bounds below; any
// other kind can still be fitted if the caller supplies a full ParamSpec list.
const std::vector<IsothermType>& FittableKinds();

// True when the kind's loading depends on temperature through parameters the
// tool fits (so multi-temperature data can identify them).
bool IsTemperatureDependent(IsothermType t);

// ---------------------------------------------------------------------------
// Data summary used for default bounds.
// ---------------------------------------------------------------------------

struct DataSummary {
  int n = 0;
  double q_max = 0.0;
  double p_max = 0.0;       // max partial pressure [Pa]
  double p_min_pos = 0.0;   // min positive partial pressure [Pa]
  double c_max = 0.0;       // max gas concentration P*y/(R*T) [mol/m3]
  double c_min_pos = 0.0;
  double slope_min = 0.0;   // min q/p over points with p, q > 0
  double slope_max = 0.0;
  double t_mean = 0.0;
  int n_temperatures = 0;   // distinct temperatures (0.5 K tolerance)
};

bool PointUsable(const Point& p);

DataSummary Summarise(const std::vector<Point>& pts);

// Reference temperature: the mean data temperature, rounded to 0.01 K so it
// reads cleanly when written into a deck as T_ref / T_0.
double ReferenceTemperature(const DataSummary& s);

// ---------------------------------------------------------------------------
// Default parameter specs (initial guesses + search bounds) from the data.
// These are search-space heuristics (decades around what the data can
// resolve), not material properties.
// ---------------------------------------------------------------------------

bool DefaultSpecs(IsothermType type, const std::vector<Point>& pts,
                  std::vector<ParamSpec>& out, double& t_ref, std::string& err);

FitResult Fit(IsothermType type, const std::vector<Point>& points_in,
              std::vector<ParamSpec> specs, double t_ref = 0.0,
              const FitOptions& opt = {});

// Evaluate a kind at arbitrary points with the given parameter values (any
// parameter not listed keeps its engine default). Returns false on an unknown
// parameter name.
bool Evaluate(IsothermType type, const std::vector<std::pair<std::string, double>>& values,
              const std::vector<Point>& pts, std::vector<double>& q, std::string& err);

}  // namespace isofit
}  // namespace adsim
