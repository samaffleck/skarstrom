/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

#include "adsim/math/user_variables.h"

#include <string>
#include <cstdint>
#include <memory>
#include <vector>

#include <cereal/types/string.hpp>
#include <cereal/access.hpp>
#include <cereal/cereal.hpp>

namespace adsim {

class ExpressionParam;

// Collects the parameters whose cached value is non-finite (inf/NaN) during an
// audit walk (see EvalContext::audit and UserData::FindNonFiniteExpressions).
// The walker sets `scope` to a human-readable name before visiting each
// sub-tree; every recorded entry copies it.
struct NonFiniteExpressionAudit {
    struct Entry {
        const ExpressionParam* param = nullptr;
        std::string scope;  // where the parameter lives, e.g. "Bed geometry"
        std::string expr;   // the expression text
        double      value = 0.0;
    };
    std::string        scope;
    std::vector<Entry> entries;

    // Adds `p` unless it was already recorded (the walker may visit a
    // parameter twice: once under a specific scope, once under a coarser one).
    void Record(const ExpressionParam& p);
};

// Simulation state made available to expression evaluation.
// Fill whichever fields are meaningful for the parameter being evaluated.
struct EvalContext {
    const UserVariableRegistry* registry = nullptr;
    double T = 298.15;  // fluid temperature [K]
    double P = 1e5;     // pressure [Pa]
    double c = 40.0;    // total molar concentration [mol/m³]
    double y = 0.0;     // mole fraction [-]
    double v = 0.0;     // velocity [m/s]
    // Audit mode. When set, Evaluate() does NOT re-evaluate anything: it only
    // records the parameter in *audit if IsNonFiniteConstant(). This lets the
    // EvaluateExpressions cascade double as a read-only visitor over every
    // ExpressionParam the deck owns.
    NonFiniteExpressionAudit* audit = nullptr;
};

// Wraps a double value plus its expression text. The expression text is the
// single source of truth shown in the GUI; the cached double is what the
// solver consumes. Internal numeric assignments (operator=(double)) keep both
// in sync by re-stringifying the new value into expr_str_.
class ExpressionParam {
public:
    ExpressionParam() : value_(0.0), expr_str_(FormatScalar(0.0)) {}
    ExpressionParam(double init) : value_(init), expr_str_(FormatScalar(init)) {}

    // Assigning a plain number replaces expr_str_ with the formatted value so
    // the GUI display tracks the new value.
    ExpressionParam& operator=(double v) {
        value_ = v;
        expr_str_ = FormatScalar(v);
        last_error_.clear();
        has_inline_unit_ = false;
        return *this;
    }

    operator double()       { return value_; }
    operator double() const { return value_; }

    double&       Get()       { return value_; }
    double        Get() const { return value_; }

    // True when expr_str_ contains anything beyond a plain numeric literal
    // (so optimisation passes etc. can distinguish "user typed a formula" from
    // "this is just a number 1.5e5"). expr_str_ is always non-empty after
    // default construction, so emptiness is not a useful check.
    bool               HasExpression() const;
    const std::string& GetExpr()       const { return expr_str_; }
    const std::string& LastError()     const { return last_error_; }

    // True when the most recent Evaluate() consumed a trailing `[unit]` suffix
    // and matched a known symbol. False for plain numeric expressions or for
    // expressions without a bracketed unit.
    bool               HasInlineUnit()      const { return has_inline_unit_; }
    UnitCategory       InlineUnitCategory() const { return inline_unit_category_; }

    // Direct mutable access to the canonical expression text (a base-unit
    // number, or a formula / inline-[unit] expression the core parser
    // understands). value_ is always kept equal to Evaluate(expr_str_), so this
    // is what the solver consumes.
    std::string&       GetExprMutable()      { return expr_str_; }

    // Index into the GUI unit set (from GetUnits/GetUnitsForCategory) that this
    // field is currently displayed/entered in. Persisted so the user's unit
    // choice survives save/load. Lives here (not in transient GUI state) to
    // avoid keying per-field UI state by raw pointer. 0 == SI base unit.
    int&               DisplayUnitIndex()       { return display_unit_index_; }
    int                DisplayUnitIndex() const { return display_unit_index_; }

    // Transient GUI edit buffer (NOT serialized). The unit widget shows the
    // value in the user-selected display unit here and converts to/from the
    // canonical base-unit expr_str_ on commit. Kept as a member — tied to this
    // object's lifetime — rather than a pointer-keyed GUI map, so it can never
    // dangle when the owning container reallocates.
    std::string&       UiBuffer()            { return ui_buffer_; }

    void SetExpr(const std::string& s) { expr_str_ = s; last_error_.clear(); has_inline_unit_ = false; }
    void ClearExpr()                   { expr_str_.clear(); last_error_.clear(); has_inline_unit_ = false; }

    // Evaluate the expression and store the result in value_.
    // Returns true on success. No-op and returns true when no expression is set.
    // In audit mode (ctx.audit set) nothing is evaluated; see EvalContext::audit.
    bool Evaluate(const EvalContext& ctx);

    // Per-cell evaluation inside the solver residual (the constant / Arrhenius
    // mass-transfer coefficient). Evaluate() re-parses expr_str_ on every call
    // (strtod for a number, a fresh mu::Parser for a formula), so the residual
    // instead compiles once per run: CompileForRun() evaluates the expression;
    // one that does not use T, P, c, y or v is then a fixed value for the run,
    // and one that does keeps a single parser bound to the user variables'
    // values at compile time. EvaluateAt() returns the value at ctx's state and
    // falls back to Evaluate() when expr_str_ no longer matches the compiled
    // text (or nothing was compiled).
    void   CompileForRun(const UserVariableRegistry* registry);
    double EvaluateAt(const EvalContext& ctx);

    // True when the expression references a run-time state variable (T, P, c,
    // y, v). Such expressions are evaluated per cell during the run, so their
    // value outside the run (placeholder state) is not meaningful.
    bool UsesStateVariables() const;
    // True when CompileForRun resolved this parameter to a state-independent
    // value: EvaluateAt returns the same number for every EvalContext.
    bool IsRunConstant() const;

    // True when the cached value is inf/NaN and the expression does not depend
    // on run-time state, i.e. the run would use a non-finite constant.
    bool IsNonFiniteConstant() const;

    // Format a double as a shortest-round-trip decimal string.
    static std::string FormatScalar(double v);

private:
    double       value_     = 0.0;
    std::string  expr_str_;
    int          display_unit_index_ = 0;  // GUI display unit (0 == SI base)
    std::string  ui_buffer_;              // transient GUI edit buffer — not serialized
    std::string  last_error_;             // transient — not serialized
    bool         has_inline_unit_ = false; // transient — set by Evaluate
    UnitCategory inline_unit_category_ = UnitCategory::Dimensionless;

    // Transient compiled form (see CompileForRun). The parser binds addresses
    // inside its own RunCache, so a copy starts uncompiled.
    struct RunCache;
    struct RunCacheHolder {
        std::unique_ptr<RunCache> ptr;
        RunCacheHolder();
        ~RunCacheHolder();
        RunCacheHolder(const RunCacheHolder&);
        RunCacheHolder& operator=(const RunCacheHolder&);
        RunCacheHolder(RunCacheHolder&&) noexcept;
        RunCacheHolder& operator=(RunCacheHolder&&) noexcept;
    };
    RunCacheHolder run_cache_;

    friend class cereal::access;
    template<class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
        ar(cereal::make_nvp("Value",      value_),
           cereal::make_nvp("Expression", expr_str_));
        if (version >= 2) {
            ar(cereal::make_nvp("DisplayUnitIndex", display_unit_index_));
        }
    }
};

} // namespace adsim

CEREAL_CLASS_VERSION(adsim::ExpressionParam, 2);
