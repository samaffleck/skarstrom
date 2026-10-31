#include "adsim/math/expression_param.h"
#include "adsim/math/unit_parser.h"

#include <muParser.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>

namespace adsim {

bool ExpressionParam::HasExpression() const {
    if (expr_str_.empty()) return false;
    const char* s = expr_str_.c_str();
    // Skip leading whitespace
    while (*s == ' ' || *s == '\t' || *s == '\n') ++s;
    if (*s == '\0') return false;
    char* end = nullptr;
    std::strtod(s, &end);
    if (end == s) return true;  // didn't consume anything → not a number
    // Skip trailing whitespace
    while (*end == ' ' || *end == '\t' || *end == '\n') ++end;
    return *end != '\0';  // text after the number → it's an expression
}

std::string ExpressionParam::FormatScalar(double v) {
    // Find the shortest %g precision that round-trips exactly back to v.
    char buf[32];
    for (int p = 1; p <= 17; ++p) {
        std::snprintf(buf, sizeof(buf), "%.*g", p, v);
        if (std::strtod(buf, nullptr) == v) return std::string(buf);
    }
    return std::string(buf);
}

void NonFiniteExpressionAudit::Record(const ExpressionParam& p) {
    for (const auto& e : entries)
        if (e.param == &p) return;
    entries.push_back({&p, scope, p.GetExpr(), p.Get()});
}

namespace {
// Names Evaluate() binds to the run-time state (EvalContext T, P, c, y, v).
constexpr const char* kStateVariables[] = {"T", "P", "c", "y", "v"};
} // namespace

bool ExpressionParam::UsesStateVariables() const {
    if (!HasExpression()) return false;
    try {
        mu::Parser p;
        p.DefineConst("R",  8.31446261815324);
        p.DefineConst("PI", 3.14159265358979);
        std::string math_expr = expr_str_;
        UnitXform xform;
        ParseUnitSuffix(expr_str_, math_expr, xform, nullptr);
        p.SetExpr(math_expr);
        // GetUsedVar parses without evaluating and reports undefined names too,
        // so user variables need not be defined here.
        const auto& used = p.GetUsedVar();
        for (const char* name : kStateVariables)
            if (used.count(name)) return true;
    } catch (const mu::Parser::exception_type&) {
        // Unparseable: Evaluate() reports the error; not state-dependent.
    }
    return false;
}

bool ExpressionParam::IsNonFiniteConstant() const {
    return !std::isfinite(value_) && !UsesStateVariables();
}

bool ExpressionParam::Evaluate(const EvalContext& ctx) {
    if (ctx.audit) {
        if (IsNonFiniteConstant()) ctx.audit->Record(*this);
        return true;
    }

    // Fast path: when expr_str_ is empty or holds a plain numeric literal,
    // skip the entire muParser setup. value_ is already correct after
    // operator=(double) / ctor, but expr_str_ can be set externally (SetExpr,
    // GUI binding via GetExprMutable, or load() from JSON), so sync value_
    // from the string. One strtod call is orders of magnitude cheaper than
    // constructing a mu::Parser + parse + eval.
    if (!HasExpression()) {
        if (!expr_str_.empty()) {
            char* end = nullptr;
            const double parsed = std::strtod(expr_str_.c_str(), &end);
            if (end != expr_str_.c_str()) value_ = parsed;
        }
        last_error_.clear();
        has_inline_unit_ = false;
        return true;
    }

    try {
        mu::Parser p;

        // Built-in physical constants
        p.DefineConst("R",  8.31446261815324);
        p.DefineConst("PI", 3.14159265358979);

        // Simulation state variables — local copies so p can take their address
        double T = ctx.T;
        double P = ctx.P;
        double c = ctx.c;
        double y = ctx.y;
        double v = ctx.v;
        p.DefineVar("T", &T);
        p.DefineVar("P", &P);
        p.DefineVar("c", &c);
        p.DefineVar("y", &y);
        p.DefineVar("v", &v);

        // User-defined variables
        std::unordered_map<std::string, double> user_vars;
        if (ctx.registry) {
            user_vars = ctx.registry->BuildVarMap();
        }
        for (auto& [name, val] : user_vars) {
            p.DefineVar(name.c_str(), &val);
        }

        // Strip optional trailing [unit] suffix (COMSOL-style: "10 [bar]")
        std::string math_expr = expr_str_;
        UnitXform xform;
        UnitCategory parsed_cat = UnitCategory::Dimensionless;
        bool has_unit = ParseUnitSuffix(expr_str_, math_expr, xform, &parsed_cat);

        p.SetExpr(math_expr);
        double raw = p.Eval();
        value_ = has_unit ? xform.ToBase(raw) : raw;
        has_inline_unit_ = has_unit;
        if (has_unit) inline_unit_category_ = parsed_cat;
        last_error_.clear();
        return true;
    } catch (const mu::Parser::exception_type& e) {
        last_error_ = e.GetMsg();
        return false;
    }
}

struct ExpressionParam::RunCache {
    std::string expr;                    // expr_str_ this was compiled from
    std::unique_ptr<mu::Parser> parser;  // null: state-independent, value_ is fixed
    // Storage the parser's variables are bound to.
    double T = 0.0, P = 0.0, c = 0.0, y = 0.0, v = 0.0;
    std::vector<std::pair<std::string, double>> user_vars;
    UnitXform xform;
    bool has_unit = false;
};

ExpressionParam::RunCacheHolder::RunCacheHolder() = default;
ExpressionParam::RunCacheHolder::~RunCacheHolder() = default;
ExpressionParam::RunCacheHolder::RunCacheHolder(const RunCacheHolder&) {}
ExpressionParam::RunCacheHolder& ExpressionParam::RunCacheHolder::operator=(const RunCacheHolder&) {
    ptr.reset();
    return *this;
}
ExpressionParam::RunCacheHolder::RunCacheHolder(RunCacheHolder&&) noexcept = default;
ExpressionParam::RunCacheHolder& ExpressionParam::RunCacheHolder::operator=(RunCacheHolder&&) noexcept = default;

void ExpressionParam::CompileForRun(const UserVariableRegistry* registry) {
    run_cache_.ptr.reset();
    EvalContext ctx;
    ctx.registry = registry;
    // A failed evaluation stays uncompiled, so EvaluateAt keeps going through
    // Evaluate() and reports the error the same way.
    if (!Evaluate(ctx)) return;

    auto cache = std::make_unique<RunCache>();
    cache->expr = expr_str_;
    if (UsesStateVariables()) {
        try {
            auto p = std::make_unique<mu::Parser>();
            p->DefineConst("R",  8.31446261815324);
            p->DefineConst("PI", 3.14159265358979);
            p->DefineVar("T", &cache->T);
            p->DefineVar("P", &cache->P);
            p->DefineVar("c", &cache->c);
            p->DefineVar("y", &cache->y);
            p->DefineVar("v", &cache->v);
            if (registry) {
                // Same name -> value resolution as Evaluate (last duplicate wins).
                for (const auto& [name, val] : registry->BuildVarMap())
                    cache->user_vars.emplace_back(name, val);
            }
            for (auto& [name, val] : cache->user_vars)
                p->DefineVar(name.c_str(), &val);

            std::string math_expr = expr_str_;
            cache->has_unit = ParseUnitSuffix(expr_str_, math_expr, cache->xform, nullptr);
            p->SetExpr(math_expr);
            cache->parser = std::move(p);
        } catch (const mu::Parser::exception_type&) {
            return;
        }
    }
    run_cache_.ptr = std::move(cache);
}

bool ExpressionParam::IsRunConstant() const {
    const RunCache* rc = run_cache_.ptr.get();
    return rc && rc->expr == expr_str_ && !rc->parser;
}

double ExpressionParam::EvaluateAt(const EvalContext& ctx) {
    RunCache* rc = run_cache_.ptr.get();
    if (!rc || rc->expr != expr_str_) {
        Evaluate(ctx);
        return value_;
    }
    if (!rc->parser) return value_;

    rc->T = ctx.T;
    rc->P = ctx.P;
    rc->c = ctx.c;
    rc->y = ctx.y;
    rc->v = ctx.v;
    try {
        const double raw = rc->parser->Eval();
        value_ = rc->has_unit ? rc->xform.ToBase(raw) : raw;
    } catch (const mu::Parser::exception_type& e) {
        last_error_ = e.GetMsg();
    }
    return value_;
}

} // namespace adsim
