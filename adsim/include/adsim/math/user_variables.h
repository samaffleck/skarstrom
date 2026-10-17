/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <cstdint>

#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/access.hpp>
#include <cereal/cereal.hpp>

namespace adsim {

enum class UnitCategory {
    // ── Existing (do NOT reorder — serialised as int) ──────────────────────
    Dimensionless       = 0,
    Time                = 1,
    Length              = 2,
    Temperature         = 3,
    Pressure            = 4,
    Energy              = 5,
    SpecificMolarEnergy = 6,
    InverseTime         = 7,

    // ── Extended ─────────────────────────────────────────────────────────
    Velocity                  = 9,
    VolumetricFlowRate        = 10,
    MolarFlowRate             = 11,
    MassFlowRate              = 12,
    MoleFraction              = 13,
    Viscosity                 = 14,
    Density                   = 15,
    MolarConcentration        = 16,
    HeatCapacity              = 17,
    MolarMass                 = 18,
    ThermalConductivity       = 19,
    HeatTransferCoefficient   = 20,
    PressureGradient          = 21,
    AdsorbedConcentration     = 22,
    Mass                      = 23,
    Volume                    = 24,
    MolarAmount               = 25,
    Area                      = 26,
    Diffusivity               = 27,
    SpecificEnergy            = 28,

    Custom              = 8,   // keep last for combo display
};

// Human-readable name for a UnitCategory — used in error messages (front ends
// may have their own richer labels).
const char* UnitCategoryName(UnitCategory cat);

struct UserVariable {
    std::string  name;
    double       value         = 0.0;
    std::string  value_expr;           // expression string; empty when value is a plain number
    std::string  description;
    UnitCategory unit_category = UnitCategory::Dimensionless;

private:
    friend class cereal::access;
    template<class Archive>
    void save(Archive& ar, std::uint32_t const /*version*/) const {
        ar(cereal::make_nvp("Name",         name),
           cereal::make_nvp("Value",        value),
           cereal::make_nvp("ValueExpr",    value_expr),
           cereal::make_nvp("Description",  description),
           cereal::make_nvp("UnitCategory", static_cast<int>(unit_category)));
    }
    template<class Archive>
    void load(Archive& ar, std::uint32_t const version) {
        if (version >= 3) {
            int cat_int = 0;
            ar(cereal::make_nvp("Name",         name),
               cereal::make_nvp("Value",        value),
               cereal::make_nvp("ValueExpr",    value_expr),
               cereal::make_nvp("Description",  description),
               cereal::make_nvp("UnitCategory", cat_int));
            unit_category = static_cast<UnitCategory>(cat_int);
        } else if (version >= 2) {
            ar(cereal::make_nvp("Name",        name),
               cereal::make_nvp("Value",       value),
               cereal::make_nvp("ValueExpr",   value_expr),
               cereal::make_nvp("Description", description));
            unit_category = UnitCategory::Dimensionless;
        } else {
            int cat_int = 0;
            ar(cereal::make_nvp("Name",         name),
               cereal::make_nvp("Value",        value),
               cereal::make_nvp("UnitCategory", cat_int),
               cereal::make_nvp("Description",  description));
            value_expr.clear();
            unit_category = UnitCategory::Dimensionless;
        }
    }
};

class UserVariableRegistry {
public:
    void Add(UserVariable var);
    void Remove(const std::string& name);
    bool Contains(const std::string& name) const;
    UserVariable* Find(const std::string& name);
    const UserVariable* Find(const std::string& name) const;
    const std::vector<UserVariable>& GetAll() const { return vars_; }
    std::vector<UserVariable>&       GetAll()       { return vars_; }

    // Build a flat name→value map for expression evaluation
    std::unordered_map<std::string, double> BuildVarMap() const;

    // Returns true if name clashes with a built-in (T, P, c, y, v, R, PI, etc.)
    static bool IsReservedName(const std::string& name);

private:
    std::vector<UserVariable> vars_;

    friend class cereal::access;
    template<class Archive>
    void serialize(Archive& ar, std::uint32_t const /*version*/) {
        ar(cereal::make_nvp("Variables", vars_));
    }
};

} // namespace adsim

CEREAL_CLASS_VERSION(adsim::UserVariable,         3);
CEREAL_CLASS_VERSION(adsim::UserVariableRegistry, 1);
