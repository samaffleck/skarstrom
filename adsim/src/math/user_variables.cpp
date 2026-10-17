#include "adsim/math/user_variables.h"

#include <algorithm>

namespace adsim {

const char* UnitCategoryName(UnitCategory cat) {
    switch (cat) {
        case UnitCategory::Dimensionless:           return "Dimensionless";
        case UnitCategory::Time:                    return "Time";
        case UnitCategory::Length:                  return "Length";
        case UnitCategory::Temperature:             return "Temperature";
        case UnitCategory::Pressure:                return "Pressure";
        case UnitCategory::Energy:                  return "Energy";
        case UnitCategory::SpecificMolarEnergy:     return "Specific Molar Energy";
        case UnitCategory::InverseTime:             return "Inverse Time";
        case UnitCategory::Velocity:                return "Velocity";
        case UnitCategory::VolumetricFlowRate:      return "Volumetric Flow Rate";
        case UnitCategory::MolarFlowRate:           return "Molar Flow Rate";
        case UnitCategory::MassFlowRate:            return "Mass Flow Rate";
        case UnitCategory::MoleFraction:            return "Mole Fraction";
        case UnitCategory::Viscosity:               return "Viscosity";
        case UnitCategory::Density:                 return "Density";
        case UnitCategory::MolarConcentration:      return "Molar Concentration";
        case UnitCategory::HeatCapacity:            return "Heat Capacity";
        case UnitCategory::MolarMass:               return "Molar Mass";
        case UnitCategory::ThermalConductivity:     return "Thermal Conductivity";
        case UnitCategory::HeatTransferCoefficient: return "Heat Transfer Coefficient";
        case UnitCategory::PressureGradient:        return "Pressure Gradient";
        case UnitCategory::AdsorbedConcentration:   return "Adsorbed Concentration";
        case UnitCategory::Mass:                    return "Mass";
        case UnitCategory::Volume:                  return "Volume";
        case UnitCategory::MolarAmount:             return "Molar Amount";
        case UnitCategory::Area:                    return "Area";
        case UnitCategory::Diffusivity:             return "Diffusivity";
        case UnitCategory::SpecificEnergy:          return "Specific Energy";
        case UnitCategory::Custom:                  return "Custom";
    }
    return "Unknown";
}

void UserVariableRegistry::Add(UserVariable var) {
    auto it = std::find_if(vars_.begin(), vars_.end(),
                           [&](const UserVariable& v) { return v.name == var.name; });
    if (it != vars_.end()) {
        *it = std::move(var);
    } else {
        vars_.push_back(std::move(var));
    }
}

void UserVariableRegistry::Remove(const std::string& name) {
    vars_.erase(
        std::remove_if(vars_.begin(), vars_.end(),
                       [&](const UserVariable& v) { return v.name == name; }),
        vars_.end());
}

bool UserVariableRegistry::Contains(const std::string& name) const {
    return Find(name) != nullptr;
}

UserVariable* UserVariableRegistry::Find(const std::string& name) {
    auto it = std::find_if(vars_.begin(), vars_.end(),
                           [&](const UserVariable& v) { return v.name == name; });
    return it != vars_.end() ? &(*it) : nullptr;
}

const UserVariable* UserVariableRegistry::Find(const std::string& name) const {
    auto it = std::find_if(vars_.begin(), vars_.end(),
                           [&](const UserVariable& v) { return v.name == name; });
    return it != vars_.end() ? &(*it) : nullptr;
}

std::unordered_map<std::string, double> UserVariableRegistry::BuildVarMap() const {
    std::unordered_map<std::string, double> map;
    map.reserve(vars_.size());
    for (const auto& v : vars_) {
        map[v.name] = v.value;
    }
    return map;
}

bool UserVariableRegistry::IsReservedName(const std::string& name) {
    static const std::vector<std::string> kReserved = {
        "T", "P", "c", "y", "v", "R", "PI",
        // muparser built-in functions that look like identifiers
        "sin", "cos", "tan", "asin", "acos", "atan", "atan2",
        "sinh", "cosh", "tanh",
        "exp", "log", "log2", "log10", "sqrt", "abs",
        "min", "max", "pow", "sign"
    };
    return std::find(kReserved.begin(), kReserved.end(), name) != kReserved.end();
}

} // namespace adsim
