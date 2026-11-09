// test_cereal_class_versions.cpp
//
// Pins the cereal class version of every serialised engine class.
//
// CEREAL_CLASS_VERSION specialises cereal::detail::Version<T>. A translation
// unit that does not see the specialisation uses the primary template and
// gets version 0. The serialize templates are inline in the public headers,
// so the version must be declared in the same header: otherwise which
// version a save writes (and which branch a load takes) depends on which
// instantiation the linker keeps.
//
//   * HeadersDeclareVersions: this file includes only public headers, so a
//     version declared in a .cpp reads 0 here and the check fails.
//   * SavedDecksWriteDeclaredVersions: saving each bundled deck writes the
//     same versions (polymorphic classes, their bases, and the plain classes
//     at fixed places in the deck). An unknown polymorphic class also fails,
//     so a new class must be added to the table.

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <cereal/archives/json.hpp>
#include <cereal/external/rapidjson/document.h>
#include <cereal/external/rapidjson/pointer.h>
#include <gtest/gtest.h>

#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/boundary_conditions/boundary_value.h"
#include "adsim/boundary_conditions/buffer_tank_source.h"
#include "adsim/boundary_conditions/connection.h"
#include "adsim/boundary_conditions/cycle_endpoint_summary.h"
#include "adsim/boundary_conditions/flow_boundary_value.h"
#include "adsim/boundary_conditions/step.h"
#include "adsim/equations/dispersion_model.h"
#include "adsim/equations/heat_transfer_model.h"
#include "adsim/equations/isotherm.h"
#include "adsim/equations/mass_transfer_model.h"
#include "adsim/equations/reaction_kinds.h"
#include "adsim/equations/wall_temperature_boundary_condition.h"
#include "adsim/math/expression_param.h"
#include "adsim/math/user_variables.h"
#include "adsim/objects/compressor.h"
#include "adsim/objects/fluid.h"
#include "adsim/objects/porous_media.h"
#include "adsim/objects/reaction.h"
#include "adsim/objects/reactor.h"
#include "adsim/objects/reactor_equations.h"
#include "adsim/objects/user_data.h"
#include "adsim/objects/valve.h"
#include "adsim/objects/void_region.h"
#include "adsim/objects/wall.h"
#include "adsim/solver/cycle.h"
#include "adsim/solver/cyclic_steady_state.h"
#include "adsim/solver/experimental_data.h"
#include "adsim/solver/ida_settings.h"
#include "adsim/solver/initial_conditions.h"
#include "adsim/solver/kpi_config.h"
#include "adsim/solver/sensitivity_study.h"
#include "adsim/solver/warning_settings.h"

namespace adsim {
namespace {

struct Pinned {
  std::string name;      // as cereal writes it in "polymorphic_name"
  std::uint32_t seen;    // what this translation unit sees
  std::uint32_t pinned;  // the version the class writes
};

#define ADSIM_PIN(T, V) Pinned{#T, cereal::detail::Version<T>::version, V}

// Bump the number here together with the class's CEREAL_CLASS_VERSION.
const std::vector<Pinned>& PinnedVersions() {
  static const std::vector<Pinned> table = {
      // objects
      ADSIM_PIN(adsim::UserData, 4),
      ADSIM_PIN(adsim::Fluid, 4),
      ADSIM_PIN(adsim::Bed, 4),
      ADSIM_PIN(adsim::BedWall, 4),
      ADSIM_PIN(adsim::VoidRegion, 2),
      ADSIM_PIN(adsim::Valve, 1),
      ADSIM_PIN(adsim::Adsorbent, 5),
      ADSIM_PIN(adsim::AdsorptionSource, 1),
      ADSIM_PIN(adsim::AdsorbentDirectory, 1),
      ADSIM_PIN(adsim::MixedAdsorbentLayer, 7),
      ADSIM_PIN(adsim::StructuredGeometry, 2),
      ADSIM_PIN(adsim::AdsorbentLayers, 1),
      ADSIM_PIN(adsim::EnergySource, 2),
      ADSIM_PIN(adsim::ReactorEquations, 2),
      ADSIM_PIN(adsim::LayerEquations, 1),
      ADSIM_PIN(adsim::EffectivenessModel, 1),
      ADSIM_PIN(adsim::StoichiometryEntry, 1),
      ADSIM_PIN(adsim::Reaction, 1),
      ADSIM_PIN(adsim::ReactionSource, 1),
      ADSIM_PIN(adsim::ConstantPumpPowerModel, 2),
      ADSIM_PIN(adsim::PolyPumpPowerModel, 2),
      ADSIM_PIN(adsim::IsothermalCompressorEnergyModel, 1),
      ADSIM_PIN(adsim::AdiabaticCompressorEnergyModel, 1),
      ADSIM_PIN(adsim::ConstantPumpEfficiencyModel, 2),
      ADSIM_PIN(adsim::VacuumPumpEfficiencyModel, 2),
      ADSIM_PIN(adsim::PolyPumpEfficiencyModel, 3),
      ADSIM_PIN(adsim::IPump, 3),
      ADSIM_PIN(adsim::Compressor, 2),
      ADSIM_PIN(adsim::VacuumPump, 2),
      // math
      ADSIM_PIN(adsim::UserVariable, 3),
      ADSIM_PIN(adsim::UserVariableRegistry, 1),
      ADSIM_PIN(adsim::ExpressionParam, 2),
      // equations
      ADSIM_PIN(adsim::DispersionModelConstant, 1),
      ADSIM_PIN(adsim::HeatDispersionModelConstant, 1),
      ADSIM_PIN(adsim::DispersionModelRuthven, 1),
      ADSIM_PIN(adsim::DispersionModelTaylorAris, 1),
      ADSIM_PIN(adsim::HeatDispersionModelRuthven, 1),
      ADSIM_PIN(adsim::HeatDispersionModelYagiandKunii, 1),
      ADSIM_PIN(adsim::HeatDispersionModelChannel, 1),
      ADSIM_PIN(adsim::ConstantParticleHeatTransferModel, 1),
      ADSIM_PIN(adsim::ConstantInternalHeatTransferModel, 1),
      ADSIM_PIN(adsim::InternalHeatTransferModel_Dixon, 1),
      ADSIM_PIN(adsim::InternalHeatTransferModel_Leva, 1),
      ADSIM_PIN(adsim::ConstantExternalHeatTransferModel, 1),
      ADSIM_PIN(adsim::ExternalHeatTransferModel_Holman, 1),
      ADSIM_PIN(adsim::FilmParticleHeatTransferModel, 1),
      ADSIM_PIN(adsim::InertIsotherm, 2),
      ADSIM_PIN(adsim::HenryIsotherm, 2),
      ADSIM_PIN(adsim::LangmuirIsotherm, 2),
      ADSIM_PIN(adsim::AntiLangmuirIsotherm, 2),
      ADSIM_PIN(adsim::FreundlichIsotherm, 2),
      ADSIM_PIN(adsim::SipsIsotherm, 2),
      ADSIM_PIN(adsim::Sips2Isotherm, 2),
      ADSIM_PIN(adsim::ExtendedSips2Isotherm, 2),
      ADSIM_PIN(adsim::DualSiteSipsIsotherm, 2),
      ADSIM_PIN(adsim::RedlichPetersonIsotherm, 2),
      ADSIM_PIN(adsim::TothIsotherm, 2),
      ADSIM_PIN(adsim::Toth2Isotherm, 2),
      ADSIM_PIN(adsim::ArrheniusSipsIsotherm, 2),
      ADSIM_PIN(adsim::ExtendedLangmuirFreundlichIsotherm, 1),
      ADSIM_PIN(adsim::GabIsotherm, 1),
      ADSIM_PIN(adsim::GabYoungIsotherm, 1),
      ADSIM_PIN(adsim::TothWadstIsotherm, 1),
      ADSIM_PIN(adsim::TothHumidIsotherm, 1),
      ADSIM_PIN(adsim::ExtendedLangmuirIsotherm, 2),
      ADSIM_PIN(adsim::DualSiteLangmuirIsotherm, 2),
      ADSIM_PIN(adsim::ExtendedDualSiteLangmuirIsotherm, 2),
      ADSIM_PIN(adsim::LangmuirIsotherm_Concentration, 2),
      ADSIM_PIN(adsim::ExtendedLangmuirIsotherm_Concentration, 2),
      ADSIM_PIN(adsim::DualSiteLangmuirIsotherm_Concentration, 2),
      ADSIM_PIN(adsim::ExtendedDualSiteLangmuirIsotherm_Concentration, 2),
      ADSIM_PIN(adsim::MassTransferCoefficientConstant, 3),
      ADSIM_PIN(adsim::MassTransferCoefficientMacroporeKnudsen, 3),
      ADSIM_PIN(adsim::MassTransferCoefficientMacropore, 3),
      ADSIM_PIN(adsim::MassTransferCoefficient_Haynes, 2),
      ADSIM_PIN(adsim::MassTransferCoefficientMacroporeCapacity, 2),
      ADSIM_PIN(adsim::MassTransferCoefficientSeriesResistance, 2),
      ADSIM_PIN(adsim::MassTransferResistances, 1),
      ADSIM_PIN(adsim::MassTransferCoefficientArrhenius, 3),
      ADSIM_PIN(adsim::LinearDrivingForce, 1),
      ADSIM_PIN(adsim::LinearDrivingForceSecondOrder, 1),
      ADSIM_PIN(adsim::DualRateLinearDrivingForce, 1),
      ADSIM_PIN(adsim::IReactionRate, 1),
      ADSIM_PIN(adsim::PowerLawRate, 1),
      ADSIM_PIN(adsim::ReversiblePowerLawRate, 1),
      ADSIM_PIN(adsim::LhhwRate, 1),
      ADSIM_PIN(adsim::ExpressionRate, 1),
      ADSIM_PIN(adsim::LhhwAdsorptionTerm, 1),
      ADSIM_PIN(adsim::WallTemperatureBoundaryConditionFixed, 1),
      ADSIM_PIN(adsim::WallTemperatureBoundaryConditionZeroFlux, 1),
      // boundary conditions
      ADSIM_PIN(adsim::ConstantBoundaryValue, 2),
      ADSIM_PIN(adsim::LinearBoundaryValue, 4),
      ADSIM_PIN(adsim::ExponentialBoundaryValue, 4),
      ADSIM_PIN(adsim::LogLinearBoundaryValue, 2),
      ADSIM_PIN(adsim::DoubleExponentialBoundaryValue, 1),
      ADSIM_PIN(adsim::VelocityBoundaryValue, 2),
      ADSIM_PIN(adsim::VolumetricFlowBoundaryValue, 2),
      ADSIM_PIN(adsim::MolarFlowBoundaryValue, 2),
      ADSIM_PIN(adsim::MassFlowBoundaryValue, 2),
      ADSIM_PIN(adsim::IBoundaryStream, 5),
      ADSIM_PIN(adsim::BoundaryStream, 1),
      ADSIM_PIN(adsim::ValveBoundaryStream, 2),
      ADSIM_PIN(adsim::PumpBoundaryStream, 1),
      ADSIM_PIN(adsim::PressureBoundaryCondition, 5),
      ADSIM_PIN(adsim::FlowBoundaryCondition, 1),
      ADSIM_PIN(adsim::ClosedBoundaryCondition, 1),
      ADSIM_PIN(adsim::HeatingMedium, 1),
      ADSIM_PIN(adsim::Step, 6),
      ADSIM_PIN(adsim::BufferTankSource, 3),
      ADSIM_PIN(adsim::CycleEndpointSummary, 1),
      ADSIM_PIN(adsim::Connection, 1),
      // solver
      ADSIM_PIN(adsim::Cycle, 2),
      ADSIM_PIN(adsim::ExperimentalData, 1),
      ADSIM_PIN(adsim::CyclicSteadyState, 2),
      ADSIM_PIN(adsim::AccelerateCyclicSteadyState, 0),
      ADSIM_PIN(adsim::IDASettings, 7),
      ADSIM_PIN(adsim::InitialConditions, 3),
      ADSIM_PIN(adsim::WarningSettings, 5),
      ADSIM_PIN(adsim::KPIEndpointRef, 1),
      ADSIM_PIN(adsim::ComponentKPIConfig, 2),
      ADSIM_PIN(adsim::KPIConfig, 4),
      ADSIM_PIN(adsim::ReactionKPIConfig, 0),
      ADSIM_PIN(adsim::sensitivity::Input, 0),
      ADSIM_PIN(adsim::sensitivity::RunPolicy, 0),
      ADSIM_PIN(adsim::sensitivity::Study, 0),
  };
  return table;
}

#undef ADSIM_PIN

std::uint32_t PinnedVersion(const std::string& name) {
  for (const auto& p : PinnedVersions())
    if (p.name == name) return p.pinned;
  ADD_FAILURE() << name << " is not in the pinned version table";
  return 0;
}

TEST(CerealClassVersions, HeadersDeclareVersions) {
  for (const auto& p : PinnedVersions()) {
    EXPECT_EQ(p.seen, p.pinned)
        << p.name << ": the public headers give version " << p.seen << ", expected " << p.pinned
        << ". Put its CEREAL_CLASS_VERSION in the header that defines its serialize template.";
  }
}

// Base class of each polymorphic class, as written in the derived class's
// "value0" (cereal::base_class); interfaces without a declared version write 0.
std::uint32_t PinnedBaseVersion(const std::string& derived) {
  static const std::map<std::string, std::string> bases = {
      {"adsim::BoundaryStream", "adsim::IBoundaryStream"},
      {"adsim::ValveBoundaryStream", "adsim::IBoundaryStream"},
      {"adsim::PumpBoundaryStream", "adsim::IBoundaryStream"},
      {"adsim::Compressor", "adsim::IPump"},
      {"adsim::VacuumPump", "adsim::IPump"},
      {"adsim::PowerLawRate", "adsim::IReactionRate"},
      {"adsim::ReversiblePowerLawRate", "adsim::IReactionRate"},
      {"adsim::LhhwRate", "adsim::IReactionRate"},
      {"adsim::ExpressionRate", "adsim::IReactionRate"},
  };
  auto it = bases.find(derived);
  return it == bases.end() ? 0 : PinnedVersion(it->second);
}

// cereal writes a class's version only where the class first appears in an
// archive, so absent fields are skipped.
void CheckVersion(const rapidjson::Value& v, std::uint32_t expected, const std::string& where) {
  if (!v.IsObject() || !v.HasMember("cereal_class_version")) return;
  EXPECT_EQ(v["cereal_class_version"].GetUint(), expected) << where;
}

void WalkPolymorphic(const rapidjson::Value& v, const std::string& deck, int& seen) {
  if (v.IsObject()) {
    if (v.HasMember("polymorphic_name") && v.HasMember("ptr_wrapper") && v["ptr_wrapper"].HasMember("data")) {
      const std::string name = v["polymorphic_name"].GetString();
      const auto& data = v["ptr_wrapper"]["data"];
      if (data.IsObject() && data.HasMember("cereal_class_version")) ++seen;
      CheckVersion(data, PinnedVersion(name), deck + ": " + name);
      if (data.IsObject() && data.HasMember("value0"))
        CheckVersion(data["value0"], PinnedBaseVersion(name), deck + ": base of " + name);
    }
    for (auto it = v.MemberBegin(); it != v.MemberEnd(); ++it) WalkPolymorphic(it->value, deck, seen);
  } else if (v.IsArray()) {
    for (const auto& e : v.GetArray()) WalkPolymorphic(e, deck, seen);
  }
}

// Plain members whose place in a deck varies (inside polymorphic boundary
// conditions and streams), found by member name.
void WalkNamed(const rapidjson::Value& v, const std::string& deck, int& seen) {
  static const std::map<std::string, std::string> named = {
      {"Boundary_Condition_Connection", "adsim::Connection"},
      {"Boundary_Connection", "adsim::Connection"},
  };
  if (v.IsObject()) {
    for (auto it = v.MemberBegin(); it != v.MemberEnd(); ++it) {
      auto cls = named.find(it->name.GetString());
      if (cls != named.end()) {
        if (it->value.IsObject() && it->value.HasMember("cereal_class_version")) ++seen;
        CheckVersion(it->value, PinnedVersion(cls->second), deck + ": " + cls->first + " (" + cls->second + ")");
      }
      WalkNamed(it->value, deck, seen);
    }
  } else if (v.IsArray()) {
    for (const auto& e : v.GetArray()) WalkNamed(e, deck, seen);
  }
}

std::vector<std::filesystem::path> ExampleDecks() {
  const auto dir = std::filesystem::current_path() / ".." / "examples" / "decks";
  std::vector<std::filesystem::path> decks;
  for (const auto& e : std::filesystem::directory_iterator(dir))
    if (e.path().extension() == ".json") decks.push_back(e.path());
  std::sort(decks.begin(), decks.end());
  return decks;
}

TEST(CerealClassVersions, SavedDecksWriteDeclaredVersions) {
  // Plain (non-polymorphic) members at fixed places in a deck.
  const std::vector<std::pair<const char*, const char*>> fixed = {
      {"/value0", "adsim::UserData"},
      {"/value0/Fluid", "adsim::Fluid"},
      {"/value0/Adsorbent_Directory", "adsim::AdsorbentDirectory"},
      {"/value0/Adsorbent_Directory/Adsorbents/0/value", "adsim::Adsorbent"},
      {"/value0/Adsorbent_Directory/Adsorbents/0/value/Adsorption_Source", "adsim::AdsorptionSource"},
      {"/value0/Bed", "adsim::Bed"},
      {"/value0/Bed/Wall", "adsim::BedWall"},
      {"/value0/Bed/Layers", "adsim::AdsorbentLayers"},
      {"/value0/Bed/Layers/Layers/0", "adsim::MixedAdsorbentLayer"},
      {"/value0/Bed/Layers/Layers/0/Structured_Geometry", "adsim::StructuredGeometry"},
      {"/value0/Bed/Layers/Layers/0/Layer_Equations", "adsim::LayerEquations"},
      {"/value0/Bed/Equations", "adsim::ReactorEquations"},
      {"/value0/Bed/Bottom_Void", "adsim::VoidRegion"},
      {"/value0/Initial_Conditions", "adsim::InitialConditions"},
      {"/value0/Cycle", "adsim::Cycle"},
      {"/value0/Cycle/Step_List/0", "adsim::Step"},
      {"/value0/Cycle/Step_List/0/Heating_Medium", "adsim::HeatingMedium"},
      {"/value0/Solver_Settings", "adsim::IDASettings"},
      {"/value0/Solver_Settings/Cyclic_Steady_State", "adsim::CyclicSteadyState"},
      {"/value0/Solver_Settings/Accelerate_Cyclic_Steady_State", "adsim::AccelerateCyclicSteadyState"},
      {"/value0/Experimental_Data", "adsim::ExperimentalData"},
      {"/value0/KPI_Config", "adsim::KPIConfig"},
      {"/value0/Variable_Registry", "adsim::UserVariableRegistry"},
  };

  const auto decks = ExampleDecks();
  ASSERT_FALSE(decks.empty());
  int polymorphic_seen = 0;
  int named_seen = 0;
  for (const auto& deck : decks) {
    SCOPED_TRACE(deck.filename().string());
    UserData data;
    ASSERT_NO_THROW(LoadFromJson(data, deck.string()));
    std::ostringstream oss;
    {
      cereal::JSONOutputArchive ar(oss);
      ar(data);
    }
    rapidjson::Document doc;
    doc.Parse<rapidjson::kParseFullPrecisionFlag | rapidjson::kParseNanAndInfFlag>(oss.str().c_str());
    ASSERT_FALSE(doc.HasParseError());

    for (const auto& [pointer, cls] : fixed) {
      const rapidjson::Value* v = rapidjson::Pointer(pointer).Get(doc);
      ASSERT_NE(v, nullptr) << pointer;
      ASSERT_TRUE(v->IsObject() && v->HasMember("cereal_class_version")) << pointer << " has no version";
      CheckVersion(*v, PinnedVersion(cls), std::string(pointer) + " (" + cls + ")");
    }
    // The fields gated on the versions of issue-prone classes are written.
    const auto* source = rapidjson::Pointer("/value0/Adsorbent_Directory/Adsorbents/0/value/Adsorption_Source").Get(doc);
    EXPECT_TRUE(source->HasMember("Mixing_Rule"));
    const auto* layer = rapidjson::Pointer("/value0/Bed/Layers/Layers/0").Get(doc);
    EXPECT_TRUE(layer->HasMember("Momentum_Model"));

    WalkPolymorphic(doc, deck.filename().string(), polymorphic_seen);
    WalkNamed(doc, deck.filename().string(), named_seen);
  }
  EXPECT_GT(polymorphic_seen, 0);
  EXPECT_GT(named_seen, 0);
}

}  // namespace
}  // namespace adsim
