#pragma once

// STL includes
#include <algorithm>
#include <cctype>
#include <filesystem>
#include <unordered_map>
#include <vector>
#include <memory>
#include <cmath>

// AdSim includes
#include "../src/pch.h"
#include "adsim/objects/user_data.h"
#include "adsim/solver/coupled_solver.h"
#include "adsim/solver/kpis.h"
#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/boundary_conditions/boundary_value.h"
#include "adsim/boundary_conditions/flow_boundary_value.h"
#include "adsim/boundary_conditions/connection.h"
#include "adsim/equations/isotherm.h"
#include "adsim/math/constants.h"

// GTest includes
#include <gtest/gtest.h>

namespace adsim {

  // Build KPI endpoint refs for a set of (step, bed-end) faces by selecting every
  // endpoint (all streams of a flow BC, or the pressure boundary) on that face.
  // Mirrors the historical whole-face KPI selection; used by the paper tests that
  // build KPIParams directly.
  inline std::vector<KPIEndpointRef> KpiEndpointsForFaces(
      const Cycle& cycle,
      const std::vector<std::pair<int, BoundaryLocation>>& faces) {
    std::vector<KPIEndpointRef> out;
    for (const auto& [step, loc] : faces) {
      if (step < 0 || step >= static_cast<int>(cycle.step_list_.size())) continue;
      for (const auto& ep : cycle.step_list_[step].GetEndpoints())
        if (ep.location == loc) out.push_back({step, ep.name, ep.location});
    }
    return out;
  }

  inline bool is_strictly_increasing(const std::vector<double>& t) {
    if (t.size() < 2) 
      return true;
    for (size_t i = 1; i < t.size(); ++i) {
      if (!(t[i] > t[i - 1])) return false;
    }
    return true;
  }

  inline double interp_linear(const std::vector<double>& t,
                              const std::vector<double>& y,
                              double t_query)
  {
    auto it = std::upper_bound(t.begin(), t.end(), t_query);
    size_t idx1 = static_cast<size_t>(std::distance(t.begin(), it));
    if (idx1 == 0) return y.front();
    if (idx1 >= t.size()) return y.back();
    
    size_t i0 = idx1 - 1;
    size_t i1 = idx1;
    double t0 = t[i0], t1 = t[i1];
    double y0 = y[i0], y1 = y[i1];
    double w = (t_query - t0) / (t1 - t0);
    return y0 + w * (y1 - y0);
  }

  class UserDataTest : public ::testing::Test {
  public:
    UserData data;
    std::filesystem::path test_path = std::filesystem::current_path() / "Results";
    void SetUp() override {
      data.bed.wall_.thickness_ = 0.0175;
      data.bed.wall_.heat_capacity_ = 537;
      data.bed.wall_.density_ = 7800;
      data.bed.wall_.thermal_conductivity_ = 16.2;
    }

    void AddComponents(const std::vector<std::string>& components) {
      for (const std::string& component : components)
        data.fluid.AddComponent(component);
    }

    // Adds a layer to the bed and a corresponding adsorbent material to the directory.
    // Returns {LayerID, AdsorbentID} for the newly created layer and adsorbent.
    // The adsorbent is assigned 100% composition in the layer.
    std::pair<LayerID, AdsorbentID> AddLayer(const std::string& layer_name) {
      // Create adsorbent material in directory
      auto ads_id = data.adsorbent_directory.Add(layer_name);

      // Create layer in bed
      auto layer_id = data.bed.layers_.Add(layer_name);
      auto& layer = data.bed.layers_.Get(layer_id);
      layer.momentum_model_ = MomentumModel::KarmenCozeny;  // fixture predates Ergun: keep its Darcy / Kozeny-Carman pressure drop
      layer.domain_.length = 1.0;
      layer.particle_diameter_ = 0.0011;
      layer.eb_ = 0.37;
      layer.domain_.n_cells = 20;
      layer.composition_[ads_id] = 1.0;
      data.adsorbent_directory.Get(ads_id).ep_ = 0.37;

      return {layer_id, ads_id};
    }

    MixedAdsorbentLayer& GetLayer(LayerID layer_id) {
      return data.bed.layers_.Get(layer_id);
    }

    Adsorbent& GetAdsorbent(AdsorbentID ads_id) {
      return data.adsorbent_directory.Get(ads_id);
    }

    void SetHOA(AdsorbentID ads_id, const std::string& component, double HOA) {
      GetAdsorbent(ads_id).adsorption_source_.heat_of_adsorption_[component] = HOA;
    }

    void AddIsotherm(const std::vector<std::string>& components, AdsorbentID ads_id, std::unique_ptr<IIsotherm> isotherm) {
      GetAdsorbent(ads_id).adsorption_source_.AddIsotherm(components, std::move(isotherm));
    }

    void AddLinearDrivingForceModel(AdsorbentID ads_id, const std::string& component, double mtc) {
      auto ldf = std::make_unique<LinearDrivingForce>();
      auto ldf_mtc = std::make_unique<MassTransferCoefficientConstant>();
      ldf_mtc->constant_mtc_ = mtc;
      ldf->mtc_model_ = std::move(ldf_mtc);
      GetAdsorbent(ads_id).adsorption_source_.mass_transfer_models_[component] = std::move(ldf);
    }

    void AddLinearDrivingForceModelArrhenius(AdsorbentID ads_id, const std::string& component, double pre_exp_factor, double activation_energy) {
      auto ldf = std::make_unique<LinearDrivingForce>();
      auto ldf_mtc = std::make_unique<MassTransferCoefficientArrhenius>();
      ldf_mtc->pre_exp_factor_ = pre_exp_factor;
      ldf_mtc->activation_energy_ = activation_energy;
      ldf->mtc_model_ = std::move(ldf_mtc);
      GetAdsorbent(ads_id).adsorption_source_.mass_transfer_models_[component] = std::move(ldf);
    }

    void AddLinearDrivingForceModel_Haynes(AdsorbentID ads_id, const std::string& component,
      double tortuosity, double ep) {
      auto ldf = std::make_unique<LinearDrivingForce>();
      auto ldf_mtc = std::make_unique<MassTransferCoefficient_Haynes>();
      ldf_mtc->tortuosity_ = tortuosity;
      ldf_mtc->macropore_void_fraction_ = ep;
      ldf->mtc_model_ = std::move(ldf_mtc);
      GetAdsorbent(ads_id).adsorption_source_.mass_transfer_models_[component] = std::move(ldf);
    }

    void AddLinearDrivingForceModelSecondOrder(AdsorbentID ads_id, const std::string& component, double mtc) {
      auto ldf = std::make_unique<LinearDrivingForceSecondOrder>();
      auto ldf_mtc = std::make_unique<MassTransferCoefficientConstant>();
      ldf_mtc->constant_mtc_ = mtc;
      ldf->mtc_model_ = std::move(ldf_mtc);
      GetAdsorbent(ads_id).adsorption_source_.mass_transfer_models_[component] = std::move(ldf);
    }

    void AddLinearDrivingForceModel_MacroporeKnudsen(AdsorbentID ads_id, const std::string& component,
      double pore_radius, double tortuosity, double ep) {
      auto ldf = std::make_unique<LinearDrivingForce>();
      auto ldf_mtc = std::make_unique<MassTransferCoefficientMacroporeKnudsen>();
      ldf_mtc->pore_radius_ = pore_radius;
      ldf_mtc->tortuosity_ = tortuosity;
      ldf_mtc->macropore_void_fraction_ = ep;
      ldf->mtc_model_ = std::move(ldf_mtc);
      GetAdsorbent(ads_id).adsorption_source_.mass_transfer_models_[component] = std::move(ldf);
    }

    void AddLinearDrivingForceModel_Macropore(AdsorbentID ads_id, const std::string& component,
      double tortuosity, double ep) {
      auto ldf = std::make_unique<LinearDrivingForce>();
      auto ldf_mtc = std::make_unique<MassTransferCoefficientMacropore>();
      ldf_mtc->tortuosity_ = tortuosity;
      ldf_mtc->macropore_void_fraction_ = ep;
      ldf->mtc_model_ = std::move(ldf_mtc);
      GetAdsorbent(ads_id).adsorption_source_.mass_transfer_models_[component] = std::move(ldf);
    }

    void AddConstantInternalHeatTransferCoefficientModel(double htc, LayerID layer_id) {
      auto ihtc = std::make_unique<ConstantInternalHeatTransferModel>();
      ihtc->heat_transfer_coefficient_ = htc;
      GetLayer(layer_id).equations_.internal_heat_transfer_model_ = std::move(ihtc);
    }

    void AddConstantExternalHeatTransferCoefficientModel(double htc, LayerID layer_id) {
      auto ihtc = std::make_unique<ConstantExternalHeatTransferModel>();
      ihtc->heat_transfer_coefficient_ = htc;
      GetLayer(layer_id).equations_.external_heat_transfer_model_ = std::move(ihtc);
    }

    void AddInitialConditions(double P0, double T0, double Tw0, const std::unordered_map<std::string, double>& yi0, double u0=0.0) {
      data.ic.pressure = P0;
      data.ic.temperature = T0;
      data.ic.wall_temperature = Tw0;
      data.ic.velocity = u0;
      for (const auto& [component, yi] : yi0) {
        data.ic.mole_fraction[component] = yi;
      }
    }

    // Boundary condition helpers
    std::unique_ptr<IFlowBoundaryValue> CreateVelocityBC(double vel) {
      auto velocity = std::make_unique<VelocityBoundaryValue>();
      static_cast<ConstantBoundaryValue*>(velocity->profile_.get())->constant_value_ = vel;
      return std::move(velocity);
    }

    std::unique_ptr<IFlowBoundaryValue> CreateVolumeFlowBC(double volume_flow) {
      auto flow = std::make_unique<VolumetricFlowBoundaryValue>();
      static_cast<ConstantBoundaryValue*>(flow->profile_.get())->constant_value_ = volume_flow;
      return std::move(flow);
    }

    std::unique_ptr<IFlowBoundaryValue> CreateMolarFlowBC(double molar_flow) {
      auto flow = std::make_unique<MolarFlowBoundaryValue>();
      static_cast<ConstantBoundaryValue*>(flow->profile_.get())->constant_value_ = molar_flow;
      return std::move(flow);
    }

    std::unique_ptr<IFlowBoundaryValue> CreateMassFlowBC(double mass_flow) {
      auto flow = std::make_unique<MassFlowBoundaryValue>();
      static_cast<ConstantBoundaryValue*>(flow->profile_.get())->constant_value_ = mass_flow;
      return std::move(flow);
    }

    std::unique_ptr<BoundaryStream> CreateNoValveInletStream(std::unique_ptr<IFlowBoundaryValue> flowBC,
      const std::unordered_map<std::string, double>& yi,
      double temperature) {
      auto stream = std::make_unique<BoundaryStream>();
      stream->boundary_flow_ = std::move(flowBC);
      for (const auto& [component, y] : yi) {
        stream->mole_fraction[component] = y;
      }
      stream->temperature_ = temperature;
      return std::move(stream);
    }

    std::unique_ptr<ValveBoundaryStream> CreateValveInletStream(std::unique_ptr<IBoundaryValue> upstreamPressure,
      double CV,
      const std::unordered_map<std::string, double>& yi,
      double temperature) {
      auto stream = std::make_unique<ValveBoundaryStream>();
      stream->boundary_value_ = std::move(upstreamPressure);
      for (const auto& [component, y] : yi) {
        stream->mole_fraction[component] = y;
      }
      stream->temperature_ = temperature;
      stream->valve_.cv_ = CV;
      return std::move(stream);
    }

    std::unique_ptr<ValveBoundaryStream> CreateValveOutletStream(std::unique_ptr<IBoundaryValue> downstreamPressure,
      double CV) {
      auto stream = std::make_unique<ValveBoundaryStream>();
      stream->boundary_value_ = std::move(downstreamPressure);
      stream->valve_.cv_ = CV;
      return std::move(stream);
    }

    std::unique_ptr<BoundaryStream> CreateNoValveOutletStream(std::unique_ptr<IFlowBoundaryValue> flowBC) {
      auto stream = std::make_unique<BoundaryStream>();
      stream->boundary_flow_ = std::move(flowBC);
      return std::move(stream);
    }

    std::unique_ptr<ConstantBoundaryValue> CreateConstantBoundaryValue(double value) {
      auto bv = std::make_unique<ConstantBoundaryValue>();
      bv->constant_value_ = value;
      return std::move(bv);
    }

    std::unique_ptr<LinearBoundaryValue> CreateLinearBoundaryValue(double start_value, double end_value, double ramp_time) {
      auto bv = std::make_unique<LinearBoundaryValue>();
      bv->start_value_ = start_value;
      bv->end_value_ = end_value;
      bv->ramp_time_ = ramp_time;
      return std::move(bv);
    }

    std::unique_ptr<LinearBoundaryValue> CreateLinearBoundaryValue(double end_value, double ramp_time) {
      auto bv = std::make_unique<LinearBoundaryValue>();
      bv->end_value_ = end_value;
      bv->ramp_time_ = ramp_time;
      return std::move(bv);
    }

    std::unique_ptr<ExponentialBoundaryValue> CreateExponentialBoundaryValue(double end_value, double ramp_time) {
      auto bv = std::make_unique<ExponentialBoundaryValue>();
      bv->end_value_ = end_value;
      bv->step_time_ = ramp_time;
      return std::move(bv);
    }

    std::unique_ptr<PressureBoundaryCondition> CreatePressureInletBoundary(std::unique_ptr<IBoundaryValue> boundary_pressure, 
      const std::unordered_map<std::string, double>& yi, double Tin) {
      auto pressureBC = std::make_unique<PressureBoundaryCondition>();
      pressureBC->pressure_ = std::move(boundary_pressure);
      pressureBC->temperature_ = Tin;
      for (const auto& [component, y] : yi) {
        pressureBC->mole_fraction_[component] = y;
      }
      return std::move(pressureBC);
    }

    void AddStep(const std::string& step_name, double step_time) {
      data.cycle.AddStep(step_name);
      data.cycle.GetStep(step_name).step_time_ = step_time;
    }

    void ConfigureCycle(int number_of_cycles) {
      data.cycle.number_of_cycles_ = number_of_cycles;
    }

    // Configure KPI settings for cyclic tests.
    // feed_steps / product_steps: step indices (0-based)
    // feed_locs / product_locs: BoundaryLocation for each corresponding step
    // recycled_product_steps / recycled_product_locs: optional recycled product streams
    // product_comp_idx: index into the components vector for the product component
    void ConfigureKPI(
      int product_comp_idx,
      const std::vector<int>& feed_steps,
      const std::vector<BoundaryLocation>& feed_locs,
      const std::vector<int>& product_steps,
      const std::vector<BoundaryLocation>& product_locs,
      const std::vector<int>& recycled_product_steps = {},
      const std::vector<BoundaryLocation>& recycled_product_locs = {}
    ) {
      auto& kpi_cfg = data.kpi_config;
      kpi_cfg.SyncToComponents(data.fluid.components_);

      // Look up the per-component entry corresponding to product_comp_idx.
      if (product_comp_idx < 0 ||
          product_comp_idx >= static_cast<int>(data.fluid.components_.size())) {
        return;
      }
      const std::string& comp_name = data.fluid.components_[product_comp_idx];
      ComponentKPIConfig* kpi = kpi_cfg.FindByComponent(comp_name);
      if (!kpi) return;

      kpi->tracked = true;
      kpi->feed_endpoints.clear();
      kpi->product_endpoints.clear();
      for (size_t i = 0; i < feed_steps.size(); ++i)
        AppendEndpointsAt(kpi->feed_endpoints, feed_steps[i], feed_locs[i]);
      for (size_t i = 0; i < product_steps.size(); ++i)
        AppendEndpointsAt(kpi->product_endpoints, product_steps[i], product_locs[i]);
      // Recycled product streams are re-entry endpoints; fold them into the
      // product list so the signed sum subtracts them (matches the historical
      // "product - recycled product" recovery definition). Recycle is no longer
      // a separate role.
      for (size_t i = 0; i < recycled_product_steps.size(); ++i)
        AppendEndpointsAt(kpi->product_endpoints, recycled_product_steps[i], recycled_product_locs[i]);
    }

    // Append every KPI endpoint on `step`'s boundary at `loc` (all streams of a
    // flow BC, or the pressure boundary). Mirrors the old whole-face selection.
    void AppendEndpointsAt(std::vector<adsim::KPIEndpointRef>& out, int step, BoundaryLocation loc) {
      if (step < 0 || step >= static_cast<int>(data.cycle.step_list_.size())) return;
      for (const auto& ep : data.cycle.step_list_[step].GetEndpoints()) {
        if (ep.location == loc) out.push_back({step, ep.name, ep.location});
      }
    }

    // Boundary condition helper functions
    void AddFlowBC(std::vector<std::unique_ptr<IBoundaryStream>> streams, BoundaryLocation location, BoundaryDirection direction, const std::string& step_name) {
      auto flowBC = std::make_unique<FlowBoundaryCondition>();
      flowBC->location_ = location;
      flowBC->direction_ = direction;
      for (auto& stream : streams) {
        flowBC->streams_.push_back(std::move(stream));
      }
      Step& step = data.cycle.GetStep(step_name);
      if (location == BoundaryLocation::TOP) {
        step.top_boundary_condition_ = std::move(flowBC);
      } else {
        step.bottom_boundary_condition_ = std::move(flowBC);
      }
    }

    void AddClosedBC(BoundaryLocation location, const std::string& step_name) {
      auto closed = std::make_unique<ClosedBoundaryCondition>();
      closed->location_ = location;
      Step& step = data.cycle.GetStep(step_name);
      if (location == BoundaryLocation::TOP) {
        step.top_boundary_condition_ = std::move(closed);
      } else {
        step.bottom_boundary_condition_ = std::move(closed);
      }
    }

    void AddPressureOutletBC(std::unique_ptr<IBoundaryValue> boundary_pressure, BoundaryLocation location, const std::string& step_name, 
      const std::unordered_map<std::string, double>& mole_fractions = {}) {
      auto pressureBC = std::make_unique<PressureBoundaryCondition>();
      pressureBC->pressure_ = std::move(boundary_pressure);
      pressureBC->location_ = location;
      pressureBC->direction_ = BoundaryDirection::OUTLET;
      // Set mole fractions for backflow handling
      for (const auto& [component, y] : mole_fractions) {
        pressureBC->mole_fraction_[component] = y;
      }
      Step& step = data.cycle.GetStep(step_name);
      if (location == BoundaryLocation::TOP) {
        step.top_boundary_condition_ = std::move(pressureBC);
      } else {
        step.bottom_boundary_condition_ = std::move(pressureBC);
      }
    }

    void AddPressureInletBC(std::unique_ptr<IBoundaryValue> boundary_pressure, 
      const std::unordered_map<std::string, double>& yi, double Tin, 
      BoundaryLocation location, const std::string& step_name) {
      auto pressureBC = std::make_unique<PressureBoundaryCondition>();
      pressureBC->pressure_ = std::move(boundary_pressure);
      pressureBC->temperature_ = Tin;
      for (const auto& [component, y] : yi) {
        pressureBC->mole_fraction_[component] = y;
      }
      pressureBC->location_ = location;
      pressureBC->direction_ = BoundaryDirection::INLET;
      Step& step = data.cycle.GetStep(step_name);
      if (location == BoundaryLocation::TOP) {
        step.top_boundary_condition_ = std::move(pressureBC);
      } else {
        step.bottom_boundary_condition_ = std::move(pressureBC);
      }
    }

    void AddBC(std::unique_ptr<IBoundaryCondition> bc, BoundaryLocation location, const std::string& step_name) {
      bc->location_ = location;
      Step& step = data.cycle.GetStep(step_name);
      if (location == BoundaryLocation::TOP) {
        step.top_boundary_condition_ = std::move(bc);
      } else {
        step.bottom_boundary_condition_ = std::move(bc);
      }
    }

    bool RunSolver(std::string& message) {
      if (data.cycle.step_list_.empty()) {
        message = "No steps in cycle";
        return false;
      }

      CoupledSolver solver{data};
      std::string error_message{};
      int retval = solver.IntegrateCycle(data.cycle, false);
      
      bool is_successful = (retval >= 0);
      message = (is_successful ? "Success" : "Failed with code: " + std::to_string(retval));
      return is_successful;
    }

    void SaveResultsAsExample(const std::string& test_name, CoupledSolver& solver) {
      // Example deck names are attribute names of the Python module
      // skarstrom.examples (scripts/gen_examples_module.py), so they must be
      // identifiers: AuthorYear_Case, letters, digits and underscores.
      const bool identifier = !test_name.empty() && std::isalpha(static_cast<unsigned char>(test_name[0])) &&
                              std::all_of(test_name.begin(), test_name.end(), [](char c) {
                                return std::isalnum(static_cast<unsigned char>(c)) || c == '_';
                              });
      if (!identifier) {
        ADD_FAILURE() << "example deck name '" << test_name << "' is not a Python identifier";
        return;
      }
      std::filesystem::path dir_path = test_path / test_name;
      std::filesystem::path examples_path = std::filesystem::current_path() / "../examples/decks";
      
      std::filesystem::create_directories(examples_path);
      
      std::string json_filename = test_name + ".json";
      std::filesystem::path json_path = examples_path / json_filename;
      
      SaveToJson(data, json_path.string());    
    }

    void SaveTest(const std::string& test_name, CoupledSolver& solver) {
      // Create directory structure: Results/test_name/
      std::filesystem::path results_dir = test_path; // Already "Results"
      std::filesystem::path test_dir = results_dir / test_name;
      std::filesystem::path csv_results_dir = test_dir / "Results";
      
      std::filesystem::create_directories(results_dir);
      std::filesystem::create_directories(test_dir);
      std::filesystem::create_directories(csv_results_dir);
      
      // Save JSON file: Results/test_name/test_name.json
      std::string json_filename = test_name + ".json";
      std::filesystem::path json_path = test_dir / json_filename;
      
      SaveToJson(data, json_path.string());

      // Save CSV files: Results/test_name/Results/*.csv
      solver.SaveLogs(csv_results_dir.string());
    }

    void SaveTest(const std::string& test_name, CoupledSolverChuncked& solver) {
      // Create directory structure: Results/test_name/
      std::filesystem::path results_dir = test_path; // Already "Results"
      std::filesystem::path test_dir = results_dir / test_name;
      std::filesystem::path csv_results_dir = test_dir / "Results";
      
      std::filesystem::create_directories(results_dir);
      std::filesystem::create_directories(test_dir);
      std::filesystem::create_directories(csv_results_dir);
      
      // Save JSON file: Results/test_name/test_name.json
      std::string json_filename = test_name + ".json";
      std::filesystem::path json_path = test_dir / json_filename;
      
      SaveToJson(data, json_path.string());

      // Save CSV files: Results/test_name/Results/*.csv
      solver.SaveLogs(csv_results_dir.string());
    }

    // Save only the UserData configuration to Results/test_name/test_name.json
    // (the same path LoadSimulation reads), without any solver logs. Tests
    // must only ever load files saved under their own test_name — sharing
    // files between tests races under parallel ctest.
    void SaveSimulation(const std::string& test_name) {
      std::filesystem::path test_dir = test_path / test_name;
      std::filesystem::create_directories(test_dir);
      SaveToJson(data, (test_dir / (test_name + ".json")).string());
    }

    void LoadSimulation(const std::string& test_name) {
      // Load JSON file: Results/test_name/test_name.json
      std::filesystem::path results_dir = test_path; // Already "Results"
      std::filesystem::path test_dir = results_dir / test_name;
      std::string json_filename = test_name + ".json";
      std::filesystem::path json_path = test_dir / json_filename;
      
      std::ifstream json_file(json_path);
      if (!json_file.is_open()) {
        throw std::runtime_error("Failed to open JSON file: " + json_path.string());
      }
      
      cereal::JSONInputArchive archive(json_file);
      archive(data);
    }

    // Returns total moles of component accumulated in the bed (gas + solid phase) during a step.
    // Positive = net adsorption, negative = net desorption.
    double CalculateTotalAccumulated(const EquationLog& log, int cycle_index, int step_index, const std::string& component) const {
      double gas  = CalculateComponentAccumulation_Volume(log, cycle_index, step_index, ComponentVariableType::COMPONENT_MOLAR_CONCENTRATION, component);
      double solid = CalculateComponentAccumulation_Mass(log, cycle_index, step_index, ComponentVariableType::ADSORBED_CONCENTRATION, component);
      return gas + solid;
    }

    bool IsMassConserved(const EquationLog& log) const {
      for (int cycle_index : log.sys_data.cycle_index_set_) {
        for (int step_index : log.sys_data.step_index_set_) {
          auto yt = log.GetVariable(cycle_index, step_index, VariableType::TOTAL_MOLE_FRACTION);
          for (const auto& yt : yt.values) {
            for (int n = 1; n < yt.size() - 1; ++n) { // Dont include ghost cells!
              if (std::abs(yt[n] - 1) > 1e-3 || !std::isfinite(yt[n]))
                return false;
            }
          }          
        }
      }

      return true;
    }

    // Total-mole balance over each step: what crossed the two boundaries must
    // equal what accumulated in the gas and adsorbed phases.
    //
    // The boundary term is built by summing the per-component fluxes rather than
    // using GetTotalMolesIn. The two are not equivalent: GetTotalMolesIn
    // evaluates u*C with C taken from the ghost cell and ignores axial
    // dispersion, whereas GetTotalComponentMolesIn reproduces the discrete flux
    // the species equations actually use (upwinded convection plus the
    // dispersive flux at the boundary face). The cruder form does not close to a
    // useful tolerance, which is why this check previously logged a diagnostic
    // and then returned true regardless.
    bool IsMolesConserved(const EquationLog& log) const {
      auto& sys_data = log.sys_data;
      auto& components = sys_data.components_;
      auto& cycle_list = sys_data.cycle_index_set_;
      auto& step_list = sys_data.step_index_set_;
      constexpr double relative_tolerance_percent = 0.2;
      constexpr double absolute_tolerance_moles = 1e-5;

      for (int cycle_index : cycle_list) {
        for (int step_index : step_list) {
          double moles_in_bottom = 0.0;
          double moles_in_top = 0.0;
          double gas_moles_accumulated = 0.0;
          double solid_moles_accumulated = 0.0;
          for (const auto& comp : components) {
            moles_in_bottom += GetTotalComponentMolesIn(log, cycle_index, step_index, BoundaryLocation::BOTTOM, comp);
            moles_in_top += GetTotalComponentMolesIn(log, cycle_index, step_index, BoundaryLocation::TOP, comp);
            gas_moles_accumulated += CalculateComponentAccumulation_Volume(log, cycle_index, step_index, ComponentVariableType::COMPONENT_MOLAR_CONCENTRATION, comp);
            solid_moles_accumulated += CalculateComponentAccumulation_Mass(log, cycle_index, step_index, ComponentVariableType::ADSORBED_CONCENTRATION, comp);
          }
          double expected_moles_accumulated = moles_in_bottom + moles_in_top;

          double moles_accumulated = gas_moles_accumulated + solid_moles_accumulated;
          double error = std::abs(moles_accumulated - expected_moles_accumulated);
          // Use the larger of accumulated or expected to avoid division by very small numbers
          double denominator = std::max(std::abs(moles_accumulated), std::abs(expected_moles_accumulated));
          double percentage_error = (denominator > 1e-10) ? 100 * error / denominator : 0.0;
          // A through-flow step can move far more gas across the boundaries than
          // it accumulates (e.g. +34287 in, -34285 out, net 1.8 mol). The solver
          // integrates those flows only to ~reltol of the throughput, so also
          // accept an imbalance within 10x reltol (1e-5) of the throughput.
          const double throughput = std::abs(moles_in_bottom) + std::abs(moles_in_top);
          constexpr double throughput_relative_tolerance = 1e-5;
          if (error > absolute_tolerance_moles && percentage_error > relative_tolerance_percent
              && error > throughput_relative_tolerance * throughput) {
            std::cerr << "Mass conservation failed. Percentage error = " << percentage_error << "%, " << " (cycle " << cycle_index << ", step " << step_index << "). "
                      << "\nPercentage error = " << percentage_error << "%, "
                      << "\nAbsolute error = " << error << ", "
                      << "\nAccumulated = " << moles_accumulated << " (gas: " << gas_moles_accumulated << ", solid: " << solid_moles_accumulated << "), "
                      << "\nExpected = " << expected_moles_accumulated << " (bottom: " << moles_in_bottom << ", top: " << moles_in_top << ").\n\n";
            return false;
          }
        }
      }

      return true;
    }

    bool IsComponentMolesConserved(const EquationLog& log) const {
      auto& sys_data = log.sys_data;
      auto& components = sys_data.components_;
      auto& cycle_list = sys_data.cycle_index_set_;
      auto& step_list = sys_data.step_index_set_;
      constexpr double relative_tolerance_percent = 0.2;
      constexpr double absolute_tolerance_moles = 1e-5;
      
      for (auto component : components) {
        for (int cycle_index : cycle_list) {
          for (int step_index : step_list) {
            double moles_in_bottom = GetTotalComponentMolesIn(log, cycle_index, step_index, BoundaryLocation::BOTTOM, component);
            double moles_in_top = GetTotalComponentMolesIn(log, cycle_index, step_index, BoundaryLocation::TOP, component);
            double expected_moles_accumulated = moles_in_bottom + moles_in_top;
  
            double gas_moles_accumulated = CalculateComponentAccumulation_Volume(log, cycle_index, step_index, ComponentVariableType::COMPONENT_MOLAR_CONCENTRATION, component);
            double solid_moles_accumulated = CalculateComponentAccumulation_Mass(log, cycle_index, step_index, ComponentVariableType::ADSORBED_CONCENTRATION, component);
  
            double moles_accumulated = gas_moles_accumulated + solid_moles_accumulated;
            double error = std::abs(moles_accumulated - expected_moles_accumulated);
            // Use the larger of accumulated or expected to avoid division by very small numbers
            double denominator = std::max(std::abs(moles_accumulated), std::abs(expected_moles_accumulated));
            double percentage_error = (denominator > 1e-10) ? 100 * error / denominator : 0.0;
            // Same throughput criterion as IsMolesConserved.
            const double throughput = std::abs(moles_in_bottom) + std::abs(moles_in_top);
            constexpr double throughput_relative_tolerance = 1e-5;
            if (error > absolute_tolerance_moles && percentage_error > relative_tolerance_percent
                && error > throughput_relative_tolerance * throughput) {
                std::cerr << "Component mass conservation failed for " << component << " (cycle " << cycle_index << ", step " << step_index << "). "
                        << "\nPercentage error = " << percentage_error << "%, "
                        << "\nAbsolute error = " << error << ", "
                        << "\nAccumulated = " << moles_accumulated << " (gas: " << gas_moles_accumulated << ", solid: " << solid_moles_accumulated << "), "
                        << "\nExpected = " << expected_moles_accumulated << " (bottom: " << moles_in_bottom << ", top: " << moles_in_top << ").\n\n";
              return false;
            }   
          }
        }
      }

      return true;
    }

    // Resolve a relative data path (e.g. "data/ward_2022/Qin=100") against the
    // filesystem regardless of the current working directory. The digitised data
    // lives under <repo>/adsim/data/, so walk up from the CWD trying both the
    // path as-given and under an "adsim/" prefix. This lets the tests find their
    // data whether run via ctest (WORKING_DIRECTORY=adsim), from the repo root,
    // or from a build directory. Falls back to the original path (which then
    // errors informatively) if nothing matches.
    static std::filesystem::path ResolveDataPath(const std::string& relative) {
      namespace fs = std::filesystem;
      for (fs::path dir = fs::current_path(); ; dir = dir.parent_path()) {
        if (fs::exists(dir / relative))            return dir / relative;
        if (fs::exists(dir / "adsim" / relative))  return dir / "adsim" / relative;
        if (dir == dir.parent_path()) break; // reached filesystem root
      }
      return relative;
    }

    void SaveExperimentalData(const std::string& path) {
      data.exp_data_.DeleteAllData();

      const std::filesystem::path resolved = ResolveDataPath(path);
      std::vector<std::filesystem::path> csv_files;
      for (const auto& entry : std::filesystem::directory_iterator(resolved)) {
        if (entry.is_regular_file() && entry.path().extension() == ".csv") {
          csv_files.push_back(entry.path());
        }
      }

      std::sort(csv_files.begin(), csv_files.end());

      for (const auto& file : csv_files) {
        data.exp_data_.LoadData(file.string(), true);
      }
    }

    // Helper to extract outlet mole fraction from EquationLog
    std::pair<std::vector<double>, std::vector<double>> GetOutletMoleFraction(
      const CoupledSolver& solver, 
      const std::string& component,
      int cycle_index,
      int step_index,
      BoundaryLocation loc) const 
    {
      const auto& log = solver.GetLog();
      const auto& mole_fraction_entry = log.GetComponentFaceVariable(cycle_index, step_index,
        ComponentVariableType::MOLE_FRACTION, component);
      
      if (mole_fraction_entry.values.empty() || mole_fraction_entry.times.empty()) {
        return {{}, {}};
      }
      
      size_t cell_idx = 0;
      if (loc == BoundaryLocation::TOP) {
        cell_idx = mole_fraction_entry.values[0].size() - 1;
      }
      
      std::vector<double> times;
      std::vector<double> mole_fractions;
      
      for (size_t i = 0; i < mole_fraction_entry.times.size(); ++i) {
        if (i < mole_fraction_entry.values.size() && 
            cell_idx < mole_fraction_entry.values[i].size()) {
          times.push_back(mole_fraction_entry.times[i]);
          mole_fractions.push_back(mole_fraction_entry.values[i][cell_idx]);
        }
      }

      return {times, mole_fractions};
    }

    double CalculateExperimentalError(
      const std::vector<double>& sim_time,
      const std::vector<double>& sim_var,
      const std::vector<double>& exp_time,
      const std::vector<double>& exp_var
    ) {
      // Basic validation.
      if (sim_time.size() != sim_var.size() || exp_time.size() != exp_var.size()) {
        throw std::invalid_argument("Time and value vectors must have the same length for each series.");
      }
      if (sim_time.size() < 2) {
        throw std::invalid_argument("Simulation series must have at least two points for interpolation.");
      }
      if (exp_time.empty()) {
        throw std::invalid_argument("Experimental series is empty.");
      }
      if (!is_strictly_increasing(sim_time)) {
        throw std::invalid_argument("sim_time must be strictly increasing (no duplicates).");
      }
      if (!is_strictly_increasing(exp_time)) {
        throw std::invalid_argument("exp_time must be strictly increasing (no duplicates).");
      }

      const double tmin = sim_time.front();
      const double tmax = sim_time.back();

      double sum_sq = 0.0;
      size_t count = 0;

      auto start_it = std::lower_bound(exp_time.begin(), exp_time.end(), tmin);
      auto end_it   = std::upper_bound(exp_time.begin(), exp_time.end(), tmax);

      for (auto it = start_it; it != end_it; ++it) {
        size_t i = static_cast<size_t>(std::distance(exp_time.begin(), it));
        double t = *it;

        // Interpolate simulation at experimental time t.
        double sim_y = interp_linear(sim_time, sim_var, t);
        double e = sim_y - exp_var[i];
        sum_sq += e * e;
        ++count;
      }

      // If there is no overlap, return NaN to signal the issue.
      if (count == 0) {
        return std::numeric_limits<double>::quiet_NaN();
      }

      // RMSE
      double mse = sum_sq / static_cast<double>(count);
      return std::sqrt(mse);
    }

  };

} // end adsim namespace
