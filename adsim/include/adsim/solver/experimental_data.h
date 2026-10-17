/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// STL includes 
#include <string>
#include <vector>
#include <unordered_map>
#include <filesystem>

// adsim includes
#include "adsim/objects/reactor.h"


namespace adsim {

	class XYData {
	public:
		std::string name;
		std::vector<double> x;
		std::vector<double> y;
		std::string x_label;
		std::string y_label;

	private:
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(
				cereal::make_nvp("Name", name),
				cereal::make_nvp("x_data", x),
				cereal::make_nvp("y_data", y),
				cereal::make_nvp("x_label", x_label),
				cereal::make_nvp("y_label", y_label)
			);
		}
	};

  class ExperimentalData {
  public:
		bool LoadData(const std::filesystem::path& path, bool headers_included);
		bool AddData(const XYData& data);
		bool DeleteDataByIndex(size_t idx);
		void DeleteAllData() { data_.clear(); }		
		std::vector<XYData>& GetAllData() { return data_; }
		size_t size() const { return data_.size(); }

	private:
		std::vector<XYData> data_;
		
		friend class cereal::access;
		template <class Archive>
		void serialize(Archive& ar, std::uint32_t const version) {
			ar(
				cereal::make_nvp("Data", data_)
			);
		}

  };

} // end adsim namespace
