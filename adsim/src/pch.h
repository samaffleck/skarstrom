#pragma once

// Cereal includes
#include <cereal/cereal.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/archives/json.hpp>

// Eigen includes
#include <Eigen/Dense>
#include <Eigen/SparseLU>

// STL includes
#include <assert.h>
#include <stdexcept>
#include <vector>
#include <memory>
#include <string>
#include <unordered_map>
#include <cmath> 
#include <algorithm>
#include <atomic>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <chrono>
#include <iostream>
#include <charconv>

// CoolProp includes
#include <CoolProp.h>

// SUNDIALS includes
#include <ida/ida.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>


