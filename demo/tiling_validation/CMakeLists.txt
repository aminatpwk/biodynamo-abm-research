cmake_minimum_required(VERSION 3.19.3)

project(tiling_validation)

find_package(BioDynaMo REQUIRED)
include("${BDM_USE_FILE}")
include_directories("src")

# tiling_experiment reuses the epidemiology demo's Person + behaviors + SimParam,
# so that demo's src dir is added to the include path.
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../epidemiology/src")

# Sources and dictionary headers are listed explicitly PER TARGET on purpose:
#  * bdm_add_executable compiles all SOURCES into one shared library and
#    generates its own main() stub, so each target must contain exactly one
#    translation unit with a main(). Globbing all src/*.cc into both targets
#    would put two main()s in each library and fail to link.
#  * Keeping the header lists disjoint keeps the epidemiology (and GSL-adjacent)
#    headers out of the validation target's ROOT dictionary entirely.
set(VALIDATION_HEADERS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/tile_extraction.h")

set(EXPERIMENT_HEADERS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/tile_extraction.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/tiling_config.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../epidemiology/src/person.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../epidemiology/src/behavior.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../epidemiology/src/sim-param.h")

bdm_add_executable(tiling_validation
                   HEADERS "${VALIDATION_HEADERS}"
                   SOURCES src/tiling_validation.cc
                   LIBRARIES "${BDM_REQUIRED_LIBRARIES}")

bdm_add_executable(tiling_experiment
                   HEADERS "${EXPERIMENT_HEADERS}"
                   SOURCES src/tiling_experiment.cc
                   LIBRARIES "${BDM_REQUIRED_LIBRARIES}")
