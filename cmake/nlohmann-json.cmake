#[[
if(TARGET nlohmann_json::nlohmann_json)
    return()
endif()

set(NLOHMANNJSON_VERSION "v3.10.4")

include(FetchContent)
FetchContent_Declare(
  nlohmann_json
  GIT_REPOSITORY https://github.com/ArthurSonzogni/nlohmann_json_cmake_fetchcontent
  GIT_TAG v3.10.4)
#FetchContent_MakeAvailable(nlohmann_json)

FetchContent_GetProperties(nlohmann_json)
if(NOT json_POPULATED)
  FetchContent_Populate(nlohmann_json)
  add_subdirectory(${nlohmann_json_SOURCE_DIR} ${nlohmann_json_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

FetchContent_MakeAvailable(nlohmann_json)

]]

#
# Copyright 2020 Adobe. All rights reserved.
# This file is licensed to you under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License. You may obtain a copy
# of the License at http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under
# the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR REPRESENTATIONS
# OF ANY KIND, either express or implied. See the License for the specific language
# governing permissions and limitations under the License.
#
if(TARGET nlohmann_json::nlohmann_json)
    return()
endif()

message(STATUS "Third-party (external): creating target 'nlohmann_json::nlohmann_json'")

# nlohmann_json is a big repo for a single header, so we just download the release archive
set(NLOHMANNJSON_VERSION "v3.7.3")

include(FetchContent)
FetchContent_Declare(
    nlohmann_json
    URL "https://github.com/nlohmann/json/releases/download/${NLOHMANNJSON_VERSION}/include.zip"
    URL_HASH SHA256=87b5884741427220d3a33df1363ae0e8b898099fbc59f1c451113f6732891014
)
FetchContent_MakeAvailable(nlohmann_json)

add_library(nlohmann_json INTERFACE)
add_library(nlohmann_json::nlohmann_json ALIAS nlohmann_json)

include(GNUInstallDirs)
target_include_directories(nlohmann_json INTERFACE
    "$<BUILD_INTERFACE:${nlohmann_json_SOURCE_DIR}>/include"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

# Let the code know that we have nlohmann json
# TODO: This shouldn't be here (should be set by the caller)
target_compile_definitions(nlohmann_json INTERFACE -DLA_WITH_NLOHMANNJSON)

# Install rules
set(CMAKE_INSTALL_DEFAULT_COMPONENT_NAME nlohmann_json)
install(DIRECTORY ${nlohmann_json_SOURCE_DIR}/include/nlohmann DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(TARGETS nlohmann_json EXPORT NlohmannJson_Targets)
install(EXPORT NlohmannJson_Targets DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/nlohmann_json NAMESPACE nlohmann_json::)
export(EXPORT NlohmannJson_Targets FILE "${CMAKE_CURRENT_BINARY_DIR}/NlohmannJsonTargets.cmake")

