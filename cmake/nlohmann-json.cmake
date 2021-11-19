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



# FetchContent_Declare(
#     nlohmann_json
#     URL "https://github.com/nlohmann/json/releases/download/${NLOHMANNJSON_VERSION}/include.zip"
#     URL_HASH SHA256=6bea5877b1541d353bd77bdfbdb2696333ae5ed8f9e8cc22df657192218cad91
# )
# FetchContent_MakeAvailable(nlohmann_json)

#add_library(nlohmann_json INTERFACE)
#target_include_directories(nlohmann_json INTERFACE
#    ${nlohmann_json_SOURCE_DIR}/include)
#add_library(nlohmann::json ALIAS nlohmann_json)
