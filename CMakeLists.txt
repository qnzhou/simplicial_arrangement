cmake_minimum_required(VERSION 3.14)

project(simplicial_arrangement LANGUAGES CXX C)

include(FetchContent)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)
include(absl)
include(implicit_predicates)
include(spdlog)
include(nlohmann-json)
include(sanitizers)
include(fast_arrangement)
include(cli11)
include(mshio)


option(LIBIGL_COPYLEFT_WITH_CGAL "Build CGAL support in libigl" On)
include(libigl)

add_library(arrangement STATIC
    ${CMAKE_SOURCE_DIR}/src/PyMesh/Arrangement.cpp
    ${CMAKE_SOURCE_DIR}/src/PyMesh/FastArrangement.cpp
    ${CMAKE_SOURCE_DIR}/src/PyMesh/MeshArrangement.cpp)
target_link_libraries(arrangement PUBLIC igl::core igl_copyleft::cgal
    fast_arrangement::fast_arrangement)

option(SIMPLICIAL_ARRANGEMENT_LOOKUP_TABLE_GEN
    "Lookup table generation mode (for developers only)" OFF)
option(SIMPLICIAL_ARRANGEMENT_DETERMINISTIC
    "Reqruire deterministic result (slightly slower)" OFF)
option(SIMPLICIAL_ARRANGEMENT_NON_ROBUST
    "Use non-robust predicates (for developers only)" ON)

file(GLOB INC_FILES ${CMAKE_CURRENT_LIST_DIR}/include/*.h)
file(GLOB SRC_FILES ${CMAKE_CURRENT_LIST_DIR}/src/*.cpp)

add_library(simplicial_arrangement STATIC ${SRC_FILES} ${INC_FILES})
target_include_directories(simplicial_arrangement PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/include)
target_link_libraries(simplicial_arrangement
    PUBLIC absl::numeric absl::flat_hash_map nlohmann_json::nlohmann_json
    PRIVATE implicit_predicates::implicit_predicates spdlog::spdlog)
target_compile_features(simplicial_arrangement PRIVATE cxx_std_17)
if (SIMPLICIAL_ARRANGEMENT_LOOKUP_TABLE_GEN)
    target_compile_definitions(simplicial_arrangement PUBLIC -DLOOKUP_TABLE_GEN)
endif()
if (SIMPLICIAL_ARRANGEMENT_DETERMINISTIC)
    target_compile_definitions(simplicial_arrangement PRIVATE
        -DSIMPLICIAL_ARRANGEMENT_DETERMINISTIC)
endif()
if (SIMPLICIAL_ARRANGEMENT_NON_ROBUST)
    target_compile_definitions(simplicial_arrangement PRIVATE
        -DSIMPLICIAL_ARRANGEMENT_NON_ROBUST)
endif()
if (NOT MSVC)
    target_compile_options(simplicial_arrangement PRIVATE "-Wall" "-Werror")
endif()
add_library(simplicial_arrangement::simplicial_arrangement ALIAS
    simplicial_arrangement)

option(SIMPLICIAL_ARRANGEMENT_UNIT_TESTS "Build unit tests for simplicial arrangement" OFF)
if (SIMPLICIAL_ARRANGEMENT_UNIT_TESTS)
    include(CTest)
    enable_testing()
    include(Catch2)

    file(GLOB TEST_FILES "${CMAKE_CURRENT_LIST_DIR}/tests/*.cpp")
    add_executable(simplicial_arrangement_tests ${TEST_FILES})
    target_include_directories(simplicial_arrangement_tests PRIVATE
        "${CMAKE_CURRENT_LIST_DIR}/src/")
    target_link_libraries(simplicial_arrangement_tests PRIVATE
        simplicial_arrangement::simplicial_arrangement
        implicit_predicates::implicit_predicates Catch2::Catch2 spdlog::spdlog)
    target_compile_features(simplicial_arrangement_tests PRIVATE cxx_std_17)
    if (NOT MSVC)
        target_compile_options(simplicial_arrangement_tests PRIVATE "-Wall" "-Werror")
    else()
        target_compile_definitions(simplicial_arrangement_tests PRIVATE "-D_USE_MATH_DEFINES")
    endif()
    target_compile_definitions(simplicial_arrangement_tests PRIVATE CATCH_CONFIG_ENABLE_BENCHMARKING)
    catch_discover_tests(simplicial_arrangement_tests)

    option(SIMPLICIAL_ARRANGEMENT_TEST_COVERAGE "Track unit test coverage" OFF)
    if (SIMPLICIAL_ARRANGEMENT_TEST_COVERAGE)
        if (${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
            target_compile_options(simplicial_arrangement PRIVATE -fprofile-instr-generate -fcoverage-mapping)
            target_link_options(simplicial_arrangement PRIVATE -fprofile-instr-generate -fcoverage-mapping)
            target_compile_options(simplicial_arrangement_tests PRIVATE -fprofile-instr-generate -fcoverage-mapping)
            target_link_options(simplicial_arrangement_tests PRIVATE -fprofile-instr-generate -fcoverage-mapping)
        endif()
    endif()
    if (SIMPLICIAL_ARRANGEMENT_NON_ROBUST)
        target_compile_definitions(simplicial_arrangement_tests PRIVATE
            -DSIMPLICIAL_ARRANGEMENT_NON_ROBUST)
    endif()
endif()




set(LOOKUP_TABLE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/lookup_table")
target_compile_definitions(simplicial_arrangement PUBLIC -DLOOKUP_TABLE_PATH="${LOOKUP_TABLE_PATH}")

add_library(impl_arrangement_util STATIC "app/implicit_arrangement_util.cpp")
target_compile_features(impl_arrangement_util PRIVATE cxx_std_17)
#target_include_directories(impl_arrangement_util PUBLIC ${CMAKE_CURRENT_LIST_DIR}/include)
target_link_libraries(impl_arrangement_util PUBLIC
    simplicial_arrangement::simplicial_arrangement igl::core mshio::mshio)

add_executable(impl_arrangement "app/implicit_arrangement.cpp")
target_compile_features(impl_arrangement PRIVATE cxx_std_17)
target_link_libraries(impl_arrangement simplicial_arrangement::simplicial_arrangement impl_arrangement_util igl::core CLI11::CLI11)
if (SIMPLICIAL_ARRANGEMENT_TEST_COVERAGE)
    if (${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
        target_compile_options(impl_arrangement PRIVATE -fprofile-instr-generate -fcoverage-mapping)
        target_link_options(impl_arrangement PRIVATE -fprofile-instr-generate -fcoverage-mapping)
    endif()
endif()

add_executable(marching_tet "app/marching_tet.cpp")
target_link_libraries(marching_tet simplicial_arrangement::simplicial_arrangement impl_arrangement_util)
target_compile_features(marching_tet PRIVATE cxx_std_17)
if (SIMPLICIAL_ARRANGEMENT_TEST_COVERAGE)
    if (${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
        target_compile_options(marching_tet PRIVATE -fprofile-instr-generate -fcoverage-mapping)
        target_link_options(marching_tet PRIVATE -fprofile-instr-generate -fcoverage-mapping)
    endif()
endif()

add_executable(mesh_arrangement "app/mesh_arrangement.cpp")
target_include_directories(mesh_arrangement PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(mesh_arrangement impl_arrangement_util igl::core arrangement CLI11::CLI11)
target_compile_features(mesh_arrangement PRIVATE cxx_std_17)
if (SIMPLICIAL_ARRANGEMENT_TEST_COVERAGE)
    if (${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
        target_compile_options(mesh_arrangement PRIVATE -fprofile-instr-generate -fcoverage-mapping)
        target_link_options(mesh_arrangement PRIVATE -fprofile-instr-generate -fcoverage-mapping)
    endif()
endif()

add_executable(material_interface "app/material_interface.cpp")
target_compile_features(material_interface PRIVATE cxx_std_17)
target_link_libraries(material_interface simplicial_arrangement::simplicial_arrangement impl_arrangement_util igl::core CLI11::CLI11)
if (SIMPLICIAL_ARRANGEMENT_TEST_COVERAGE)
    if (${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
        target_compile_options(material_interface PRIVATE -fprofile-instr-generate -fcoverage-mapping)
        target_link_options(material_interface PRIVATE -fprofile-instr-generate -fcoverage-mapping)
    endif()
endif()

add_executable(dual_contouring "app/dual_contouring.cpp")
target_compile_features(dual_contouring PRIVATE cxx_std_17)
target_link_libraries(dual_contouring impl_arrangement_util igl::core CLI11::CLI11)
if (SIMPLICIAL_ARRANGEMENT_TEST_COVERAGE)
    if (${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
        target_compile_options(dual_contouring PRIVATE -fprofile-instr-generate -fcoverage-mapping)
        target_link_options(dual_contouring PRIVATE -fprofile-instr-generate -fcoverage-mapping)
    endif()
endif()

add_executable(dual_contouring_MI "app/dual_contouring_MI.cpp")
target_compile_features(dual_contouring_MI PRIVATE cxx_std_17)
target_link_libraries(dual_contouring_MI impl_arrangement_util igl::core CLI11::CLI11)
if (SIMPLICIAL_ARRANGEMENT_TEST_COVERAGE)
    if (${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
        target_compile_options(dual_contouring_MI PRIVATE -fprofile-instr-generate -fcoverage-mapping)
        target_link_options(dual_contouring_MI PRIVATE -fprofile-instr-generate -fcoverage-mapping)
    endif()
endif()

add_executable(generate_material_interface_lookup_table
    "app/generate_material_interface_lookup_table.cpp")
target_link_libraries(generate_material_interface_lookup_table PRIVATE
    simplicial_arrangement::simplicial_arrangement nlohmann_json::nlohmann_json
    spdlog::spdlog)
target_compile_features(generate_material_interface_lookup_table PRIVATE cxx_std_17)
