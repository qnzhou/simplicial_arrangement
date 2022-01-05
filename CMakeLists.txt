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

option(LIBIGL_WITH_CGAL "Build CGAL support in libigl" On)
include(libigl)

add_library(arrangement STATIC
    ${CMAKE_SOURCE_DIR}/src/PyMesh/Arrangement.cpp
    ${CMAKE_SOURCE_DIR}/src/PyMesh/FastArrangement.cpp
    ${CMAKE_SOURCE_DIR}/src/PyMesh/MeshArrangement.cpp)
target_link_libraries(arrangement PUBLIC igl::core igl::cgal
    fast_arrangement::fast_arrangement)

file(GLOB INC_FILES ${CMAKE_CURRENT_LIST_DIR}/include/*.h)
file(GLOB SRC_FILES ${CMAKE_CURRENT_LIST_DIR}/src/*.cpp)

add_library(simplicial_arrangement STATIC ${SRC_FILES} ${INC_FILES})
target_include_directories(simplicial_arrangement PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/include)
target_link_libraries(simplicial_arrangement
    PUBLIC absl::numeric absl::flat_hash_map nlohmann_json::nlohmann_json
    PRIVATE implicit_predicates::implicit_predicates spdlog::spdlog)
target_compile_features(simplicial_arrangement PRIVATE cxx_std_17)
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
endif()




set(LOOKUP_TABLE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/lookup_table")
target_compile_definitions(simplicial_arrangement PUBLIC -DLOOKUP_TABLE_PATH="${LOOKUP_TABLE_PATH}")

add_library(impl_arrangement_util STATIC "app/implicit_arrangement_util.cpp")
target_compile_features(impl_arrangement_util PRIVATE cxx_std_17)
#target_include_directories(impl_arrangement_util PUBLIC ${CMAKE_CURRENT_LIST_DIR}/include)
target_link_libraries(impl_arrangement_util PUBLIC simplicial_arrangement::simplicial_arrangement)

add_executable(impl_arrangement "app/implicit_arrangement.cpp")
target_compile_features(impl_arrangement PRIVATE cxx_std_17)
target_link_libraries(impl_arrangement simplicial_arrangement::simplicial_arrangement impl_arrangement_util)
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
target_link_libraries(mesh_arrangement impl_arrangement_util igl::core arrangement)
target_compile_features(mesh_arrangement PRIVATE cxx_std_17)
if (SIMPLICIAL_ARRANGEMENT_TEST_COVERAGE)
    if (${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
        target_compile_options(mesh_arrangement PRIVATE -fprofile-instr-generate -fcoverage-mapping)
        target_link_options(mesh_arrangement PRIVATE -fprofile-instr-generate -fcoverage-mapping)
    endif()
endif()
