if(TARGET fast_arrangement::fast_arrangement)
    return()
endif()

include(FetchContent)
include(indirect_predicates)

# cinolib
# It is a large library and it seems we do not need to build it directly.
FetchContent_Declare(
    cinolib
    GIT_REPOSITORY https://github.com/mlivesu/cinolib.git
    GIT_TAG v1.0
)
if(NOT cinolib_POPULATED)
  FetchContent_Populate(cinolib)
endif()

# arrangement
FetchContent_Declare(
    fast_arrangement
    GIT_REPOSITORY https://github.com/qnzhou/FastAndRobustMeshArrangements.git
    GIT_TAG master
)
if(NOT fast_arrangement_POPULATED)
  FetchContent_Populate(fast_arrangement)
endif()

add_library(fast_arrangement
    ${cinolib_SOURCE_DIR}/external/predicates/shewchuk.c)
target_link_libraries(fast_arrangement
    PUBLIC Eigen3::Eigen indirect_predicates::indirect_predicates)
target_include_directories(fast_arrangement
    PUBLIC
    ${cinolib_SOURCE_DIR}/include
    ${fast_arrangement_SOURCE_DIR}/code
    )
#target_compile_options(fast_arrangement PRIVATE -O0)
target_compile_definitions(fast_arrangement PUBLIC CINOLIB_USES_EXACT_PREDICATES)
if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    target_compile_definitions(fast_arrangement PUBLIC -DIS64BITPLATFORM)
endif()

add_library(fast_arrangement::fast_arrangement ALIAS fast_arrangement)