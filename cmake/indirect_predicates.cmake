if (TARGET indirect_predicates::indirect_predicates)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
    indirect_predicates
    GIT_REPOSITORY https://github.com/MarcoAttene/Indirect_Predicates.git
    GIT_TAG ca74afb8c664ed24b4482fb3c8a3175102ce3dcb
    )

FetchContent_Populate(indirect_predicates)

SET(INDIRECT_PREDICATES_INC_FILES
    ${indirect_predicates_SOURCE_DIR}/numerics.h
    ${indirect_predicates_SOURCE_DIR}/implicit_point.h
    ${indirect_predicates_SOURCE_DIR}/indirect_predicates.h)

add_library(indirect_predicates INTERFACE)
#add_library(indirect_predicates ${INDIRECT_PREDICATES_SRC_FILES} ${INDIRECT_PREDICATES_INC_FILES})
target_include_directories(indirect_predicates INTERFACE ${indirect_predicates_SOURCE_DIR}/include)
target_compile_features(indirect_predicates INTERFACE cxx_std_20)

# Compiler-specific options
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # grant IEEE 754 compliance
    target_compile_options(indirect_predicates INTERFACE "-frounding-math"
        "-Wno-format-security")
    # use intrinsic functions (CHECK WHAT TO DO FOR GCC !!!!!!!!)
    # target_compile_options(indirect_predicates PUBLIC "/Oi")
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    # grant IEEE 754 compliance
    target_compile_options(indirect_predicates INTERFACE "/fp:strict")
    # use intrinsic functions
    target_compile_options(indirect_predicates INTERFACE "/Oi")
    # turn off annoying warnings
    target_compile_options(indirect_predicates INTERFACE "/D _CRT_SECURE_NO_WARNINGS")
elseif (${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
    # From the repo README:
    # https://github.com/MarcoAttene/Indirect_Predicates#system-requirements
    #
    # Mac OSX. Unfortunately CLANG does not support direct access to the
    # floating point environment. The only way to use this software with CLANG
    # is to disable all optimizations (-O0).
    target_compile_options(indirect_predicates INTERFACE "-O0"
        "-Wno-format-security" "-Wno-unknown-pragmas" "-Wno-return-type")
endif()

add_library(indirect_predicates::indirect_predicates ALIAS indirect_predicates)
