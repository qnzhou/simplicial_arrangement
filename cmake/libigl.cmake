include_guard()

include(FetchContent)
FetchContent_Declare(
    libigl
    GIT_REPOSITORY https://github.com/libigl/libigl.git
    GIT_TAG jdumas/cmake-refactor
)
FetchContent_GetProperties(libigl)
if(libigl_POPULATED)
    return()
endif()
FetchContent_MakeAvailable(libigl)
