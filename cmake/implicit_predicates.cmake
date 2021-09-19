if (TARGET implicit_predicates::implicit_predicates)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
    implicit_predicates
    GIT_REPOSITORY git@github.com:qnzhou/implicit_predicates.git
    GIT_TAG main
    )

FetchContent_MakeAvailable(implicit_predicates)
