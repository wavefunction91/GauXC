include(FetchContent)
include(gauxc-dep-versions)

FetchContent_Declare(linalg-cmake-modules
  GIT_REPOSITORY ${GAUXC_LINALG_MODULES_REPOSITORY}
  GIT_TAG        ${GAUXC_LINALG_MODULES_REVISION}
)

FetchContent_MakeAvailable(linalg-cmake-modules)

# Ensure the module path includes the fetched content
list(PREPEND CMAKE_MODULE_PATH ${linalg-cmake-modules_SOURCE_DIR})

