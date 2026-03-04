cmake_minimum_required (VERSION 3.22)

project(dsa
  LANGUAGES CXX
  VERSION 0.0.0
)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(ENABLE_DEBUG_FLAGS "Enable extra debug compile options" OFF)
option(MAKE_TEST "Enable compilation of test binary" OFF)

# user paths to htslib if needed
set(HTSLIB_INCLUDE_DIR "" CACHE PATH "Path to directory that contains htslib/hts.h")
set(HTSLIB_LIBRARY "" CACHE FILEPATH "Absolute path to libhts.{so,dylib,a}")

# //--- find deps ---//
find_package(PkgConfig QUIET)
find_package(OpenSSL REQUIRED)

pkg_check_modules(LIBDEFLATE REQUIRED IMPORTED_TARGET libdeflate)

# ------- htslib -------
# NOTE: no version control
message(STATUS "Finding htslib")
set(HTSLIB_TARGET "")
set(_htslib_inc "")
set(_htslib_lib "")

if(HTSLIB_INCLUDE_DIR AND HTSLIB_LIBRARY)
  message(STATUS "Using user-specified htslib install")
  set(_htslib_inc "${HTSLIB_INCLUDE_DIR}")
  set(_htslib_lib "${HTSLIB_LIBRARY}")
endif()

if(_htslib_inc AND _htslib_lib)
  add_library(htslib UNKNOWN IMPORTED GLOBAL)
  add_library(htslib::htslib ALIAS htslib)
  set_target_properties(htslib PROPERTIES
    IMPORTED_LOCATION "${_htslib_lib}"
    INTERFACE_INCLUDE_DIRECTORIES "${_htslib_inc}"
  )
  set(HTSLIB_TARGET htslib::htslib)
endif()

if(NOT HTSLIB_TARGET AND PkgConfig_FOUND)
  # assuming the pkg-config module is called "htslib"
  pkg_check_modules(HTSLIB QUIET IMPORTED_TARGET htslib)
  if(HTSLIB_FOUND)
    message(STATUS "Located htslib install via pkg-config.")
    set(HTSLIB_TARGET PkgConfig::HTSLIB)
  endif()
endif()

if(NOT HTSLIB_TARGET)
  message(FATAL_ERROR
    "Could not find htslib directly or via pkg-config. "
    "Set both HTSLIB_INCLUDE_DIR and HTSLIB_LIBRARY to use user install."
  )
endif()
# //--- end deps ---//

# //--- testing ---//
if(MAKE_TEST)
  message(STATUS "Finding catch2 for test binary")
  Include(FetchContent)

  FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v3.8.1
  )

  FetchContent_MakeAvailable(Catch2)

  add_executable(test-dsa tests/test.cc)
  target_include_directories(test-dsa PRIVATE src)
  target_link_libraries(test-dsa PRIVATE
    Catch2::Catch2WithMain
  )
endif(MAKE_TEST)
# //--- end testing ---//

add_executable(dsa
  src/ref.cc
  src/mask.cc
  src/mask_loader.cc
  src/pileup_custom.cc
  src/pileup.cc
  src/pileup_batch.cc
  src/dsa.cc
)

target_compile_options(dsa PRIVATE
    -Wall
    -Wextra
    -Wpedantic
    -Wshadow
    -Wfloat-equal
    -Wnull-dereference
    -Warray-bounds
    -Wformat=2
)

if(ENABLE_DEBUG_FLAGS)
  target_compile_options(dsa PRIVATE
    -Wconversion
    -Wsign-conversion
    -Wunused
    -Wold-style-cast
    -fsanitize=address,undefined
    -fno-omit-frame-pointer
    -O0 -g
  )
  target_link_options(dsa PRIVATE
    -fsanitize=address,undefined
  )
endif()

target_include_directories(
  dsa PRIVATE src
)

target_link_libraries(dsa
  PRIVATE
  ${HTSLIB_TARGET}
  PkgConfig::LIBDEFLATE
  OpenSSL::SSL
  OpenSSL::Crypto
)

