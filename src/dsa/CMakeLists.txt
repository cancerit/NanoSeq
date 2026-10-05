cmake_minimum_required (VERSION 3.22)

project(dsa
  LANGUAGES CXX
  VERSION 0.0.0
)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(COMPILE_COMMANDS "Export compile_commands.json" OFF)
option(ENABLE_DEBUG_FLAGS "Enable extra debug compile options" OFF)
option(MAKE_TEST "Enable compilation of test binary" OFF)

if(COMPILE_COMMANDS)
  set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
endif()

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

set(DSA_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
include(${DSA_SRC_DIR}/sources.cmake)
list(TRANSFORM DSA_SOURCES PREPEND ${DSA_SRC_DIR}/)

# //--- testing ---//
if(MAKE_TEST)
  add_executable(test-dsa ${DSA_SRC_DIR}/tests/test.cc)
  target_include_directories(test-dsa PRIVATE
    ${DSA_SRC_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../../vendor/doctest-2.4.11
  )
endif(MAKE_TEST)
# //--- end testing ---//

add_executable(dsa ${DSA_SOURCES})

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
  dsa PRIVATE ${DSA_SRC_DIR}
)

target_link_libraries(dsa
  PRIVATE
  ${HTSLIB_TARGET}
  PkgConfig::LIBDEFLATE
  OpenSSL::SSL
  OpenSSL::Crypto
)

