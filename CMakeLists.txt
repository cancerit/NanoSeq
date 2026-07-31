cmake_minimum_required(VERSION 3.18)

# Set handling of downloaded content.
# Fetched archives should get extraction-time timestamps, not the archive's
# original ones, so a changed URL is correctly seen as new content.
if(POLICY CMP0135)
  cmake_policy(SET CMP0135 NEW)
endif()

# ==============================================================================
# Version
# ==============================================================================
# NANOSEQ_VERSION is set from `git describe` when building from a git
# checkout, and left empty otherwise (e.g. building from a release tarball).
find_package(Git QUIET)
set(NANOSEQ_VERSION "")
if(Git_FOUND)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --tags --always --dirty
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE NANOSEQ_VERSION_RAW
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE NANOSEQ_GIT_RESULT
    ERROR_QUIET)
  if(NANOSEQ_GIT_RESULT EQUAL 0)
    set(NANOSEQ_VERSION "${NANOSEQ_VERSION_RAW}")
  endif()
endif()

project(NanoSeq LANGUAGES C CXX)

# ==============================================================================
# Options
# ==============================================================================
option(NANOSEQ_TESTS "Register NanoSeq's ctest suite" ON)
option(NANOSEQ_INSTALL "Install NanoSeq binaries via 'cmake --install'" OFF)
option(NANOSEQ_FETCH_DEPS
  "Fetch and build htslib/libdeflate ourselves; turn OFF to use pre-installed copies found via CMAKE_PREFIX_PATH instead (e.g. on network-restricted build hosts)"
  ON)

# ==============================================================================
# Build settings
# ==============================================================================
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)
endif()

add_compile_options(-Wall)
add_compile_definitions(NANOSEQ_VERSION="${NANOSEQ_VERSION}")

# Build tree mirrors an install tree's bin/ layout so test/run1-4 (which
# expect $PREFIX/bin/<tool>) can run directly against build output.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

include(GNUInstallDirs)

# ==============================================================================
# Dependencies
# ==============================================================================
find_package(Threads REQUIRED)
find_package(ZLIB REQUIRED)
find_package(CURL REQUIRED)
find_package(OpenSSL REQUIRED)   # provides OpenSSL::Crypto

# versions.sh is the single source of truth for these two version numbers,
# also read directly by build-scripts/opt-build.sh's own htslib/libdeflate
# build (for samtools, which needs its own copy).
file(STRINGS ${CMAKE_SOURCE_DIR}/versions.sh NANOSEQ_VERSIONS_RAW)
foreach(_nanoseq_version_line ${NANOSEQ_VERSIONS_RAW})
  if(_nanoseq_version_line MATCHES "^VER_HTSLIB=(.*)$")
    set(NANOSEQ_HTSLIB_VERSION "${CMAKE_MATCH_1}")
  elseif(_nanoseq_version_line MATCHES "^VER_LIBDEFLATE=(.*)$")
    set(NANOSEQ_LIBDEFLATE_VERSION "${CMAKE_MATCH_1}")
  endif()
endforeach()
# libdeflate's own project() version has no leading "v" (unlike its git tags).
string(REGEX REPLACE "^v" "" NANOSEQ_LIBDEFLATE_VERSION_NUM "${NANOSEQ_LIBDEFLATE_VERSION}")

if(NANOSEQ_FETCH_DEPS)
  # // --- libdeflate --- //
  # Fetched and built as part of this project's own build, statically, with
  # its CLI/test targets disabled - only the libdeflate::libdeflate_static
  # target is needed here. htslib (below) links against this same copy.
  include(FetchContent)
  set(LIBDEFLATE_BUILD_SHARED_LIB OFF CACHE BOOL "" FORCE)
  set(LIBDEFLATE_BUILD_GZIP OFF CACHE BOOL "" FORCE)
  set(LIBDEFLATE_BUILD_TESTS OFF CACHE BOOL "" FORCE)
  FetchContent_Declare(libdeflate
    URL https://github.com/ebiggers/libdeflate/archive/${NANOSEQ_LIBDEFLATE_VERSION}.tar.gz)
  FetchContent_GetProperties(libdeflate)
  if(NOT libdeflate_POPULATED)
    FetchContent_Populate(libdeflate)
    add_subdirectory(${libdeflate_SOURCE_DIR} ${libdeflate_BINARY_DIR} EXCLUDE_FROM_ALL)
  endif()
  set_target_properties(libdeflate_static PROPERTIES POSITION_INDEPENDENT_CODE ON)

  # // --- htslib --- //
  # htslib has no CMake build of its own, so it's fetched and driven through
  # its autotools ./configure via ExternalProject, statically linked against
  # the libdeflate fetched above. Built by the project since nanoseq uses
  # private headers not usually made available by make install.
  include(ExternalProject)
  include(ProcessorCount)
  ProcessorCount(NANOSEQ_NPROC)
  if(NANOSEQ_NPROC EQUAL 0)
    set(NANOSEQ_NPROC 1)
  endif()

  set(HTSLIB_INSTALL_DIR ${CMAKE_BINARY_DIR}/htslib-prefix/install)
  file(MAKE_DIRECTORY ${HTSLIB_INSTALL_DIR}/include)
  ExternalProject_Add(htslib_ext
    URL https://github.com/samtools/htslib/releases/download/${NANOSEQ_HTSLIB_VERSION}/htslib-${NANOSEQ_HTSLIB_VERSION}.tar.bz2
    PREFIX ${CMAKE_BINARY_DIR}/htslib-prefix
    BUILD_IN_SOURCE TRUE
    CONFIGURE_COMMAND <SOURCE_DIR>/configure
      --disable-plugins --enable-libcurl --with-libdeflate
      --prefix=${HTSLIB_INSTALL_DIR}
      CC=${CMAKE_C_COMPILER}
      CFLAGS=-I${libdeflate_SOURCE_DIR}
      LDFLAGS=-L$<TARGET_FILE_DIR:libdeflate::libdeflate_static>
    BUILD_COMMAND make -j${NANOSEQ_NPROC}
    INSTALL_COMMAND make install
      COMMAND ${CMAKE_COMMAND} -E make_directory ${HTSLIB_INSTALL_DIR}/include/cram
      COMMAND bash -c "cp <SOURCE_DIR>/cram/*.h ${HTSLIB_INSTALL_DIR}/include/cram/"
      COMMAND ${CMAKE_COMMAND} -E copy <SOURCE_DIR>/header.h ${HTSLIB_INSTALL_DIR}/include/header.h
    BUILD_BYPRODUCTS ${HTSLIB_INSTALL_DIR}/lib/libhts.a
  )
  add_dependencies(htslib_ext libdeflate_static)

  add_library(htslib::htslib STATIC IMPORTED)
  set_target_properties(htslib::htslib PROPERTIES
    IMPORTED_LOCATION ${HTSLIB_INSTALL_DIR}/lib/libhts.a
    INTERFACE_INCLUDE_DIRECTORIES ${HTSLIB_INSTALL_DIR}/include)
  add_dependencies(htslib::htslib htslib_ext)
else()
  # NANOSEQ_FETCH_DEPS=OFF: look for pre-installed copies instead, e.g. via
  # -DCMAKE_PREFIX_PATH=/path/to/prefix pointing at a prefix built the way
  # build-scripts/opt-build.sh builds one (its --with-libdeflate htslib
  # configure needs the same libdeflate found here).
  include(FindPackageHandleStandardArgs)

  # find preinstalled libdeflate, matching to specified version
  find_package(libdeflate ${NANOSEQ_LIBDEFLATE_VERSION_NUM} EXACT CONFIG REQUIRED)

  # htslib has no CMake package config and no guaranteed pkg-config file, so
  # it's located by hand. Does not ship the internal cram/*.h or header.h
  # by default - a pre-installed copy must make both available, or paths
  # must be passed in via -D.
  find_path(HTSLIB_INCLUDE_DIR
    NAMES htslib/sam.h
    DOC "Directory containing htslib/sam.h")
  find_path(HTSLIB_CRAM_INCLUDE_DIR
    NAMES cram/cram.h
    PATHS ${HTSLIB_INCLUDE_DIR}
    DOC "htslib's cram headers")
  find_path(HTSLIB_HEADER_H_DIR
    NAMES header.h
    PATHS ${HTSLIB_INCLUDE_DIR}
    DOC "htslib's header.h")
  find_library(HTSLIB_LIBRARY
    NAMES libhts.a hts
    DOC "Static htslib archive")
  find_package_handle_standard_args(htslib
    REQUIRED_VARS HTSLIB_LIBRARY HTSLIB_INCLUDE_DIR HTSLIB_CRAM_INCLUDE_DIR HTSLIB_HEADER_H_DIR)
  if(htslib_FOUND AND NOT TARGET htslib::htslib)
    add_library(htslib::htslib STATIC IMPORTED)
    set_target_properties(htslib::htslib PROPERTIES
      IMPORTED_LOCATION "${HTSLIB_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${HTSLIB_INCLUDE_DIR}")
  endif()
  mark_as_advanced(HTSLIB_INCLUDE_DIR HTSLIB_CRAM_INCLUDE_DIR HTSLIB_HEADER_H_DIR HTSLIB_LIBRARY)
endif()

# // --- system deps via pkg-config --- //
find_library(LZMA_LIBRARY NAMES lzma REQUIRED)
find_library(BZ2_LIBRARY NAMES bz2 REQUIRED)
# dl/m are not separate findable libraries on all platforms (e.g. macOS's
# libSystem) - hence not REQUIRED
find_library(DL_LIBRARY NAMES dl)
find_library(M_LIBRARY NAMES m)

# // --- vendored gztream --- //
add_library(gzstream STATIC ${CMAKE_SOURCE_DIR}/vendor/gzstream-1.5/gzstream.C)
target_include_directories(gzstream PUBLIC ${CMAKE_SOURCE_DIR}/vendor/gzstream-1.5)
target_link_libraries(gzstream PUBLIC ZLIB::ZLIB)

# ==============================================================================
# Targets
# ==============================================================================
# Shared link requirements. Order matters for static
# archives: htslib before libdeflate.
add_library(nanoseq_common INTERFACE)
target_link_libraries(nanoseq_common INTERFACE
  htslib::htslib
)
if(TARGET libdeflate::libdeflate_static)
  target_link_libraries(nanoseq_common INTERFACE libdeflate::libdeflate_static)
else()
  target_link_libraries(nanoseq_common INTERFACE libdeflate::libdeflate_shared)
endif()
target_link_libraries(nanoseq_common INTERFACE
  gzstream
  ZLIB::ZLIB
  Threads::Threads
  CURL::libcurl
  OpenSSL::Crypto
  ${LZMA_LIBRARY}
  ${BZ2_LIBRARY}
)
if(DL_LIBRARY)
  target_link_libraries(nanoseq_common INTERFACE ${DL_LIBRARY})
endif()
if(M_LIBRARY)
  target_link_libraries(nanoseq_common INTERFACE ${M_LIBRARY})
endif()

add_executable(bamaddreadbundles src/bamaddreadbundles.cc)
target_link_libraries(bamaddreadbundles PRIVATE nanoseq_common)

add_executable(gridsearch src/gridsearch.cc)
target_link_libraries(gridsearch PRIVATE nanoseq_common)

add_executable(dsa src/dsa.cc src/pileup.cc src/bed_reader.cc src/read_bundler.cc src/writeout.cc)
target_link_libraries(dsa PRIVATE nanoseq_common)

add_executable(statistics src/statistics.cc)
target_link_libraries(statistics PRIVATE nanoseq_common)

add_executable(variantcaller src/variantcaller.cc)
target_link_libraries(variantcaller PRIVATE nanoseq_common)

add_executable(randomreadinbundle src/randomreadinbundle.cc)
target_link_libraries(randomreadinbundle PRIVATE nanoseq_common)

# C executable
add_executable(bamcov src/bamcov.c src/sam_opts.c src/sam_utils.c)
set_target_properties(bamcov PROPERTIES LINKER_LANGUAGE C)
target_link_libraries(bamcov PRIVATE nanoseq_common)

# ==============================================================================
# Install
# ==============================================================================
if(NANOSEQ_INSTALL)
  install(TARGETS
    bamaddreadbundles
    gridsearch
    dsa
    statistics
    variantcaller
    randomreadinbundle
    bamcov
    RUNTIME DESTINATION bin)
endif()

# ==============================================================================
# Tests
# ==============================================================================
if(NANOSEQ_TESTS)
  enable_testing()

  # run1-4 resolve binaries via $PREFIX/bin/<tool> and expect cwd == test/.
  # CMAKE_RUNTIME_OUTPUT_DIRECTORY (set above) already puts binaries under
  # ${CMAKE_BINARY_DIR}/bin, so pointing PREFIX there lets these run directly
  # against build output with no install step required first.
  foreach(script run1 run2 run3 run4)
    add_test(NAME ${script}
             COMMAND ${CMAKE_SOURCE_DIR}/test/${script}
             WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test)
    set_tests_properties(${script} PROPERTIES
      ENVIRONMENT "PREFIX=${CMAKE_BINARY_DIR}")
  endforeach()

  # smoke tests for binaries not tested above
  foreach(tool gridsearch bamcov statistics)
    add_test(NAME ${tool}_smoke COMMAND $<TARGET_FILE:${tool}> -h)
  endforeach()
endif()
