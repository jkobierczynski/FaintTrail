cmake_minimum_required(VERSION 3.20)
project(fainttrail VERSION 0.5.0 LANGUAGES C CXX
        DESCRIPTION "Synthetic tracking engine for faint asteroid and comet detection")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(FAINTTRAIL_BUILD_CLI   "Build the fainttrail command-line tool" ON)
option(FAINTTRAIL_BUILD_GUI   "Build the Qt 6 GUI (fainttrail-gui)"     ON)
option(FAINTTRAIL_BUILD_TESTS "Build unit tests"                         ON)
option(FAINTTRAIL_BUILD_TOOLS "Build the validation benchmark"        OFF)
option(FAINTTRAIL_NATIVE      "Optimise for the build machine (-march=native)" OFF)

find_package(OpenMP)

# ---------------------------------------------------------------- library
add_library(fainttrail
  src/fits.cpp
  src/wcs.cpp
  src/simulate.cpp
  src/preprocess.cpp
  src/search.cpp
  src/validate.cpp
  src/io.cpp
  src/astrometry.cpp
  src/regression.cpp
  src/backend.cpp
  src/backend_cpu.cpp
  src/c_api.cpp)
add_library(fainttrail::fainttrail ALIAS fainttrail)

target_include_directories(fainttrail PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)
target_compile_definitions(fainttrail PRIVATE
  FAINTTRAIL_VERSION_STRING="${PROJECT_VERSION}")

if(MSVC)
  target_compile_options(fainttrail PRIVATE /W4 /O2)
else()
  target_compile_options(fainttrail PRIVATE -Wall -Wextra -Wpedantic
    $<$<CONFIG:Release>:-O3>)
  if(FAINTTRAIL_NATIVE)
    target_compile_options(fainttrail PRIVATE -march=native)
  endif()
endif()

if(OpenMP_CXX_FOUND)
  target_link_libraries(fainttrail PRIVATE OpenMP::OpenMP_CXX)
  target_compile_definitions(fainttrail PRIVATE FAINTTRAIL_HAVE_OPENMP=1)
endif()

# ---------------------------------------------------------------- CLI
if(FAINTTRAIL_BUILD_CLI)
  add_executable(fainttrail_cli cli/main.cpp)
  set_target_properties(fainttrail_cli PROPERTIES OUTPUT_NAME fainttrail)
  target_link_libraries(fainttrail_cli PRIVATE fainttrail)
endif()

# ---------------------------------------------------------------- GUI
if(FAINTTRAIL_BUILD_GUI)
  find_package(Qt6 6.2 COMPONENTS Widgets Concurrent QUIET)
  if(Qt6_FOUND)
    set(CMAKE_AUTOMOC ON)
    add_executable(fainttrail-gui
      gui/main.cpp
      gui/mainwindow.cpp
      gui/mainwindow.hpp
      gui/imageview.cpp
      gui/imageview.hpp)
    target_link_libraries(fainttrail-gui PRIVATE fainttrail Qt6::Widgets Qt6::Concurrent)
    message(STATUS "fainttrail: building Qt ${Qt6_VERSION} GUI")
  else()
    message(STATUS "fainttrail: Qt 6 not found, GUI disabled")
  endif()
endif()

# ---------------------------------------------------------------- tests
if(FAINTTRAIL_BUILD_TESTS)
  enable_testing()
  add_executable(fainttrail_tests tests/test_main.cpp tests/test_c_api.c)
  target_link_libraries(fainttrail_tests PRIVATE fainttrail)
  target_include_directories(fainttrail_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
  add_test(NAME fainttrail_tests COMMAND fainttrail_tests)
  set_tests_properties(fainttrail_tests PROPERTIES TIMEOUT 600)
  # Real-data regression cases: small fixtures in the repository run everywhere;
  # full-frame cases are skipped (exit code 77) unless FAINTTRAIL_REGRESSION_DATA
  # points at the full data sets.
  if(FAINTTRAIL_BUILD_CLI)
    file(GLOB FAINTTRAIL_REGRESSION_CASES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/*.case)
    foreach(case_file ${FAINTTRAIL_REGRESSION_CASES})
      get_filename_component(case_name ${case_file} NAME_WE)
      add_test(NAME regress_${case_name} COMMAND fainttrail_cli regress ${case_file} --strict --quiet)
      set_tests_properties(regress_${case_name} PROPERTIES SKIP_RETURN_CODE 77 TIMEOUT 3600 LABELS regression)
    endforeach()
  endif()
endif()

# ---------------------------------------------------------------- tools
if(FAINTTRAIL_BUILD_TOOLS)
  add_executable(fainttrail_benchmark tools/benchmark.cpp)
  target_link_libraries(fainttrail_benchmark PRIVATE fainttrail)
  add_executable(fainttrail_kernel_bench tools/kernel_bench.cpp)
  target_link_libraries(fainttrail_kernel_bench PRIVATE fainttrail)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(fainttrail_kernel_bench PRIVATE OpenMP::OpenMP_CXX)
  endif()
endif()

# ---------------------------------------------------------------- install
include(GNUInstallDirs)
install(TARGETS fainttrail
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(DIRECTORY include/fainttrail DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
if(TARGET fainttrail_cli)
  install(TARGETS fainttrail_cli RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
if(TARGET fainttrail-gui)
  install(TARGETS fainttrail-gui RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
install(FILES LICENSE README.md DESTINATION ${CMAKE_INSTALL_DOCDIR})
