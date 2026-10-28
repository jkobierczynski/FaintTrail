cmake_minimum_required(VERSION 3.20)
project(fainttrail VERSION 0.9.0 LANGUAGES C CXX
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
option(FAINTTRAIL_CUDA        "Build the CUDA stacking backend (needs the CUDA Toolkit)" OFF)

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
  src/memory.cpp
  src/tiled.cpp
  src/backend.cpp
  src/backend_cpu.cpp
  src/cuda_plan.cpp
  src/registration.cpp
  src/solve.cpp
  src/process.cpp
  src/photometry.cpp
  src/c_api.cpp)
add_library(fainttrail::fainttrail ALIAS fainttrail)

target_include_directories(fainttrail PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)
target_compile_definitions(fainttrail PRIVATE
  FAINTTRAIL_VERSION_STRING="${PROJECT_VERSION}")

if(MSVC)
  target_compile_options(fainttrail PRIVATE "$<$<COMPILE_LANGUAGE:CXX>:/W4>")
else()
  target_compile_options(fainttrail PRIVATE "$<$<COMPILE_LANGUAGE:CXX>:-Wall;-Wextra;-Wpedantic>"
    $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:CXX>>:-O3>)
  if(FAINTTRAIL_NATIVE)
    target_compile_options(fainttrail PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-march=native>)
  endif()
endif()

# ---------------------------------------------------------------- CUDA backend
# Blackwell cards (RTX 50xx, sm_120) need CUDA 12.8 or newer. By default the code
# is built for the GPU in the build machine ("native"); pass
# -DCMAKE_CUDA_ARCHITECTURES="86;89;120" (for example) to build for others.
if(FAINTTRAIL_CUDA)
  if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.24)
      set(CMAKE_CUDA_ARCHITECTURES native)
    else()
      set(CMAKE_CUDA_ARCHITECTURES 75 86 89)
    endif()
  endif()
  enable_language(CUDA)
  find_package(CUDAToolkit REQUIRED)
  set(CMAKE_CUDA_STANDARD 17)
  set(CMAKE_CUDA_STANDARD_REQUIRED ON)
  target_sources(fainttrail PRIVATE src/backend_cuda.cu)
  target_compile_definitions(fainttrail PRIVATE FAINTTRAIL_HAVE_CUDA=1)
  # no fused multiply-adds, so the results stay bit-identical to the CPU
  if(CMAKE_CUDA_COMPILER_ID STREQUAL "NVIDIA")
    target_compile_options(fainttrail PRIVATE "$<$<COMPILE_LANGUAGE:CUDA>:--fmad=false;-lineinfo>")
  else()  # Clang as CUDA compiler
    target_compile_options(fainttrail PRIVATE "$<$<COMPILE_LANGUAGE:CUDA>:-ffp-contract=off>")
  endif()
  target_link_libraries(fainttrail PRIVATE CUDA::cudart_static)
  if(MSVC)
    # cudart_static.lib names the static C runtime (LIBCMT) as its default library;
    # the project uses the DLL runtime, so keep LIBCMT out of the link (LNK4098)
    target_link_options(fainttrail INTERFACE /NODEFAULTLIB:LIBCMT)
  endif()
  message(STATUS "FaintTrail: CUDA backend enabled (CUDA ${CUDAToolkit_VERSION}, architectures ${CMAKE_CUDA_ARCHITECTURES})")
endif()

if(WIN32)
  target_link_libraries(fainttrail PRIVATE psapi)
endif()
if(OpenMP_CXX_FOUND)
  target_link_libraries(fainttrail PRIVATE OpenMP::OpenMP_CXX)
  target_compile_definitions(fainttrail PRIVATE FAINTTRAIL_HAVE_OPENMP=1)
endif()

# ---------------------------------------------------------------- CLI
if(FAINTTRAIL_BUILD_CLI)
  add_executable(fainttrail_cli cli/main.cpp)
  # On Windows an executable that exports symbols (cudart_static exports
  # NvOptimusEnablementCuda) gets an import library; give it a name that does not
  # collide with the static library fainttrail.lib (LNK1149)
  set_target_properties(fainttrail_cli PROPERTIES OUTPUT_NAME fainttrail ARCHIVE_OUTPUT_NAME fainttrail_cli)
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
    # a windowed application on Windows (no console window); Qt supplies WinMain
    set_target_properties(fainttrail-gui PROPERTIES WIN32_EXECUTABLE ON)
    if(APPLE)  # an application bundle: fainttrail-gui.app
      set_target_properties(fainttrail-gui PROPERTIES
        MACOSX_BUNDLE ON
        MACOSX_BUNDLE_BUNDLE_NAME "FaintTrail"
        MACOSX_BUNDLE_GUI_IDENTIFIER "io.github.jkobierczynski.fainttrail"
        MACOSX_BUNDLE_BUNDLE_VERSION "${PROJECT_VERSION}"
        MACOSX_BUNDLE_SHORT_VERSION_STRING "${PROJECT_VERSION}")
    endif()
    if(WIN32)
      # "cmake --build build --config Release --target deploy_gui" copies the Qt
      # DLLs and plugins next to fainttrail-gui.exe, so the folder runs anywhere
      if(TARGET Qt6::windeployqt)
        set(FAINTTRAIL_WINDEPLOYQT "$<TARGET_FILE:Qt6::windeployqt>")
      else()
        find_program(FAINTTRAIL_WINDEPLOYQT windeployqt
          HINTS "${Qt6_DIR}/../../../bin" "${QT6_INSTALL_PREFIX}/bin")
      endif()
      if(FAINTTRAIL_WINDEPLOYQT)
        add_custom_target(deploy_gui
          COMMAND "${FAINTTRAIL_WINDEPLOYQT}" --release --no-translations --compiler-runtime "$<TARGET_FILE:fainttrail-gui>"
          DEPENDS fainttrail-gui
          COMMENT "Copying the Qt runtime next to fainttrail-gui.exe")
        message(STATUS "fainttrail: deploy_gui target uses ${FAINTTRAIL_WINDEPLOYQT}")
      else()
        message(STATUS "fainttrail: windeployqt not found, no deploy_gui target")
      endif()
    endif()
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
  if(FAINTTRAIL_CUDA AND FAINTTRAIL_BUILD_CLI)
    # skipped (77) on machines without a usable GPU; with one, the regression
    # cases below also run on it, since "auto" picks the CUDA backend
    add_test(NAME check_backend_cuda COMMAND fainttrail_cli check-backend --backend cuda)
    set_tests_properties(check_backend_cuda PROPERTIES SKIP_RETURN_CODE 77 TIMEOUT 600)
  endif()
  # Real-data regression cases: small fixtures in the repository run everywhere;
  # full-frame cases are skipped (exit code 77) unless FAINTTRAIL_REGRESSION_DATA
  # points at the full data sets.
  if(FAINTTRAIL_BUILD_CLI)
    file(GLOB FAINTTRAIL_REGRESSION_CASES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/*.case)
    foreach(case_file ${FAINTTRAIL_REGRESSION_CASES})
      get_filename_component(case_name ${case_file} NAME_WE)
      add_test(NAME regress_${case_name} COMMAND fainttrail_cli regress ${case_file} --strict --quiet)
      set_tests_properties(regress_${case_name} PROPERTIES SKIP_RETURN_CODE 77 TIMEOUT 3600 LABELS regression)
      # the small fixtures again as tiled searches: they must give the same result
      if(case_name MATCHES "-cut$")
        add_test(NAME regress_${case_name}_tiled COMMAND fainttrail_cli regress ${case_file} --strict --quiet --tile 48)
        set_tests_properties(regress_${case_name}_tiled PROPERTIES TIMEOUT 600 LABELS regression)
      endif()
    endforeach()
  endif()
endif()

# ---------------------------------------------------------------- tools
if(FAINTTRAIL_BUILD_TOOLS)
  add_executable(fainttrail_benchmark tools/benchmark.cpp)
  target_link_libraries(fainttrail_benchmark PRIVATE fainttrail)
  add_executable(fainttrail_false_alarm_study tools/false_alarm_study.cpp)
  target_link_libraries(fainttrail_false_alarm_study PRIVATE fainttrail)
  add_executable(fainttrail_null_flags_study tools/null_flags_study.cpp)
  target_link_libraries(fainttrail_null_flags_study PRIVATE fainttrail)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(fainttrail_null_flags_study PRIVATE OpenMP::OpenMP_CXX)
  endif()
  add_executable(fainttrail_kernel_feature_study tools/kernel_feature_study.cpp)
  target_link_libraries(fainttrail_kernel_feature_study PRIVATE fainttrail)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(fainttrail_kernel_feature_study PRIVATE OpenMP::OpenMP_CXX)
  endif()
  add_executable(fainttrail_calibration_study tools/calibration_study.cpp)
  target_link_libraries(fainttrail_calibration_study PRIVATE fainttrail)
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
install(FILES LICENSE README.md DEVELOPMENT.md BUILDING.md DESTINATION ${CMAKE_INSTALL_DOCDIR})
