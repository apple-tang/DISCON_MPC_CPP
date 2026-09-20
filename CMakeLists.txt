cmake_minimum_required(VERSION 3.18)
project(DISCON_MPC_CPP LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
option(MPC_ENABLE_TIMING "Record per-stage MPC timings" OFF)
option(MPC_BUILD_TESTS "Build controller unit and interface smoke tests" OFF)
set(QPOASES_ROOT "C:/Dev/qpOASES" CACHE PATH "qpOASES installation or source root")
find_path(QPOASES_INCLUDE_DIR qpOASES.hpp HINTS "${QPOASES_ROOT}/include" REQUIRED)
find_library(QPOASES_LIBRARY NAMES qpOASES
    HINTS "${QPOASES_ROOT}/lib" "${QPOASES_ROOT}/libs"
          "C:/Dev/qpOASES-build/libs/Release" REQUIRED)

function(configure_mpc_target target)
    target_compile_definitions(${target} PRIVATE _CRT_SECURE_NO_WARNINGS)
    target_include_directories(${target} PRIVATE "${QPOASES_INCLUDE_DIR}")
    target_link_libraries(${target} PRIVATE "${QPOASES_LIBRARY}")
    if(MPC_ENABLE_TIMING)
        target_compile_definitions(${target} PRIVATE MPC_ENABLE_TIMING=1)
    endif()
endfunction()

add_library(DISCON_MPC_CPP SHARED DISCON_MPC_CPP.cpp)
configure_mpc_target(DISCON_MPC_CPP)
set_target_properties(DISCON_MPC_CPP PROPERTIES OUTPUT_NAME "DISCON"
    WINDOWS_EXPORT_ALL_SYMBOLS OFF)

if(MPC_BUILD_TESTS)
    enable_testing()
    add_executable(controller_tests tests/controller_tests.cpp)
    configure_mpc_target(controller_tests)
    add_test(NAME controller_tests COMMAND controller_tests "${CMAKE_CURRENT_SOURCE_DIR}")
endif()
