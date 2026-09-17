cmake_minimum_required(VERSION 3.18)
project(DISCON_MPC_CPP LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(MPC_ENABLE_TIMING "Enable low-overhead per-stage MPC timing instrumentation" OFF)

add_library(DISCON_MPC_CPP SHARED DISCON_MPC_CPP.cpp)

target_compile_definitions(DISCON_MPC_CPP PRIVATE _CRT_SECURE_NO_WARNINGS)
if(MPC_ENABLE_TIMING)
    target_compile_definitions(DISCON_MPC_CPP PRIVATE MPC_ENABLE_TIMING=1)
endif()

target_include_directories(DISCON_MPC_CPP PRIVATE "C:/Dev/qpOASES/include")
target_link_libraries(DISCON_MPC_CPP PRIVATE "C:/Dev/qpOASES-build/libs/Release/qpOASES.lib")

set_target_properties(DISCON_MPC_CPP PROPERTIES
    OUTPUT_NAME "DISCON"
    WINDOWS_EXPORT_ALL_SYMBOLS OFF
)

# A cleaner long-term approach is to add qpOASES as a subdirectory or create
# an imported target. The current setup intentionally hard-wires the local
# library path to validate the minimal MPC/DLL integration as quickly as possible.
