# DISCON_MPC_CPP

Minimal C++ DLL shell for OpenFAST ServoDyn / Bladed-style `DISCON` control.

## Purpose

This directory provides a starting point for migrating the current
Fortran fracture controller toward a C++ implementation that can call
`qpOASES` for online MPC.

## Current behavior

- exports a `DISCON` symbol compatible with the current Bladed DLL interface
- reads the main `avrSWAP` inputs already used by the Fortran controller
- reads the following real-time signals from `avrSWAP`:
  - rotor speed
  - horizontal hub-height wind speed
  - tower-top fore-aft velocity
  - tower-top fore-aft displacement
- treats the current wind-speed deviation as a constant preview disturbance over the prediction horizon
- solves a minimal multi-step MPC/QP with `qpOASES`
- writes back generator torque, collective pitch, and collective pitch rate

This means it is now a **minimal working MPC prototype**, but it still uses
placeholder identified sensitivities and a short fixed prediction horizon.

## Files

- `DISCON_MPC_CPP.cpp`
  - minimal DLL entry and short-horizon MPC prototype logic
- `CMakeLists.txt`
  - standalone build file for the C++ DLL

## Next steps to turn this into MPC

1. Replace placeholder physical parameters with identified values
2. Replace constant-over-horizon wind preview with a richer preview source if desired
3. Refine identified physical parameters and weights
4. Add richer state/output constraints if desired
5. Externalize tuning parameters to an input file

## Externalized MPC tuning parameters

The current `.IN` file now externalizes:

- `FractureTime`
- `Q_OMEGA`
- `Q_X`
- `Q_V`
- `R_T`
- `R_B`
- `OMEGA_ERR_MAX`
- `TOWER_DISP_MAX`
- `TOWER_VEL_MAX`

The prediction horizon lengths and model dimensions remain compile-time constants.

## Example build (Visual Studio toolchain)

```powershell
cmd /c '\"C:\\Program Files\\Microsoft Visual Studio\\2022\\Community\\VC\\Auxiliary\\Build\\vcvars64.bat\" && ^
\"C:\\Program Files\\Microsoft Visual Studio\\2022\\Community\\Common7\\IDE\\CommonExtensions\\Microsoft\\CMake\\CMake\\bin\\cmake.exe\" -S . -B build -A x64 && ^
\"C:\\Program Files\\Microsoft Visual Studio\\2022\\Community\\Common7\\IDE\\CommonExtensions\\Microsoft\\CMake\\CMake\\bin\\cmake.exe\" --build build --config Release'
```
