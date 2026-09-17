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
- `N_PRED`
- `N_CTRL_H`
- `MPC_DT`
- `CTRL_DT`
- `QP_HOTSTART`
- `Q_OMEGA`
- `Q_X`
- `Q_V`
- `Q_TG`
- `Q_BETA`
- `R_T`
- `R_B`
- `OMEGA_ERR_MAX`
- `TOWER_DISP_MAX`
- `TOWER_VEL_MAX`
- `nWSR`

`MPC_DT` is the internal prediction-model discretization step used by the MPC horizon. It does not
need to match the incoming OpenFAST call interval. The applied generator-torque and pitch moves are
still re-limited against the physical actuator-rate bounds using the actual controller call interval
before they are written back to `avrSWAP`, so increasing `MPC_DT` for prediction does not remove
the final output protection against torque/pitch overshoot.

`CTRL_DT` is the online solve interval. With a 1 ms OpenFAST step and `CTRL_DT=0.010`, the MPC solves
on the first post-fracture call and every 10 ms thereafter. The previous torque/pitch target is held
on the nine intermediate calls, while any enabled first-order actuator dynamics still update every
1 ms. The supplied configuration uses `MPC_DT=CTRL_DT=0.010 s`.

The held MPC target is not written directly to OpenFAST. At every positive-time controller call, the
applied generator torque and collective pitch move toward that target subject to the hard physical
rates `VS_MAX_TQ_RATE` and `PC_MAX_RAT`. With a 1 ms OpenFAST step these limits are 15 N-m and
0.008 degrees per step. Repeated calls at the same simulation time do not advance actuator state.

`QP_HOTSTART=0` constructs and initializes a new `QProblem` on every solve. `QP_HOTSTART=1` uses an
`SQProblem`: the first solve is cold, later solves update the changing matrices and reuse the previous
working set. A failed hotstart is retried once from a cold start. Timing CSV rows identify mode 0
(cold), 1 (hot), or 2 (hotstart failed and cold retry used), together with the consumed `nWSR` count.

`nWSR` is the qpOASES working-set iteration budget used for each online QP solve. Increase it when
larger `N_CTRL_H` or added state constraints begin to hit `RET_MAX_NWSR_REACHED`, at the cost of
more online solver work per controller step.

## Measuring U-computation time

Timing instrumentation is compile-time optional, so the ordinary DLL has no timing calls. Configure
with `-DMPC_ENABLE_TIMING=ON` to enable it. The instrumented DLL writes
`<OpenFAST-output-root>.mpc_u_timing.csv`; gain tables are still read from the six CSV files only once
during controller initialization.

Each successful MPC solve records monotonic start/end timestamps and these durations in microseconds:

- preprocessing and gain-table interpolation
- local-model construction and workspace reset
- condensed prediction matrices `G` and `c`
- QP Hessian and gradient
- input/rate/state constraints
- qpOASES solve
- primal-solution extraction and first-command limiting
- post-solve diagnostics
- total time until `U` is ready and total function time

Rows are buffered and written in batches of 256, then flushed on the final OpenFAST controller call.
This keeps disk I/O outside the measured `total_to_u_us` interval. Use `benchmark_u_timing.py` for a
repeatable DLL-only benchmark and `analyze_u_timing.py` to produce percentile summaries and a plot.
At controller initialization, an existing timing CSV with the same output root is deleted and recreated.
If another application has locked that file, timing output is disabled for the run instead of appending
new rows to stale data; the controller startup message reports `timingCsv=disabled-file-busy`.

```powershell
python benchmark_u_timing.py --dll timing_release\DISCON_timing.dll --dt 0.001 --samples 100 --warmup 20
python analyze_u_timing.py benchmark_u_timing_samples.csv --control-period-us 10000
```

## Rebuilding 5 MW aerodynamic gain tables

The legacy six tables were built from 667 separate files covering 29 wind speeds and 23 rotor speeds;
they were not copied from one 24 m/s case. However, the legacy builder used full-model DC gains and
mixed kN/rpm table units with the SI units required by the controller equations.

The replacement workflow uses local OpenFAST C/D Jacobian entries for the three controller-model
coordinates `Omega [rad/s]`, collective pitch `[rad]`, and wind speed `[m/s]`. It never evaluates
`-C A^-1 B + D`, and it refuses to create a rectangular table when a requested point is missing.

```powershell
# 1. Select post-fracture rotor-speed crossings from one time-domain output per wind speed.
python find_5mw_rotspeed_times.py --cases-root <time-domain-root> --require-complete

# 2. Create isolated Linearize=True cases and a traceable manifest.
python generate_5mw_linearization_cases.py --schedule 5mw_linearization_schedule.csv `
  --source-baseline <5MW_Baseline> --output-root openfast_5mw_linearization_cases

# 3. After manually running those cases, extract the six SI-unit tables.
python extract_5mw_mpc_gain_tables.py --lin-root openfast_5mw_linearization_cases `
  --manifest openfast_5mw_linearization_cases\linearization_manifest.csv `
  --config-template DISCON_MPC_CPP.IN --output-dir generated_5mw_gain_tables
```

For the requested 12, 14, 16, 18, 20, 22, and 24 m/s comparison, a pre-extracted candidate set is
in `generated_5mw_gain_tables_12to24_even`. Its matching `DISCON_MPC_CPP.IN` points to the tables in
that directory. Keep the legacy and local-Jacobian results separate until closed-loop validation is
complete.

The main six CSV files were updated on 2026-09-15 from the complete 29-by-23 local-Jacobian grid.
The incorrectly labelled legacy 21 m/s row was regenerated from a fresh time-domain trajectory and
23 new OpenFAST linearizations before extraction. The replaced DC-gain tables and their matching
configuration are preserved in `legacy_gain_tables_before_si_repair_20260915`.

## Example build (Visual Studio toolchain)

```powershell
cmd /c '\"C:\\Program Files\\Microsoft Visual Studio\\2022\\Community\\VC\\Auxiliary\\Build\\vcvars64.bat\" && ^
\"C:\\Program Files\\Microsoft Visual Studio\\2022\\Community\\Common7\\IDE\\CommonExtensions\\Microsoft\\CMake\\CMake\\bin\\cmake.exe\" -S . -B build -A x64 && ^
\"C:\\Program Files\\Microsoft Visual Studio\\2022\\Community\\Common7\\IDE\\CommonExtensions\\Microsoft\\CMake\\CMake\\bin\\cmake.exe\" --build build --config Release'
```
