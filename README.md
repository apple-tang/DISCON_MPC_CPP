# DISCON_MPC_CPP

Model predictive emergency shutdown control following partial-span wind-turbine blade breakage. The controller coordinates generator torque and collective pitch using a constrained quadratic program solved by qpOASES, through ServoDyn's Bladed-style `DISCON` interface.

**Author:** Daozhi Tang, Zhejiang University

**Contact:** [daozhitang@zju.edu.cn](mailto:daozhitang@zju.edu.cn)

**Research:** *An Emergency Shutdown Framework for Wind Turbines Following Blade Breakage*

The contribution is the blade-breakage shutdown scenario and the application of MPC to coordinated deceleration and tower-vibration mitigation. This repository contains the controller, not the modified OpenFAST plant, a fracture detector, or a learned world model. World-model integration is a future research direction.

## Control behavior

1. Before `FRACTURE_TIME`, run the NREL 5-MW baseline torque and collective-pitch logic.
2. After the time trigger, latch MPC operation. Wind speed and rotor speed schedule local aerodynamic derivatives. Hold the current wind deviation over the prediction horizon.
3. Optimize absolute torque/pitch commands, penalizing state deviations and successive command changes. Apply the first command and hold the target between solves.
4. Enforce magnitude and command-rate constraints, plus a predicted rotor-speed lower bound. Independently rate-limit applied outputs at every positive-time OpenFAST callback.
5. If warm start fails, retry with a fresh QP. If the solve still fails, use the scheduled shutdown fallback and report a warning. The same final actuator limits apply.

The five-state model contains rotor-speed deviation, tower-top fore-aft displacement and velocity, and applied torque/pitch deviations. The current tuning preserves the previous default: terminal cost disabled (`gTerminalCostScale = 0`), with a scheduled LQR gain for fallback. The final prediction block has zero state weight. There is no implemented terminal invariant-set constraint or general stability guarantee. Near standstill, the rotor-speed constraint is relaxed using the existing wind-band thresholds and a 0.2 rpm floor.

The Riccati calculation retains its 500-iteration budget. A nonconverged grid point receives zero feedback gain, so fallback at that point is rate-limited torque release and feathering. The startup message reports `terminalFallbackPts`; the supplied tuning currently reaches this fallback at all 667 grid points. The controller must therefore not be described as having a verified stabilizing LQR backup under this tuning.

| Actuator | Magnitude | Maximum rate |
| --- | --- | --- |
| Generator torque | 0 to 47,402.91 N m | 15,000 N m/s |
| Collective pitch | 0 to 90 degrees | 8 degrees/s |

The actuator response is hard rate-limited. The unused optional first-order lag has been removed. Repeated calls at the same time do not advance applied outputs. Prediction step `MPC_DT` and solve interval `CTRL_DT` are independent of the elapsed physical time used for final output rate limiting.

## OpenFAST interface

Use the matching **modified OpenFAST build**. Standard OpenFAST does not provide this project's custom tower/fracture records. The caller must allocate and populate `avrSWAP` through record 1027. Numbers below are one-based; C++ indices are one less.

| Record | Input | Unit |
| --- | --- | --- |
| 1, 2 | Controller status, simulation time | -, s |
| 4 | Blade 1 pitch | rad |
| 20, 21 | Generator speed, rotor speed | rad/s |
| 27 | Horizontal hub-height wind speed | m/s |
| 1019, 1021 | Tower-top fore-aft velocity, displacement | m/s, m |
| 1023 | Post-fracture rotor inertia | kg m2 |
| 1024 | Blade mass-imbalance first moment | kg m |
| 1025, 1026 | Normalized fracture location, residual property ratio | - |
| 1027 | Fracture status: 0 inactive, 1 transitioning, 2 finalized | - |

The controller locks a valid finalized fracture inertia and rebuilds the fallback schedule. Until then, it uses rotor inertia 2.56488355e7 kg m2 plus generator inertia referred through the 97:1 gearbox. `FRACTURE_TIME` switches control; the plant's structural fracture must be configured separately at the matching time.

Outputs include torque at record 47, pitch at records 42-45, pitch rate at record 46, and contactor/override flags. Only one turbine instance per loaded DLL is supported because controller state is process-global.

## Build

Requirements: C++17, CMake 3.18+, and a matching-architecture qpOASES library. On Windows, use an x64 Native Tools Command Prompt for Visual Studio. Match the compiler runtime to the qpOASES build.

```powershell
cmake -S . -B build -A x64 `
  -DQPOASES_INCLUDE_DIR=C:/Dev/qpOASES/include `
  -DQPOASES_LIBRARY=C:/Dev/qpOASES-build/libs/Release/qpOASES.lib
cmake --build build --config Release
```

Output: `build/Release/DISCON.dll`. Override the two dependency paths or set `QPOASES_ROOT` for another installation. Set ServoDyn's `DLL_FileName` to the compiled DLL, `DLL_InFile` to `DISCON_MPC_CPP.IN`, and `DLL_ProcName` to `DISCON`. Use the DLL torque and pitch control modes in the existing turbine deck.

## Configuration

The single IN file uses `KEY = VALUE`, in any order. `#` and `!` begin comments. All 13 keys are required. Unknown/duplicate keys, invalid numbers, unsorted grids and incorrectly sized tables fail initialization. **Old positional IN files are unsupported: replace the IN file together with the new DLL.**

| Key | Supplied value | Meaning |
| --- | --- | --- |
| `FRACTURE_TIME` | `30.0` | Controller switch time, s |
| `MPC_DT`, `CTRL_DT` | `0.010`, `0.010` | Prediction step and solve interval, s |
| `N_PRED`, `N_CTRL` | `330`, `70` | Prediction steps and independent command blocks |
| `Q` | `30,200,200,1e-6,400` | State weights in the five-state order above |
| `R` | `5,60` | Torque- and pitch-command change weights |
| `NWSR` | `500` | qpOASES iteration budget per attempt |
| `QP_START` | `cold` | `cold` or `warm` |
| `DEBUG` | `0` | `1` enables debug and command traces |
| `TABLE_DIR` | `.` | CSV directory relative to the IN file |
| `SPEED_POINTS_RPM` | 1 to 12, step 0.5 | Ascending rotor-speed grid |
| `WIND_POINTS_MS` | 2 to 30, step 1 | Ascending wind-speed grid |

The supplied prediction horizon is 3.3 s; the independent control horizon is 0.7 s. These values preserve the local controller configuration, not necessarily every setting used in the paper's reported experiments.

### Cold start, warm start and debug

Edit these two keys in the same IN file:

| Run | `QP_START` | `DEBUG` |
| --- | --- | --- |
| Cold start on every solve | `cold` | `0` |
| Warm start after the first cold solve | `warm` | `0` |
| Debug either solver mode | `cold` or `warm` | `1` |

Cold start constructs a new `QProblem`. Warm start uses `SQProblem` to update changing matrices and reuse the previous working set. Debug is independent of solver initialization.

Debug writes `<output-root>.mpc_debug.log` and `<output-root>.mpc_command_trace.csv`. The command trace distinguishes demanded targets from applied outputs and records states, fallback flags and one-step prediction diagnostics. It can support actuator inspection and future model-learning data preparation.

## Gain tables

The six headerless CSV files have 29 wind rows and 23 speed columns, ordered by the IN-file grids. They are loaded once at initialization and interpolated with boundary clamping.

| File | Derivative | Unit |
| --- | --- | --- |
| `GTomega.csv` | Aerodynamic torque / rotor speed | N m / (rad/s) |
| `GTbeta.csv` | Aerodynamic torque / collective pitch | N m / rad |
| `GTu.csv` | Aerodynamic torque / wind speed | N m / (m/s) |
| `GFomega.csv` | Rotor thrust / rotor speed | N / (rad/s) |
| `GFbeta.csv` | Rotor thrust / collective pitch | N / rad |
| `GFu.csv` | Rotor thrust / wind speed | N / (m/s) |

The filename suffix `u` denotes wind speed. These are local OpenFAST C/D Jacobian sensitivities, not stable DC gains from `-C A^-1 B + D`. The active tables were repaired on 2026-09-15, including the regenerated 21 m/s row. Cleanup preserves the six active CSV files byte for byte; grid provenance is retained under `docs/`.

## Timing and tests

Optional timing instrumentation is independent of `DEBUG`:

```powershell
cmake -S . -B build -DMPC_ENABLE_TIMING=ON -DMPC_BUILD_TESTS=ON
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
python tools/benchmark_u_timing.py --dll build/Release/DISCON.dll --samples 100
python tools/analyze_u_timing.py benchmark_u_timing_samples.csv --control-period-us 10000
```

The instrumented controller writes `<output-root>.mpc_u_timing.csv` and exports `MPC_GET_LAST_TIMING` / `MPC_GET_LAST_SOLVE_INFO`. Solve modes are 0 (cold), 1 (warm), and 2 (warm failed, cold retry used). Buffered rows flush on the final callback. A DLL-only benchmark uses synthetic inputs and does not establish aeroelastic stability or shutdown performance.

Tests cover configuration validation, cold/warm solves, fallback, callback timing and actuator limits. Full OpenFAST closed-loop runs remain necessary after rebuilding.

## Files

```text
DISCON_MPC_CPP.cpp       Controller and DISCON entry point
DISCON_MPC_CPP.IN        Single named configuration
CMakeLists.txt          DLL and optional tests
GT*.csv, GF*.csv         Active SI-unit gain tables
docs/                   Gain-table provenance
tests/                  Controller regression checks
tools/                  Timing benchmark and analysis
```

Historical controller copies, build trees, archived tables, linearization decks and old timing outputs are excluded from the maintained source tree. They are not runtime dependencies.
