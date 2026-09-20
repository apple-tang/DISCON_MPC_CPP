# Cleanup validation, 2026-09-20

The maintained controller was reduced to the current fixed-reference, rate-limited actuator path, with cold/warm QP initialization and optional debug/timing output. The active six gain tables and the supplied numerical tuning are unchanged.

## Checks performed

- MSVC 19.44 C++17 syntax check passed.
- CMake configuration and the Release `controller_tests` executable built successfully with timing enabled. The production DLL was not rebuilt or replaced.
- CTest passed: invalid/missing configuration handling, 13 malformed configuration cases, independent coupled-system LQR Bellman residual, cold/warm solves, 10 ms solve gating, magnitude/rate limits, repeated-time callbacks, forced solver failure, debug trace creation, final callback and the supplied 330/70 horizon.
- A separate synthetic 103-callback comparison against the pre-cleanup controller produced identical applied torque, pitch and return statuses. Maximum observed differences were 0 N m and 0 rad. This comparison covers baseline control, the fracture-time switch and a short post-fracture interval using the supplied cold-start tuning; it is not an OpenFAST aeroelastic run.
- The timing analysis and benchmark scripts passed Python syntax checks. Timing performance was not remeasured.
- Git confirmed no changes to the six active gain-table CSV files.

## Behavioral fixes

- Invalid configuration stops initialization without writing actuator commands; parsing exceptions are contained at the DLL boundary.
- The final controller callback flushes logs and releases solver state without solving another QP or changing outputs.
- The Riccati update now calculates `B' P` and `B' P A` in separate buffers. Previously, an in-place matrix product could reuse already-overwritten entries. Nonfinite Riccati iterates are rejected.

The original and cleaned controllers both reported 667 nonconverged fallback grid points with the supplied tuning and the existing 500-iteration budget. These points use zero feedback gain, yielding rate-limited torque release and feathering on QP failure. This is not a verified stabilizing LQR backup.

## Before a physical comparison

Rebuild the DLL, deploy it together with the new named IN file, and rerun the intended modified-OpenFAST cases. Changing solver mode or tuning also needs closed-loop verification. No stability or real-time performance claim is inferred from the synthetic tests.
