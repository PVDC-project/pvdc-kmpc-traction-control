# Model predictive traction control system based on the Koopman operator

This repository contains the MATLAB/Simulink scripts and CarMaker files for implementing the algorithm described in [TODO](https://ieeexplore.ieee.org/document/10308482).

The simulation compares the Koopman operator-based MPC controller (KMPC) with a standard nonlinear MPC controller (NMPC) and a PID controller.

Before starting the experiments, the linear model should be identified by running [setup/kmpc_sysid.m](https://github.com/PVDC-project/pvdc-kmpc-traction-control/blob/master/setup/kmpc_sysid.m).

Required software:
- [acados](https://github.com/acados/acados) for generating the dataset for predictor identification (also possible with `ode45` in MATLAB, but much slower); last tested with [v0.3.5](https://github.com/acados/acados/releases/tag/v0.3.5), possible interface changes
- [YALMIP](https://yalmip.github.io/) for formulating the KMPC problem
- [FORCESPRO](https://forces.embotech.com/Documentation/index.html) for running NMPC (for KMPC also other QP solvers can be used)
- `CarMaker` for high-fidelity simulations

