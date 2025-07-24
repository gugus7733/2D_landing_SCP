# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

# User instructions :
When answering a programming question :
You are an AI programming assistant. Follow the user's requirements carefully and to the letter. First, think step-by-step and describe your plan for what to build, written out in great detail. Then, output the code in a single code block. Minimize any other prose. Make your code human readable, dont save 3 characters on a variable name if it renders it unintelligible.
Make pretty code. You dont need smileys or comments every line, but align series of statements, use explicit variable names, and make a breathable code.

When answering non technical / non programming questions or day to day queries : answer with as little verbose as possible, to extract only the useful information. And avoid smileys, as they don't add much. Only use chapters, list and dots, as well as bold text when it is needed the most, not every other sentence. Do not talk to me if I was a toddler as well please. And also try to avoid as much as possible the confirmation biais and congratulating every sentence of the user like LLMs tend yo do : do not try to be gentle, be factual.
You can crack a joke or two, but try to keep a professional tone.

Here is a list of rules to follow when giving code in MATLAB :
Please respect these rules as much as possible.

- If you see any incoherence or bug that needs fixing on the way, even if it is independent of the task given, fix it robustly.
- Take into account the modularity. I want to be able to plug in or out any part of any of the sections of the pipeline without breaking the rest. Use parameters as much as possible. No magic numbers or hard coded parameters, they come from where they are defined and passed as parameters throughout the pipeline. Some functions will have self defined parameters that do not need adjustments in the short term, but if it is bound to be changed, it needs to be defined in the main config file and passed as parameter throughout the pipeline to the user function.
- Make your code human readable, dont save 3 characters on a variable name if it renders it unintelligible (instead of tstominc, use trust_omeg_increase, instead of scpsolhs, use scp_sol_history).
- Make pretty code. You dont need smileys or comments every line, but align series of statements, use explicit variable names, and make a breathable code. Make clear defined sections, I don't mind cartridges of comments and I like aligned '=' signs and # comments.
- Whenever you feel like i have to tell you things twice or have to repeat myself, please modify this CLAUDE.md in consequence, so that your default behavior just keeps improving (fill the # Additionnal good practices and code formalism section).
- If you find inconsistencies in the style of the code or the guidelines / best practices, please modify the code accordingly and then modify the CLAUDE.md file to make sure the new for formalism is respected in the future. This is to make the we have coherent code throughout the project and converge towards a coherent, clear, readable, usable and maintanable code from top to bottom. (for example it could be : prefer using " "" " instead of " '' ", or that all functions should list each argument one one new line, etc...)
- If you see spaces, comments or new lines that seem to be here for a reason (for example in variables definitions, functions calls or definitions), please let them be, they make the code more readable.
- Most of the time, even in a debugging scenario, you have to avoid logging stuff in the console. If necessary, save the values in a structure that is readable later, or print the output of a large computation / simulation. The console space is precious, I am tired of running after your spamming fprintf all around the code.
- Parameters are defined once, then passed along the pipeline until the user function. Assume the user function will always receive its arguments. Do not use blocks like "if parameter exists then read it, otherwise redefine it with random obscure default values". Just use the damn interface that is well defined and tracable.
- When adapting code from a previous implementation, do not modify the previous logic or values, try to preserve them as much as possible. Any deviation from the old implementation will then have to be patched later, which is even more work down the line !


# Additionnal good practices and code formalism : 
- Function parameters should be documented with clear INPUTS/OUTPUTS sections
- For physics equations, you can keep a short yet explicit name, like for example F is for force M for moment, v for speed, psi tet phi eta delta for angles, etc.
This allows much more readable equations, which is opposite to the more classical code where explicit is better. For example, rotation matrices are always T_A_B = rotation matrix from reference B to reference A. This might seems unreadable but this is very efficient for physics equations code. You have to find the right equilibrium between dark implicit variable and clear short physics convetions based names, this is not an easy exercize as sometimes the two are hard to decypher !
- Parameter name consistency must be maintained across all functions
- Parameters are passed explicitly through the call chain - no global access or default fallbacks
- State vectors should be managed through state_manager_2d() utility functions when possible
- All console output (fprintf) should be removed from core functions - use logging structures instead
- Avoid using fprintf at all costs !
- Trust regions should persist between successive SCP calls for improved convergence
- Reference trajectories must be properly time-shifted using actual elapsed time, not assumed timesteps
- Mode transitions (normal ↔ fine computation) require special handling to maintain trajectory continuity


## Project Overview

This is a 2D rocket soft landing simulation using Model Predictive Control (MPC) with Successive Convexification (SCP). The project simulates a Falcon 9-like rocket terminal landing phase with realistic aerodynamics, thrust vector control, and fuel consumption.

## Running the Simulation

**Main execution:**
```matlab
main_mpc_2d.m
```

This is the primary script that runs the complete 2D MPC simulation. It can be run in three modes:
- **MPC Control (default):** Uses SCP optimization for trajectory planning
- **Manual Control:** Set `disable_control = true` in main_mpc_2d.m for predefined thrust/gimbal profiles
- **Control Replay:** Set `replay_control = true` and `t_replay_control = X` to replay previously recorded controls until time X

**Control Replay System:**
The replay system allows recording control inputs during a simulation run and replaying them in subsequent runs:
1. **First run (recording):** Set `replay_control = false`, run simulation normally. Controls are automatically saved to `control_replay.mat`
2. **Replay run:** Set `replay_control = true` and `t_replay_control = X` to replay recorded controls until time X seconds
3. **Hybrid mode:** After replay time expires, simulation switches to normal MPC mode

Replay parameters in `main_mpc_2d.m`:
- `replay_control = false/true` - Enable/disable replay mode  
- `t_replay_control = X` - Time (seconds) until which replay is active
- `control_replay_filename = 'filename.mat'` - File to save/load control log

**Dependencies:**
- Requires `FigureManager` class (external dependency, not in codebase)
- MATLAB Optimization Toolbox (`quadprog`)
- Aerospace Toolbox (`atmosisa`)

## Core Architecture

### State Vector
7-DOF state: `[x, y, vx, vy, theta, omega, mass]`
- Position: (x,y) in inertial frame
- Velocity: (vx,vy) in inertial frame  
- Attitude: theta (pitch angle), omega (angular velocity)
- Mass: current vehicle mass

### Control Vector
2-DOF control: `[thrust, gimbal_angle]`
- Thrust magnitude (constrained by T_min/T_max)
- Engine gimbal angle delta (±8° limit)

### Key Files and Functions

**Core simulation files:**
- `main_mpc_2d.m` - Main simulation script with MPC, manual, and replay modes
- `run_scp_2d.m` - SCP optimization solver
- `simulate_step_2d.m` - Nonlinear dynamics propagation
- `build_subproblem_2d.m` - Convex QP formulation
- `get_initial_reference_2d.m` - Reference trajectory initialization
- `control_replay_utils.m` - Utility functions for control replay system

**Support files:**
- `visualize_results_2d.m` - Simulation results visualization
- `validate_linearization.m` - Linearization accuracy validation
- `computeInitialSlope.m` - Initial attitude computation

Legacy folder contains the files from an older implementation. Do not modify these files.
You can draw inspiration from them for the implementation.



