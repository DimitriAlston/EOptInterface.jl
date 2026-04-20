module EOptInterface

# This is the package entry point.
# It does three simple things:
# 1. load the source files with `include(...)`;
# 2. decide which names are public through `export`;
# 3. present one cleaner surface to users, examples, and notebooks.
#
# If you are new to the package, you usually do not start by reading this file
# line by line. Instead, use it as a map:
# - `userfuncs.jl` contains the low-level ODE/DAE registration layer;
# - `mpcutils.jl` contains naming, warm start, logging, and debugging helpers;
# - `trackingmpc.jl` contains the high-level tracking MPC builder and online
#   solve loop;
# - `dmcutils.jl` contains the lighter DMC path based on step responses.
#
# The export list below is important because it tells you which functions are
# meant to be used directly by outside users. Many internal helpers are not
# exported on purpose.

# ---- imports ----
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

using JuMP

# ---- includes ----
include("basefuncs.jl")
include("userfuncs.jl")
include("mpcutils.jl")
include("trackingmpc.jl")
include("dmcutils.jl")

# ---- exports ----
export build_state_trajs_from_vars!, build_tracking_mpc, check_ic_sync,
       canonical_system_parameter, canonical_system_symbol,
       canonical_system_unknown, canonicalize_system_symbols,
       collect_subsystems, compute_step_prediction_errors, current_state_map,
       control_traj,
       decision_vars, default_mpc_accepted_statuses, dump_all_constraints,
       display_mpc_name,
       dump_single_var_affines, DMCRegistration, find_conflicts,
       find_constraints_with_vars, full_solutions, get_ic_constraint!,
       group_state_symbols, is_accepted_mpc_status, log_mpc_state!,
       make_mpc_log, make_state_index_map, MPCControlSpec, MPCLog,
       MPCOutputSpec, objective_value_or_nan, pretty_subsystems,
       prepare_tracking_mpc_step!, print_conflicts, print_tracking_status,
       record_mpc_metrics!,
       record_mpc_prediction!, register_daesystem, register_dmcsystem,
       register_nlsystem, register_odesystem, seed_mpc_log!,
       resolve_mpc_base_name, sanitize_mpc_name,
       set_first_move_bounds!, shift_warm_start!, solve_tracking_mpc!,
       stage_param_traj, state_traj,
       split_mtk_state_path, state_trajectory_base_name,
       summarize_dmc_registration, summarize_ic_uniqueness,
       summarize_mpc_solve, sync_state_trajectories!, TrackingMPCConfig,
       TrackingMPCController, update_dmc_state!, update_stage_parameter!,
       update_tracking_targets!, warm_start_with_constant!, worst_residuals,
       zero_small_values!

end # module
