#
# Smallest tracking MPC example
#
# This script shows the basic public tracking MPC path:
# 1. build a tiny MTK model,
# 2. build one MPC controller,
# 3. solve one MPC step,
# 4. print the result.
#
# This file is meant to be the first tracking MPC script a new user reads.
# It is intentionally small.
# The goal is not to show every option.
# The goal is to show the shortest complete path from:
# - one MTK model,
# - to one MPC controller,
# - to one solved MPC move.

import Pkg

Pkg.activate(joinpath(@__DIR__, ".."))

using EOptInterface
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using JuMP

try
    using Ipopt
catch err
    error("tracking_mpc_demo.jl requires Ipopt in the active environment. Add it with `Pkg.add(\"Ipopt\")` and rerun.")
end

# The plant model has:
# - one state `x`,
# - one manipulated variable `u`,
# - one disturbance-like parameter `d`.
# This is enough to show the MPC workflow without extra flowsheet details.
ModelingToolkit.@parameters u = 0.0 d = 1.0
ModelingToolkit.@variables x(t)
@named sys = ODESystem([D(x) ~ -0.8 * x + u + d], t, [x], [u, d])

# Build the JuMP model used by the MPC controller.
model = Model(Ipopt.Optimizer)
set_silent(model)

# Basic MPC settings.
# `PH` is the prediction horizon length.
# `CH` is the control horizon length.
# `dt` is the sample time used inside the prediction model.
cfg = TrackingMPCConfig(
    PH = 5,
    CH = 2,
    dt = 1.0,
    integrator = "IE",
    system_kind = :ode,
    state_lower = -5.0,
    state_upper = 5.0,
)

# Build the reusable controller once.
# This step is the expensive "setup" step.
# The controller stores all JuMP variables and constraints so later solves can
# reuse the same optimization model.
ctrl = build_tracking_mpc(
    model,
    sys;
    control_specs = [
        MPCControlSpec(
            sym = sys.u,
            lower = 0.0,
            upper = 2.0,
            delta_max = 0.4,
            move_weight = 0.1,
            first_move_weight = 1.0,
        ),
    ],
    output_specs = [
        MPCOutputSpec(
            sym = sys.x,
            setpoint = 1.0,
            track_weight = 1.0,
            terminal_weight = 2.0,
            lower_soft = 0.9,
            upper_soft = 1.1,
            slack_weight = 100.0,
        ),
    ],
    stage_param_defaults = Dict(sys.d => fill(1.0, cfg.PH + 1)),
    config = cfg,
)

# Solve one MPC step from the current state and previous control.
# Here we say:
# - the plant is currently at `x = 0`;
# - the previously applied control is `u = 0`.
# The returned result contains:
# - the optimized control sequence;
# - the predicted state trajectory;
# - and scalar objective terms.
result = solve_tracking_mpc!(ctrl, Dict(sys.x => 0.0), Dict(sys.u => 0.0))

println("status: ", result.status)
println("u sequence: ", result.controls[sys.u])
println("x prediction: ", result.predictions[sys.x])
println("metrics: ", result.metrics)
