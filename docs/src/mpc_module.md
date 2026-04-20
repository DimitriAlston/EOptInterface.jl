# MPC Module

`EOptInterface` exposes two public MPC workflows.

- `build_tracking_mpc(...)`
  Use this when your prediction model is a `ModelingToolkit` ODE or DAE system.

- `register_dmcsystem(...)`
  Use this when your prediction model is already a numeric step response and you want a classic DMC block.

This page is written for public users who are comfortable reading beginner-to-intermediate Julia code.

## Five-Minute Mental Model

The package separates MPC work into three layers.

### Layer 1: Model Registration

This layer turns a dynamic model into JuMP variables and constraints.

Examples:

- `register_odesystem(...)`
- `register_daesystem(...)`
- `register_dmcsystem(...)`

If you are new, think of this as:

`dynamic model` -> `optimization-ready prediction equations`

### Layer 2: Controller Builder

This layer adds the controller-specific pieces:

- manipulated inputs
- tracked outputs
- setpoints
- move penalties
- soft constraints
- control-horizon hold rules

The main public function here is:

- `build_tracking_mpc(...)`

### Layer 3: Online Solve Loop

This layer updates the reusable controller before each MPC solve:

- copy in the newest plant state
- copy in the previously applied input
- optionally update disturbance previews
- optionally update setpoints
- solve the JuMP model
- read back the first move and the predicted trajectories

The main public functions here are:

- `prepare_tracking_mpc_step!(...)`
- `solve_tracking_mpc!(...)`

## Which API Should You Use?

| Situation | Best starting API |
| --- | --- |
| You already have a `ModelingToolkit` ODE or DAE plant | `build_tracking_mpc(...)` |
| You already have a numeric step response | `register_dmcsystem(...)` |
| You want to log states, predictions, and objective terms during a closed-loop run | `make_mpc_log(...)` and the related logging helpers |

Rule of thumb:

- If your plant is mechanistic and nonlinear, start with tracking MPC.
- If your plant is already reduced to step-response data, start with DMC.

## Direct-Transcription Tracking MPC

The modern MPC path in this package is direct transcription.

That means the package:

1. creates one decision trajectory for each selected state,
2. creates one decision trajectory for each control input,
3. creates any stage-wise parameter trajectories you want to preview,
4. writes the discretized model equations into JuMP,
5. adds the tracking objective and MPC-specific constraints.

The same JuMP model is then reused at every control step.

That reuse is important.
It is the reason the online MPC loop is fast enough to be practical.
You build once, then only update values between solves.

## Public Tracking-MPC Workflow

The normal order is:

1. define the plant with `ModelingToolkit`,
2. define controls with `MPCControlSpec`,
3. define outputs with `MPCOutputSpec`,
4. choose horizons and solver settings with `TrackingMPCConfig`,
5. build the controller with `build_tracking_mpc(...)`,
6. solve online with `solve_tracking_mpc!(...)`.

## Minimal Example

```julia
using EOptInterface
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using JuMP, Ipopt

@parameters u = 0.0 d = 1.0
@variables x(t)
@named sys = ODESystem([D(x) ~ -0.8 * x + u + d], t, [x], [u, d])

model = Model(Ipopt.Optimizer)
set_silent(model)

cfg = TrackingMPCConfig(
    PH = 5,
    CH = 2,
    dt = 1.0,
    integrator = "IE",
    system_kind = :ode,
    state_lower = -5.0,
    state_upper = 5.0,
)

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

result = solve_tracking_mpc!(ctrl, Dict(sys.x => 0.0), Dict(sys.u => 0.0))
```

What comes back in `result`:

- `result.controls`
  The optimized control trajectories.

- `result.predictions`
  The predicted state trajectories.

- `result.metrics`
  Scalar objective terms, such as tracking cost and move cost.

## What `MPCControlSpec` Means

Each `MPCControlSpec` describes one manipulated variable.

Important fields:

- `sym`
  Which `ModelingToolkit` parameter or variable is the control input.

- `lower`, `upper`
  Hard bounds on the control.

- `delta_max`
  Maximum change allowed between consecutive control moves.

- `move_weight`
  Penalty on move size after the first move.

- `first_move_weight`
  Penalty on the difference between the first optimized move and the previously applied move.

## What `MPCOutputSpec` Means

Each `MPCOutputSpec` describes one tracked output.

Important fields:

- `sym`
  Which state or algebraic quantity is being tracked.

- `setpoint`
  The desired target value.

- `track_weight`
  Weight on normal tracking error over the horizon.

- `terminal_weight`
  Extra weight on the final prediction stage.

- `lower_soft`, `upper_soft`
  Optional soft-zone bounds.

- `slack_weight`
  Cost on violating the soft-zone bounds.

## What `TrackingMPCConfig` Means

This struct holds the shared controller settings.

Important fields:

- `PH`
  Prediction horizon length in stages.

- `CH`
  Number of free control moves.
  After this point the controller holds the last move constant.

- `dt`
  Sample time used by the prediction model.

- `integrator`
  Low-level transcription rule such as `"IE"` or `"RK4"`.

- `system_kind`
  `:ode` or `:dae`.

- `state_lower`, `state_upper`
  Default state bounds used when building state trajectories.

## Online Update Pattern

Inside a callback or a receding-horizon loop, the usual order is:

```julia
state_values = current_state_map(integ, sys)
update_stage_parameter!(ctrl, sys.d, disturbance_forecast)
update_tracking_targets!(ctrl; setpoints = Dict(sys.x => 1.2))
result = solve_tracking_mpc!(ctrl, state_values, Dict(sys.u => integ.ps[sys.u]))
u_apply = result.controls[sys.u][1]
```

This is the core idea:

- the controller object stays the same,
- only the data inside the controller changes.

If you want live terminal output while the loop is running, call:

```julia
result = solve_tracking_mpc!(ctrl, state_values, previous_controls; show_status = true)
```

Scenario-specific examples can layer extra status lines on top of this without
changing the package-level API. The NDMC example does exactly that.

## When To Use `prepare_tracking_mpc_step!`

`solve_tracking_mpc!(...)` is the easiest full workflow.

Use `prepare_tracking_mpc_step!(...)` when you want more manual control.

For example:

- you want to inspect the JuMP model before solving,
- you want to call `optimize!(ctrl.model)` yourself,
- you want to change solver settings between solves.

## DAE Path

Set:

```julia
TrackingMPCConfig(system_kind = :dae)
```

This routes the controller through `register_daesystem(...)` instead of `register_odesystem(...)`.

Use this when algebraic equations must remain inside the prediction model.

The smallest low-level DAE registration example is:

- `examples/dae_registration_demo.jl`

The matching notebook is:

- `notebooks/dae_registration_simple.ipynb`

## DMC Path

`register_dmcsystem(...)` is the step-response path.

Use it when:

- you already have a step response,
- you do not want to carry a full mechanistic `ModelingToolkit` model,
- you want a lighter classic DMC formulation.

The smallest example is:

- `examples/dmc_registration_demo.jl`

## Logging Helpers

For closed-loop work, the public logging helpers are:

- `make_mpc_log(...)`
- `seed_mpc_log!(...)`
- `log_mpc_state!(...)`
- `record_mpc_prediction!(...)`
- `record_mpc_metrics!(...)`

These helpers do not solve MPC by themselves.
They only record what happened during a run.

That separation is intentional.
It keeps the controller logic and the reporting logic independent.

Live status printing is separate again:

- `solve_tracking_mpc!(...; show_status=true)` prints a compact generic line;
- examples or notebooks can add scenario-specific wrappers when they need extra
  fields that are not part of the general API.

## Suggested Reading Order

If you are completely new, read in this order:

1. `examples/tracking_mpc_demo.jl`
2. this page
3. `src/trackingmpc.jl`
4. `src/mpcutils.jl`
5. `examples/ndmc_conductivity_mpc_demo.jl`

If you want the DMC path instead:

1. `examples/dmc_registration_demo.jl`
2. `src/dmcutils.jl`

If you want the low-level model-registration layer:

1. `src/userfuncs.jl`
2. `examples/ode_model.jl`
3. `examples/dae_registration_demo.jl`

## Practical Notes

- `CH` means the number of free moves, not the number of total stages.
- Soft bounds only matter when `slack_weight > 0`.
- Stage parameters are stored as fixed JuMP variables so they can be updated online without rebuilding the optimization model.
- The logging helpers work for both direct-transcription MPC and DMC.

## One-Sentence Summary

If you have a mechanistic `ModelingToolkit` plant, start with `build_tracking_mpc(...)`; if you have a numeric step response, start with `register_dmcsystem(...)`; in both cases, the package is designed so you build once and update online many times.
