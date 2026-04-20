# MPC Module

`EOptInterface` now exposes two reusable MPC paths distilled from the active
`Research` codebase:

- `build_tracking_mpc(...)` for direct-transcription ODE or DAE MPC built from
  ModelingToolkit systems.
- `register_dmcsystem(...)` for classic DMC blocks driven by numeric step
  responses.

The direct-transcription layer is the common pattern behind the modern NDMC,
ASM3, and ADM1 controllers: build one JuMP model once, update initial
conditions and previews online, warm-start, then re-solve at each control step.

## Which API to use?

| Use case | Primary API |
| --- | --- |
| Nonlinear or mechanistic MTK plant model | `build_tracking_mpc`, `prepare_tracking_mpc_step!`, `solve_tracking_mpc!` |
| Linear step-response controller / DMC | `register_dmcsystem`, `update_dmc_state!` |
| Closed-loop logging | `make_mpc_log`, `seed_mpc_log!`, `log_mpc_state!`, `record_mpc_prediction!`, `record_mpc_metrics!` |

## Direct-Transcription Workflow

1. Define a ModelingToolkit plant model.
2. Choose manipulated inputs with `MPCControlSpec`.
3. Choose tracked outputs and optional soft-zone bounds with `MPCOutputSpec`.
4. Pass known disturbance previews or other stage-wise parameters through
   `stage_param_defaults`.
5. At run time, update initial conditions, previews, and setpoints, then solve.

The builder automatically creates:

- one trajectory per MTK state,
- one trajectory per manipulated input,
- first-move anchors against the previously applied input,
- optional move-delta constraints,
- output tracking penalties and optional soft-zone slacks,
- control-horizon hold constraints (`u[k] = u[CH]` for `k > CH`).

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

The runnable version of this example lives in `examples/tracking_mpc_demo.jl`.
The larger conductivity-surrogate case-study example lives in
`examples/ndmc_conductivity_mpc_demo.jl`, with a matching quick notebook in
`notebooks/ndmc_conductivity_mpc_simple.ipynb`.

## Online Update Pattern

Inside a closed-loop callback or receding-horizon loop, the usual order is:

```julia
state_values = current_state_map(integ, sys)
update_stage_parameter!(ctrl, sys.d, disturbance_forecast)
update_tracking_targets!(ctrl; setpoints = Dict(sys.x => 1.2))
result = solve_tracking_mpc!(ctrl, state_values, Dict(sys.u => integ.ps[sys.u]))
u_apply = result.controls[sys.u][1]
```

If you want to split preparation from solve, call
`prepare_tracking_mpc_step!` first and `optimize!(ctrl.model)` yourself.

## DAE Path

Set `TrackingMPCConfig(system_kind = :dae)` to route the builder through
`register_daesystem(...)`. This is the same pattern needed for the ADM1-style
high-order case where algebraic constraints must stay inside the prediction
model.

The minimal low-level DAE registration example lives in
`examples/dae_registration_demo.jl`, with a matching quick notebook in
`notebooks/dae_registration_simple.ipynb`. That pair is useful when you want
to inspect the raw `register_daesystem(...)` behavior directly rather than
through `build_tracking_mpc(...)`.

## DMC Path

`register_dmcsystem(...)` stays as the lean DMC interface for step-response
controllers. See `examples/dmc_registration_demo.jl` when the controller model
is already a numeric step response rather than a ModelingToolkit system.

## Notes

- `CH` is interpreted as the number of free moves; the remaining horizon is
  held constant.
- Soft bounds require a positive `slack_weight`.
- `stage_param_defaults` are stored as fixed JuMP variables so forecasts can be
  updated without rebuilding the optimization model.
- The logging helpers are intentionally generic and work for both
  direct-transcription MPC and DMC workflows.
