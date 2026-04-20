#
# NDMC conductivity MPC example
#
# This script shows a larger but still public MPC example:
# 1. build the NDMC plant,
# 2. build one reusable tracking MPC controller,
# 3. run a closed-loop simulation with callbacks,
# 4. save the results to CSV files.
#
# This file is intentionally larger than `tracking_mpc_demo.jl`.
# It is meant to be the "next step" example after the smallest demo.
# The extra size comes from three things:
# - a larger plant model,
# - a real closed-loop simulation,
# - logging and CSV export for later analysis.
#
# Live status output comes in two layers:
# - the generic `solve_tracking_mpc!(...; show_status=true)` line from the
#   package itself,
# - and an NDMC-specific line that adds disturbance and prediction details.

import Pkg

Pkg.activate(@__DIR__)
Pkg.instantiate()

using EOptInterface
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using JuMP, Ipopt
using OrdinaryDiffEq, DiffEqCallbacks
using DataFrames, CSV
using Printf

const NDMC_LEGACY_DT = 20.0
const NDMC_LEGACY_PSPAN = 400.0
const NDMC_LEGACY_MSPAN = 60.0
const NDMC_LEGACY_PH = Int(round(NDMC_LEGACY_PSPAN / NDMC_LEGACY_DT))
const NDMC_LEGACY_CH = Int(round(NDMC_LEGACY_MSPAN / NDMC_LEGACY_DT))

const NDMC_M_N = 14.0067
const NDMC_R_AOMAX = 0.67
const NDMC_K_OAO = 0.3
const NDMC_X_AO = 0.505
const NDMC_PHI_OAO = 2.5
const NDMC_SOTE = 0.1
const NDMC_COS = 9.1
const NDMC_VOLUME_TOTAL = 1000.0
const NDMC_ZONE_VOLUME = NDMC_VOLUME_TOTAL / 4.0
const NDMC_FLOW = NDMC_VOLUME_TOTAL / 240.0 / 4.0
const NDMC_LA0 = 149.6
const NDMC_A = 60.2
const NDMC_B = 0.229
const NDMC_SMOOTH_POS_EPS = 1e-8

const GENERATED_DIR = joinpath(@__DIR__, "generated")

# Smooth positive-part helper. This avoids a hard kink at zero.
smooth_nonnegative(x) = 0.5 * (x + sqrt(x * x + NDMC_SMOOTH_POS_EPS^2))

"""
    ndmc_rate_term(cO)

Return the approximate nitrification rate used by the NDMC example.

This helper isolates the oxygen-dependent biology term so the plant equations
read more clearly. New readers can treat this as the "how fast do reactions
consume ammonium right now?" calculation.
"""
function ndmc_rate_term(cO)
    cO_pos = smooth_nonnegative(cO)
    return NDMC_R_AOMAX * (cO_pos / (NDMC_K_OAO + cO_pos)) * NDMC_X_AO * 0.001 / 60.0 / NDMC_M_N
end

"""
    ndmc_conductivity_sink(cO)

Convert the nitrification rate into the conductivity sink term used in the
conductivity balances.

The NDMC states are conductivity-like quantities, so the biological reaction
rate is converted into the matching conductivity units here.
"""
function ndmc_conductivity_sink(cO)
    rate = ndmc_rate_term(cO)
    return (NDMC_LA0 - (NDMC_A + NDMC_B * NDMC_LA0) * sqrt(rate)) * rate * 1e3
end

"""
    ndmc_oxygen_sink(cO)

Return the oxygen-consumption term used in the dissolved-oxygen balance.

This uses the same biological rate logic as `ndmc_rate_term(...)`, but keeps
the result in oxygen units for the aeration state equation.
"""
function ndmc_oxygen_sink(cO)
    cO_pos = smooth_nonnegative(cO)
    return (NDMC_R_AOMAX * (cO_pos / (NDMC_K_OAO + cO_pos))) * NDMC_PHI_OAO * NDMC_X_AO / 60.0
end

"""
    ndmc_k_vector()

Return the calibrated mixing/transport coefficients used by the NDMC plant.

The script stores the high- and medium-loading values from the legacy case and
then averages them to get one public demonstration parameter set.
"""
function ndmc_k_vector()
    k1_hi = 2.37570423107917e-03
    k2_hi = 1.41065048786101e-03
    k3_hi = 1.50410563455869e-03
    k4_hi = 9.47661458622374e-01

    k1_me = 1.38465080722995e-03
    k2_me = 2.91428968468164e-03
    k3_me = 2.57662674484102e-03
    k4_me = 1.89896205118280e+00

    k1 = 0.5 * (k1_hi + k1_me)
    k2 = 0.5 * (k2_hi + k2_me)
    k3 = 0.5 * (k3_hi + k3_me)
    k4 = 0.5 * (k4_hi + k4_me)
    return (1000.0 / 0.38) .* [k1, k2, k3, k4]
end

"""
    NDMCMPCConfig

Configuration bundle for the public NDMC closed-loop example.

Keeping the settings in one struct makes the example easier to read: all case
assumptions, bounds, horizons, and disturbance settings live in one place
instead of being scattered across the script.
"""
Base.@kwdef struct NDMCMPCConfig
    simulation_span::Tuple{Float64, Float64} = (0.0, 4000.0)
    dt::Float64 = NDMC_LEGACY_DT
    PH::Int = NDMC_LEGACY_PH
    CH::Int = NDMC_LEGACY_CH
    setpoint::Float64 = 280.0
    move_weight::Float64 = 0.0
    first_move_weight::Float64 = 0.0
    terminal_weight::Float64 = 0.0
    Q_init::Float64 = 168.0
    Q_min::Float64 = 0.0
    Q_max::Float64 = 800.0
    delta_max::Union{Nothing, Float64} = nothing
    use_disturbance_forecast::Bool = false
    disturbance_start::Float64 = 2100.0
    disturbance_stop::Float64 = 2250.0
    Cs::Float64 = 285.0
    Cin_spike::NTuple{3, Float64} = (285.0, 285.0, 320.0)
    save_dt::Float64 = 10.0
    initial_state::NTuple{5, Float64} = (280.0, 280.0, 280.0, 280.0, 0.0)
    state_lower::NTuple{5, Float64} = (0.0, 0.0, 0.0, 0.0, 0.0)
    state_upper::NTuple{5, Float64} = (600.0, 600.0, 600.0, 600.0, 20.0)
    ipopt_tol::Float64 = 1e-6
    ipopt_max_cpu_time::Float64 = 60.0
    show_detailed_status::Bool = true
    show_generic_status::Bool = false
    status_digits::Int = 3
    status_prefix::String = "[NDMC MPC]"
    status_stream::Symbol = :stdout
end

"""
    load_config_from_env()

Build one `NDMCMPCConfig` from environment variables.

This is a convenience wrapper for two common use cases:
- run the example with the public defaults;
- override horizons or simulation length in CI or quick local smoke tests.
"""
function load_config_from_env()
    # Read the case settings from environment variables.
    # This lets the same public example run both:
    # - as a normal script with defaults,
    # - and as a shorter smoke test in CI or quick local checks.
    dt = parse(Float64, get(ENV, "NDMC_DT", string(NDMC_LEGACY_DT)))
    pspan = parse(Float64, get(ENV, "NDMC_PSPAN", string(NDMC_LEGACY_PSPAN)))
    mspan = parse(Float64, get(ENV, "NDMC_MSPAN", string(NDMC_LEGACY_MSPAN)))
    t_end = parse(Float64, get(ENV, "NDMC_T_END", "4000.0"))
    ch = max(1, Int(round(mspan / dt)))
    ph = max(ch, Int(round(pspan / dt)))
    return NDMCMPCConfig(
        simulation_span = (0.0, t_end),
        dt = dt,
        PH = ph,
        CH = ch,
        setpoint = parse(Float64, get(ENV, "NDMC_SETPOINT", "280.0")),
        move_weight = parse(Float64, get(ENV, "NDMC_MOVE_WEIGHT", "0.0")),
        first_move_weight = parse(Float64, get(ENV, "NDMC_FIRST_MOVE_WEIGHT", "0.0")),
        terminal_weight = parse(Float64, get(ENV, "NDMC_TERMINAL_WEIGHT", "0.0")),
        Q_init = parse(Float64, get(ENV, "NDMC_Q_INIT", "168.0")),
        Q_min = parse(Float64, get(ENV, "NDMC_Q_MIN", "0.0")),
        Q_max = parse(Float64, get(ENV, "NDMC_Q_MAX", "800.0")),
        delta_max = haskey(ENV, "NDMC_DELTA_MAX") ? parse(Float64, ENV["NDMC_DELTA_MAX"]) : nothing,
        use_disturbance_forecast = lowercase(strip(get(ENV, "NDMC_USE_DISTURBANCE_FORECAST", "false"))) in ("1", "true", "yes", "on"),
        disturbance_start = parse(Float64, get(ENV, "NDMC_DIST_START", "2100.0")),
        disturbance_stop = parse(Float64, get(ENV, "NDMC_DIST_STOP", "2250.0")),
        Cs = parse(Float64, get(ENV, "NDMC_CS", "285.0")),
        Cin_spike = (
            parse(Float64, get(ENV, "NDMC_CIN1", "285.0")),
            parse(Float64, get(ENV, "NDMC_CIN2", "285.0")),
            parse(Float64, get(ENV, "NDMC_CIN3", "320.0")),
        ),
        save_dt = parse(Float64, get(ENV, "NDMC_SAVE_DT", string(dt))),
        initial_state = (
            parse(Float64, get(ENV, "NDMC_C1_INIT", "280.0")),
            parse(Float64, get(ENV, "NDMC_C2_INIT", "280.0")),
            parse(Float64, get(ENV, "NDMC_C3_INIT", "280.0")),
            parse(Float64, get(ENV, "NDMC_CMIX_INIT", "280.0")),
            parse(Float64, get(ENV, "NDMC_CO_INIT", "0.0")),
        ),
        state_lower = (
            parse(Float64, get(ENV, "NDMC_C_LO", "0.0")),
            parse(Float64, get(ENV, "NDMC_C_LO", "0.0")),
            parse(Float64, get(ENV, "NDMC_C_LO", "0.0")),
            parse(Float64, get(ENV, "NDMC_CMIX_LO", "0.0")),
            parse(Float64, get(ENV, "NDMC_CO_LO", "0.0")),
        ),
        state_upper = (
            parse(Float64, get(ENV, "NDMC_C_HI", "600.0")),
            parse(Float64, get(ENV, "NDMC_C_HI", "600.0")),
            parse(Float64, get(ENV, "NDMC_C_HI", "600.0")),
            parse(Float64, get(ENV, "NDMC_CMIX_HI", "600.0")),
            parse(Float64, get(ENV, "NDMC_CO_HI", "20.0")),
        ),
        ipopt_tol = parse(Float64, get(ENV, "NDMC_IPOPT_TOL", "1e-6")),
        ipopt_max_cpu_time = parse(Float64, get(ENV, "NDMC_IPOPT_MAX_CPU_TIME", "60.0")),
        show_detailed_status = lowercase(strip(get(ENV, "NDMC_SHOW_DETAILED_STATUS", "true"))) in ("1", "true", "yes", "on"),
        show_generic_status = lowercase(strip(get(ENV, "NDMC_SHOW_GENERIC_STATUS", "false"))) in ("1", "true", "yes", "on"),
        status_digits = parse(Int, get(ENV, "NDMC_STATUS_DIGITS", "3")),
        status_prefix = get(ENV, "NDMC_STATUS_PREFIX", "[NDMC MPC]"),
        status_stream = begin
            stream = Symbol(lowercase(strip(get(ENV, "NDMC_STATUS_STREAM", "stdout"))))
            stream in (:stdout, :stderr) || error("NDMC_STATUS_STREAM must be `stdout` or `stderr`.")
            stream
        end,
    )
end

"""
    _ndmc_status_io(cfg)

Return the IO stream used by the NDMC live-status printer.
"""
function _ndmc_status_io(cfg::NDMCMPCConfig)
    if cfg.status_stream === :stdout
        return stdout
    elseif cfg.status_stream === :stderr
        return stderr
    end
    error("Unsupported NDMC status stream $(cfg.status_stream). Use :stdout or :stderr.")
end

"""
    _fmt_ndmc_scalar(x; digits=3)

Format one scalar for NDMC live-status output.
"""
function _fmt_ndmc_scalar(x; digits::Int = 3)
    if x isa Real && isfinite(float(x))
        return @sprintf("%.*f", digits, float(x))
    end
    return "NaN"
end

"""
    print_ndmc_mpc_status(io, cfg, sys, result, state_values, q_prev, q_apply, cin_now, t_now, step)

Print one NDMC-specific live status line after an MPC trigger.

This complements the generic `show_status=true` support in
`solve_tracking_mpc!(...)` by adding scenario-specific fields such as the
influent shock values and the next/terminal `C3` prediction.
"""
function print_ndmc_mpc_status(io::IO,
                               cfg::NDMCMPCConfig,
                               sys,
                               result,
                               state_values::AbstractDict,
                               q_prev::Real,
                               q_apply::Real,
                               cin_now,
                               t_now::Real,
                               step::Integer)
    accepted = is_accepted_mpc_status(result.status; accepted_statuses = default_mpc_accepted_statuses())
    pred_c3 = get(result.predictions, sys.C3, Float64[])
    pred_c3_next = length(pred_c3) >= 2 ? pred_c3[2] : NaN
    pred_c3_terminal = isempty(pred_c3) ? NaN : pred_c3[end]

    parts = String[
        cfg.status_prefix,
        "step=" * string(step),
        "t=" * _fmt_ndmc_scalar(t_now; digits = cfg.status_digits),
        "status=" * string(result.status),
        "accepted=" * string(accepted),
        "obj=" * _fmt_ndmc_scalar(get(result.metrics, :objective, NaN); digits = cfg.status_digits),
        "C1=" * _fmt_ndmc_scalar(get(state_values, sys.C1, NaN); digits = cfg.status_digits) * " (sp=" * _fmt_ndmc_scalar(cfg.setpoint; digits = cfg.status_digits) * ")",
        "C2=" * _fmt_ndmc_scalar(get(state_values, sys.C2, NaN); digits = cfg.status_digits) * " (sp=" * _fmt_ndmc_scalar(cfg.setpoint; digits = cfg.status_digits) * ")",
        "C3=" * _fmt_ndmc_scalar(get(state_values, sys.C3, NaN); digits = cfg.status_digits) * " (sp=" * _fmt_ndmc_scalar(cfg.setpoint; digits = cfg.status_digits) * ")",
        "Q_prev=" * _fmt_ndmc_scalar(q_prev; digits = cfg.status_digits),
        "Q_apply=" * _fmt_ndmc_scalar(q_apply; digits = cfg.status_digits),
        "Cin1=" * _fmt_ndmc_scalar(cin_now[1]; digits = cfg.status_digits),
        "Cin2=" * _fmt_ndmc_scalar(cin_now[2]; digits = cfg.status_digits),
        "Cin3=" * _fmt_ndmc_scalar(cin_now[3]; digits = cfg.status_digits),
        "pred_C3_next=" * _fmt_ndmc_scalar(pred_c3_next; digits = cfg.status_digits),
        "pred_C3_terminal=" * _fmt_ndmc_scalar(pred_c3_terminal; digits = cfg.status_digits),
    ]

    println(io, string(first(parts), " ", join(parts[2:end], " | ")))
    flush(io)
    return nothing
end

@mtkmodel NDMCPlant begin
    @parameters begin
        Cin1 = 285.0
        Cin2 = 285.0
        Cin3 = 320.0
        Q_air = 168.0
        k1 = 1.0
        k2 = 1.0
        k3 = 1.0
        k4 = 1.0
    end
    @variables begin
        C1(t) = 280.0
        C2(t) = 280.0
        C3(t) = 280.0
        Cmix(t) = 280.0
        cO(t) = 0.0
    end
    @equations begin
        D(C1) ~ (k1 * (Cmix - C1) + NDMC_FLOW * Cin1 - NDMC_FLOW * C1) / NDMC_ZONE_VOLUME - ndmc_conductivity_sink(cO)
        D(C2) ~ (k2 * (Cmix - C2) + NDMC_FLOW * Cin2 - NDMC_FLOW * C2) / NDMC_ZONE_VOLUME - ndmc_conductivity_sink(cO)
        D(C3) ~ (k3 * (Cmix - C3) + NDMC_FLOW * Cin3 - NDMC_FLOW * C3) / NDMC_ZONE_VOLUME - ndmc_conductivity_sink(cO)
        D(Cmix) ~ (k4 * (C1 + C2 + C3 - 3.0 * Cmix)) / NDMC_ZONE_VOLUME - ndmc_conductivity_sink(cO)
        D(cO) ~ ((NDMC_SOTE * 0.2967 * Q_air) / NDMC_COS / NDMC_VOLUME_TOTAL) * (NDMC_COS - cO) - ndmc_oxygen_sink(cO)
    end
end

"""
    build_ndmc_system(cfg)

Construct the ModelingToolkit system used by both the simulated plant and the
prediction model.

The example deliberately uses the same equations for plant and controller so
new users can focus on the MPC workflow instead of plant-model mismatch.
"""
function build_ndmc_system(cfg::NDMCMPCConfig)
    # Build the MTK plant used by the controller and the ODE simulation.
    # The same MTK system is used in two roles:
    # - as the prediction model inside MPC,
    # - and as the simulated plant in the ODE solve.
    k = ndmc_k_vector()
    @mtkbuild sys = NDMCPlant(
        Cin1 = cfg.Cs,
        Cin2 = cfg.Cs,
        Cin3 = cfg.Cs,
        Q_air = cfg.Q_init,
        k1 = k[1],
        k2 = k[2],
        k3 = k[3],
        k4 = k[4],
    )
    return sys
end

"""
    ndmc_initial_state(sys, cfg)

Return the initial state dictionary in ModelingToolkit key form.

The controller builders in `EOptInterface` expect state updates keyed by the
symbolic MTK states, so this helper performs that small translation once.
"""
function ndmc_initial_state(sys, cfg::NDMCMPCConfig)
    # Build the initial state dictionary in MTK key form.
    # Using MTK keys here keeps the script consistent with the public MPC API.
    return Dict(
        sys.C1 => cfg.initial_state[1],
        sys.C2 => cfg.initial_state[2],
        sys.C3 => cfg.initial_state[3],
        sys.Cmix => cfg.initial_state[4],
        sys.cO => cfg.initial_state[5],
    )
end

"""
    disturbance_triplet(t_now, cfg)

Return the influent conductivity values active at time `t_now`.

Outside the disturbance window the plant stays at its baseline feed. Inside the
window, the configured shock values are applied.
"""
function disturbance_triplet(t_now::Real, cfg::NDMCMPCConfig)
    # Return the current influent conductivity triplet.
    # Outside the disturbance window the plant stays at the baseline feed.
    if cfg.disturbance_start <= t_now < cfg.disturbance_stop
        return cfg.Cin_spike
    end
    return (cfg.Cs, cfg.Cs, cfg.Cs)
end

"""
    build_ndmc_controller(sys; cfg)

Build the reusable NDMC tracking MPC controller.

This is the one-time setup step. It creates the JuMP model, registers the plant
equations, and adds the tracking objective and control constraints. The online
loop later reuses this controller at every trigger.
"""
function build_ndmc_controller(sys; cfg::NDMCMPCConfig)
    # Build the reusable tracking MPC controller once.
    # This is the one-time setup step.
    # After this, online MPC solves only update values and re-optimize.
    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "tol", cfg.ipopt_tol)
    set_optimizer_attribute(model, "max_cpu_time", cfg.ipopt_max_cpu_time)
    set_optimizer_attribute(model, "print_level", 0)
    set_silent(model)

    mpc_cfg = TrackingMPCConfig(
        PH = cfg.PH,
        CH = cfg.CH,
        dt = cfg.dt,
        integrator = "IE",
        system_kind = :ode,
        state_lower = Dict(
            sys.C1 => cfg.state_lower[1],
            sys.C2 => cfg.state_lower[2],
            sys.C3 => cfg.state_lower[3],
            sys.Cmix => cfg.state_lower[4],
            sys.cO => cfg.state_lower[5],
        ),
        state_upper = Dict(
            sys.C1 => cfg.state_upper[1],
            sys.C2 => cfg.state_upper[2],
            sys.C3 => cfg.state_upper[3],
            sys.Cmix => cfg.state_upper[4],
            sys.cO => cfg.state_upper[5],
        ),
        rhs0 = 0.0,
        track_start_index = 2,
        lower_clip = 0.0,
    )

    ctrl = build_tracking_mpc(
        model,
        sys;
        control_specs = [
            MPCControlSpec(
                sym = sys.Q_air,
                lower = cfg.Q_min,
                upper = cfg.Q_max,
                delta_max = cfg.delta_max,
                move_weight = cfg.move_weight,
                first_move_weight = cfg.first_move_weight,
                base_name = "Q_air",
            ),
        ],
        output_specs = [
            MPCOutputSpec(sym = sys.C1, setpoint = cfg.setpoint, track_weight = 1.0, terminal_weight = cfg.terminal_weight, base_name = "C1"),
            MPCOutputSpec(sym = sys.C2, setpoint = cfg.setpoint, track_weight = 1.0, terminal_weight = cfg.terminal_weight, base_name = "C2"),
            MPCOutputSpec(sym = sys.C3, setpoint = cfg.setpoint, track_weight = 1.0, terminal_weight = cfg.terminal_weight, base_name = "C3"),
        ],
        stage_param_defaults = Dict(
            sys.Cin1 => fill(cfg.Cs, cfg.PH + 1),
            sys.Cin2 => fill(cfg.Cs, cfg.PH + 1),
            sys.Cin3 => fill(cfg.Cs, cfg.PH + 1),
        ),
        config = mpc_cfg,
    )

    return ctrl
end

"""
    set_plant_disturbance!(integ, sys, cfg, t_now)

Write the currently active disturbance into the simulated ODE integrator.

This updates the plant model only. The MPC disturbance preview is handled
separately by `update_disturbance_forecast!(...)`.
"""
function set_plant_disturbance!(integ, sys, cfg::NDMCMPCConfig, t_now::Real)
    # Write the current disturbance into the ODE integrator parameters.
    # This changes the simulated plant only.
    # The controller preview is updated separately.
    cin1, cin2, cin3 = disturbance_triplet(t_now, cfg)
    integ.ps[sys.Cin1] = cin1
    integ.ps[sys.Cin2] = cin2
    integ.ps[sys.Cin3] = cin3
    return nothing
end

"""
    update_disturbance_forecast!(ctrl, sys, cfg, t_now)

Update the stage-by-stage disturbance preview seen by the controller.

If `cfg.use_disturbance_forecast` is `false`, the controller assumes the
current disturbance will stay constant across the horizon. If it is `true`,
this helper builds a full preview over the prediction horizon.
"""
function update_disturbance_forecast!(ctrl::TrackingMPCController, sys, cfg::NDMCMPCConfig, t_now::Real)
    # Update the fixed stage preview used by the MPC model.
    # If `use_disturbance_forecast=false`, the controller only sees the current
    # disturbance value and holds it constant across the horizon.
    # If it is true, the controller gets a stage-by-stage preview.
    cin_now = disturbance_triplet(t_now, cfg)
    for (sym, idx) in zip((sys.Cin1, sys.Cin2, sys.Cin3), 1:3)
        vals = if cfg.use_disturbance_forecast
            [
                disturbance_triplet(t_now + (k - 1) * cfg.dt, cfg)[idx]
                for k in 1:ctrl.N
            ]
        else
            fill(cin_now[idx], ctrl.N)
        end
        update_stage_parameter!(ctrl, sym, vals)
    end
    return ctrl
end

"""
    solve_ndmc_step!(ctrl, sys, cfg, state_values, q_prev, t_now)

Solve one MPC step for the NDMC example.

This wrapper performs the NDMC-specific online work in the correct order:
1. refresh the disturbance preview,
2. solve the generic tracking MPC problem,
3. optionally print generic and/or NDMC-specific live status lines,
4. return the first control move that should be applied to the plant.
"""
function solve_ndmc_step!(ctrl::TrackingMPCController,
                          sys,
                          cfg::NDMCMPCConfig,
                          state_values::AbstractDict,
                          q_prev::Real,
                          t_now::Real;
                          step::Union{Nothing, Integer} = nothing)
    # Update the preview, solve the MPC, and return the first move.
    # This helper is the NDMC-specific wrapper around the generic public
    # function `solve_tracking_mpc!(...)`.
    update_disturbance_forecast!(ctrl, sys, cfg, t_now)
    cin_now = disturbance_triplet(t_now, cfg)
    result = solve_tracking_mpc!(
        ctrl,
        state_values,
        Dict(sys.Q_air => q_prev);
        lower_clip = 0.0,
        show_status = cfg.show_generic_status,
        status_io = _ndmc_status_io(cfg),
        status_time = t_now,
        status_output_syms = [sys.C1, sys.C2, sys.C3],
        status_control_syms = [sys.Q_air],
        status_digits = cfg.status_digits,
        status_prefix = cfg.status_prefix,
    )
    q_apply = result.controls[sys.Q_air][1]
    if abs(q_apply - q_prev) < 1e-8
        q_apply = q_prev
    end
    if cfg.show_detailed_status
        print_ndmc_mpc_status(
            _ndmc_status_io(cfg),
            cfg,
            sys,
            result,
            state_values,
            q_prev,
            q_apply,
            cin_now,
            t_now,
            something(step, -1),
        )
    end
    return (q_apply = q_apply, result = result)
end

"""
    metric_column(logctx, key)

Return one logged metric column, or `NaN` values if that metric was not stored.

This keeps the CSV export section simple and avoids repeated missing-key checks.
"""
function metric_column(logctx::MPCLog, key::Symbol)
    return get(logctx.Metrichist, key, fill(NaN, length(logctx.pred_times)))
end

"""
    main()

Run the full public NDMC closed-loop example and save the results to CSV.

This is the script-level entry point. It builds the plant and controller,
runs the receding-horizon simulation with callbacks, and exports the histories
used by the notebooks and plots.
"""
function main()
    mkpath(GENERATED_DIR)

    # Build the plant, the initial state, and the controller.
    cfg = load_config_from_env()
    sys = build_ndmc_system(cfg)
    state0 = ndmc_initial_state(sys, cfg)
    ctrl = build_ndmc_controller(sys; cfg = cfg)
    mpc_step = Ref(0)

    # Solve once at the initial state so the plant starts with a consistent
    # first applied move and the logs already contain one prediction snapshot.
    initial = solve_ndmc_step!(ctrl, sys, cfg, state0, cfg.Q_init, cfg.simulation_span[1]; step = mpc_step[])
    mpc_step[] += 1
    initial_q = initial.q_apply
    initial_cin = disturbance_triplet(cfg.simulation_span[1], cfg)

    track_keys = Dict(sym => Symbol("track_" * ctrl.output_names[sym]) for sym in (sys.C1, sys.C2, sys.C3))
    move_key = Symbol("move_" * ctrl.control_names[sys.Q_air])
    move1_key = Symbol("move1_" * ctrl.control_names[sys.Q_air])

    logctx = make_mpc_log(
        sys;
        control_keys = [sys.Q_air],
        predicted_keys = [sys.C1, sys.C2, sys.C3, sys.Q_air],
        metric_keys = vcat(collect(values(track_keys)), [move_key, move1_key, :objective]),
    )
    record_mpc_prediction!(logctx, cfg.simulation_span[1], initial.result.predictions)
    record_mpc_metrics!(logctx, initial.result.metrics; missing = :skip)
    seed_mpc_log!(logctx, state0, cfg.simulation_span[1]; control_values = Dict(sys.Q_air => initial_q))

    # These are the parameter values carried by the simulated plant.
    p0 = Dict(
        sys.Q_air => initial_q,
        sys.Cin1 => initial_cin[1],
        sys.Cin2 => initial_cin[2],
        sys.Cin3 => initial_cin[3],
    )

    prob_cl = ODEProblem(sys, merge(state0, p0), cfg.simulation_span)

    # Separate trigger grids are used for control updates and logging. This lets
    # us log at one cadence even if we later choose a different control period.
    times_ctrl = collect((cfg.simulation_span[1] + cfg.dt):cfg.dt:cfg.simulation_span[2])
    times_log = collect((cfg.simulation_span[1] + cfg.save_dt):cfg.save_dt:cfg.simulation_span[2])
    disturbance_times = unique(filter(t -> cfg.simulation_span[1] <= t <= cfg.simulation_span[2], [cfg.disturbance_start, cfg.disturbance_stop]))

    function mpc_affect!(integ)
        # Order matters here:
        # 1. write the currently active disturbance into the plant,
        # 2. read the latest plant state,
        # 3. solve MPC,
        # 4. apply only the first move from the optimized trajectory.
        set_plant_disturbance!(integ, sys, cfg, integ.t)
        state_values = current_state_map(integ, sys)
        q_prev = integ.ps[sys.Q_air]
        solved = solve_ndmc_step!(ctrl, sys, cfg, state_values, q_prev, integ.t; step = mpc_step[])
        mpc_step[] += 1
        q_apply = solved.q_apply
        integ.ps[sys.Q_air] = q_apply
        record_mpc_prediction!(logctx, integ.t, solved.result.predictions)
        record_mpc_metrics!(logctx, solved.result.metrics; missing = :skip)
        return nothing
    end

    function disturbance_affect!(integ)
        # This callback keeps the plant disturbance synchronized at the exact
        # window boundaries even between control triggers.
        set_plant_disturbance!(integ, sys, cfg, integ.t)
        return nothing
    end

    function log_affect!(integ)
        # Log the plant state after the most recent control/disturbance updates.
        log_mpc_state!(logctx, integ; control_values = Dict(sys.Q_air => integ.ps[sys.Q_air]), missing_control = :hold)
        return nothing
    end

    mpc_cb = PresetTimeCallback(times_ctrl, mpc_affect!; save_positions = (false, false))
    dist_cb = isempty(disturbance_times) ? nothing : PresetTimeCallback(disturbance_times, disturbance_affect!; save_positions = (false, false))
    log_cb = PresetTimeCallback(times_log, log_affect!; save_positions = (false, false))
    cbs = dist_cb === nothing ? CallbackSet(mpc_cb, log_cb) : CallbackSet(dist_cb, mpc_cb, log_cb)

    # Run the closed-loop ODE simulation.
    solve(
        prob_cl,
        Rodas5P();
        adaptive = true,
        abstol = 1e-9,
        reltol = 1e-7,
        callback = cbs,
        tstops = sort(unique(vcat(times_ctrl, times_log, disturbance_times))),
        saveat = times_log,
    )

    # Save the state and control history.
    closed_loop_df = DataFrame(
        time = logctx.ts,
        C1 = logctx.Xhist[sys.C1],
        C2 = logctx.Xhist[sys.C2],
        C3 = logctx.Xhist[sys.C3],
        Cmix = logctx.Xhist[sys.Cmix],
        cO = logctx.Xhist[sys.cO],
        Q_air = logctx.Controlhist[sys.Q_air],
    )
    CSV.write(joinpath(GENERATED_DIR, "ndmc_conductivity_eopt_closed_loop.csv"), closed_loop_df)

    # Save the applied move and objective terms from each MPC solve.
    applied_df = DataFrame(
        time = logctx.pred_times,
        Q_applied = [traj[1] for traj in logctx.Predhist[sys.Q_air]],
        track_C1 = metric_column(logctx, track_keys[sys.C1]),
        track_C2 = metric_column(logctx, track_keys[sys.C2]),
        track_C3 = metric_column(logctx, track_keys[sys.C3]),
        move_Q_air = metric_column(logctx, move_key),
        move1_Q_air = metric_column(logctx, move1_key),
        objective = metric_column(logctx, :objective),
    )
    CSV.write(joinpath(GENERATED_DIR, "ndmc_conductivity_eopt_applied_control.csv"), applied_df)

    # Save the public JuMP base names used by this example.
    name_map_df = DataFrame(
        category = ["control", "output", "output", "output"],
        symbol = [string(sys.Q_air), string(sys.C1), string(sys.C2), string(sys.C3)],
        base_name = [
            ctrl.control_names[sys.Q_air],
            ctrl.output_names[sys.C1],
            ctrl.output_names[sys.C2],
            ctrl.output_names[sys.C3],
        ],
    )
    CSV.write(joinpath(GENERATED_DIR, "ndmc_conductivity_eopt_name_map.csv"), name_map_df)

    println("Saved NDMC example outputs to:")
    println("  ", joinpath(GENERATED_DIR, "ndmc_conductivity_eopt_closed_loop.csv"))
    println("  ", joinpath(GENERATED_DIR, "ndmc_conductivity_eopt_applied_control.csv"))
    println("  ", joinpath(GENERATED_DIR, "ndmc_conductivity_eopt_name_map.csv"))
end

main()
