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

import Pkg

Pkg.activate(@__DIR__)
Pkg.instantiate()

using EOptInterface
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using JuMP, Ipopt
using OrdinaryDiffEq, DiffEqCallbacks
using DataFrames, CSV

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

function ndmc_rate_term(cO)
    cO_pos = smooth_nonnegative(cO)
    return NDMC_R_AOMAX * (cO_pos / (NDMC_K_OAO + cO_pos)) * NDMC_X_AO * 0.001 / 60.0 / NDMC_M_N
end

function ndmc_conductivity_sink(cO)
    rate = ndmc_rate_term(cO)
    return (NDMC_LA0 - (NDMC_A + NDMC_B * NDMC_LA0) * sqrt(rate)) * rate * 1e3
end

function ndmc_oxygen_sink(cO)
    cO_pos = smooth_nonnegative(cO)
    return (NDMC_R_AOMAX * (cO_pos / (NDMC_K_OAO + cO_pos))) * NDMC_PHI_OAO * NDMC_X_AO / 60.0
end

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
end

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
    )
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

function disturbance_triplet(t_now::Real, cfg::NDMCMPCConfig)
    # Return the current influent conductivity triplet.
    # Outside the disturbance window the plant stays at the baseline feed.
    if cfg.disturbance_start <= t_now < cfg.disturbance_stop
        return cfg.Cin_spike
    end
    return (cfg.Cs, cfg.Cs, cfg.Cs)
end

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

function solve_ndmc_step!(ctrl::TrackingMPCController, sys, cfg::NDMCMPCConfig, state_values::AbstractDict, q_prev::Real, t_now::Real)
    # Update the preview, solve the MPC, and return the first move.
    # This helper is the NDMC-specific wrapper around the generic public
    # function `solve_tracking_mpc!(...)`.
    update_disturbance_forecast!(ctrl, sys, cfg, t_now)
    result = solve_tracking_mpc!(
        ctrl,
        state_values,
        Dict(sys.Q_air => q_prev);
        lower_clip = 0.0,
    )
    q_apply = result.controls[sys.Q_air][1]
    return (q_apply = q_apply, result = result)
end

function metric_column(logctx::MPCLog, key::Symbol)
    return get(logctx.Metrichist, key, fill(NaN, length(logctx.pred_times)))
end

function main()
    mkpath(GENERATED_DIR)

    # Build the plant, the initial state, and the controller.
    cfg = load_config_from_env()
    sys = build_ndmc_system(cfg)
    state0 = ndmc_initial_state(sys, cfg)
    ctrl = build_ndmc_controller(sys; cfg = cfg)

    initial = solve_ndmc_step!(ctrl, sys, cfg, state0, cfg.Q_init, cfg.simulation_span[1])
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

    p0 = Dict(
        sys.Q_air => initial_q,
        sys.Cin1 => initial_cin[1],
        sys.Cin2 => initial_cin[2],
        sys.Cin3 => initial_cin[3],
    )

    prob_cl = ODEProblem(sys, merge(state0, p0), cfg.simulation_span)

    times_ctrl = collect((cfg.simulation_span[1] + cfg.dt):cfg.dt:cfg.simulation_span[2])
    times_log = collect((cfg.simulation_span[1] + cfg.save_dt):cfg.save_dt:cfg.simulation_span[2])
    disturbance_times = unique(filter(t -> cfg.simulation_span[1] <= t <= cfg.simulation_span[2], [cfg.disturbance_start, cfg.disturbance_stop]))

    function mpc_affect!(integ)
        # MPC callback: update disturbance, solve MPC, and apply the first move.
        set_plant_disturbance!(integ, sys, cfg, integ.t)
        state_values = current_state_map(integ, sys)
        q_prev = integ.ps[sys.Q_air]
        solved = solve_ndmc_step!(ctrl, sys, cfg, state_values, q_prev, integ.t)
        q_apply = solved.q_apply
        if abs(q_apply - q_prev) < 1e-8
            q_apply = q_prev
        end
        integ.ps[sys.Q_air] = q_apply
        record_mpc_prediction!(logctx, integ.t, solved.result.predictions)
        record_mpc_metrics!(logctx, solved.result.metrics; missing = :skip)
        return nothing
    end

    function disturbance_affect!(integ)
        # Disturbance-only callback for the plant model.
        set_plant_disturbance!(integ, sys, cfg, integ.t)
        return nothing
    end

    function log_affect!(integ)
        # Logging callback for states and applied control.
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
