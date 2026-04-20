using EOptInterface
using JuMP
using ModelingToolkit
using OrdinaryDiffEq
using DiffEqCallbacks
using Ipopt
using DataFrames
using Printf

"""
    MPCConfig

Configuration bundle for the UConn ASM3 MPC example.

This struct keeps the public notebook settings in one place: horizons,
tracking weights, actuator limits, and live status-printing behavior.
"""
Base.@kwdef struct MPCConfig
    PH::Int = 20
    CH::Int = 10
    Δt::Float64 = 0.2

    sp::Float64 = 3.0
    ρ::Float64 = 4e3
    Rsm::Float64 = 5.0
    R1::Float64 = 10.0
    wN::Float64 = 2.0

    w_energy::Float64 = 0.0
    α::Float64 = 0.03
    V::Float64 = 1000.0

    KLa_min::Float64 = 0.0
    KLa_max::Float64 = 400.0
    ΔKLa_max::Float64 = 200.0

    show_status::Bool = true
    status_digits::Int = 3
    status_prefix::String = "[MPC]"
    status_stream::Symbol = :stderr
    verbose::Bool = false
end

"""
    MPCController

Container for the reusable ASM3 MPC problem.

It stores the built JuMP model, references to the state and control
trajectories, the objective terms used for logging, and the most recent solved
trajectory used for warm starts.
"""
mutable struct MPCController
    cfg::MPCConfig
    sys
    model::JuMP.Model
    N::Int
    Nc::Int
    x_vars::Dict{Num, Vector{JuMP.VariableRef}}
    c_ic::Dict{Any, JuMP.ConstraintRef}
    KLa2::Vector{JuMP.VariableRef}
    KLa4::Vector{JuMP.VariableRef}
    c_d1_2::JuMP.ConstraintRef
    c_d1_4::JuMP.ConstraintRef
    sp_param::JuMP.VariableRef
    y_sym::Num
    term_tr
    term_soft
    term_d
    term_d1
    term_energy
    last_KLa2::Union{Nothing, Vector{Float64}}
    last_KLa4::Union{Nothing, Vector{Float64}}
    last_state_trajs::Union{Nothing, Dict{Num, Vector{Float64}}}
end

"""
    UConnMPCLog

Simple log container for the closed-loop ASM3 notebook.

Each vector records one quantity over time so the notebook can build plots and
tables after the simulation finishes.
"""
Base.@kwdef mutable struct UConnMPCLog
    i_state::AbstractDict
    ts::Vector{Float64} = Float64[]
    setpoints::Vector{Float64} = Float64[]
    SNH4::Vector{Float64} = Float64[]
    KLa2::Vector{Float64} = Float64[]
    KLa4::Vector{Float64} = Float64[]
    solve_times::Vector{Float64} = Float64[]
    statuses::Vector{String} = String[]
    objectives::Vector{Float64} = Float64[]
    track_terms::Vector{Float64} = Float64[]
    soft_terms::Vector{Float64} = Float64[]
    move_terms::Vector{Float64} = Float64[]
    first_move_terms::Vector{Float64} = Float64[]
    energy_terms::Vector{Float64} = Float64[]
end

"""
    _rhs_from_override(override, x, rhs0)

Read the initial value to use for symbolic variable `x`.

The helper accepts both direct symbolic keys and string-equivalent keys. This
is useful when state dictionaries come from slightly different MTK code paths.
"""
function _rhs_from_override(override, x, rhs0)
    if haskey(override, x)
        return float(override[x])
    end
    sx = string(x)
    for (k, v) in override
        if string(k) == sx
            return float(v)
        end
    end
    return rhs0
end

"""
    _mk_ic!(model, v, x; override=Dict(), rhs0=0.0)

Create or recover the initial-condition constraint for trajectory `v`.

This small wrapper keeps the controller builder readable by hiding the
"find the right initial value, then ensure the IC constraint exists" logic.
"""
function _mk_ic!(
    model::JuMP.Model,
    v::Vector{JuMP.VariableRef},
    x::Num;
    override::AbstractDict = Dict(),
    rhs0::Float64 = 0.0,
)
    rhs = _rhs_from_override(override, x, rhs0)
    return get_ic_constraint!(model, v; idx = 1, rhs_if_missing = rhs)
end

"""
    _fmt_uconn_scalar(x; digits=3)

Format one scalar for live MPC status printing.

Finite numeric values are printed with a fixed number of digits. Everything
else falls back to `"NaN"` so the log output remains readable.
"""
function _fmt_uconn_scalar(x; digits::Int = 3)
    if x isa Real && isfinite(float(x))
        return @sprintf("%.*f", digits, float(x))
    end
    return "NaN"
end

"""
    _uconn_status_io(cfg)

Return the output stream used by live MPC status printing.

The notebook lets users choose `:stdout` or `:stderr` so the status lines can
either live in the normal cell output or in the error stream.
"""
function _uconn_status_io(cfg::MPCConfig)
    if cfg.status_stream === :stdout
        return stdout
    elseif cfg.status_stream === :stderr
        return stderr
    end
    error("Unsupported UConn MPC status stream $(cfg.status_stream). Use :stdout or :stderr.")
end

"""
    print_uconn_mpc_status(io, ctrl, status, objective, tnow, y_now, setpoint, KLa2_now, KLa4_now; ...)

Print one human-readable summary line after an MPC trigger.

This is the notebook-facing status logger. It is intentionally compact so users
can watch the controller solve step by step while the closed-loop simulation
runs.
"""
function print_uconn_mpc_status(
    io::IO,
    ctrl::MPCController,
    status,
    objective,
    tnow,
    y_now,
    setpoint,
    KLa2_now,
    KLa4_now;
    digits::Int = 3,
    prefix::AbstractString = "[MPC]",
)
    accepted = is_accepted_mpc_status(status; accepted_statuses = default_mpc_accepted_statuses())
    parts = String[
        prefix,
        "t=" * _fmt_uconn_scalar(tnow; digits = digits),
        "status=" * string(status),
        "accepted=" * string(accepted),
        "obj=" * _fmt_uconn_scalar(objective; digits = digits),
        "SNH4=" * _fmt_uconn_scalar(y_now; digits = digits) * " (sp=" * _fmt_uconn_scalar(setpoint; digits = digits) * ")",
        "KLa2=" * _fmt_uconn_scalar(KLa2_now; digits = digits),
        "KLa4=" * _fmt_uconn_scalar(KLa4_now; digits = digits),
    ]
    println(io, string(first(parts), " ", join(parts[2:end], " | ")))
    flush(io)
    return nothing
end

"""
    build_controller(sys, sol; cfg=MPCConfig(), y_sym=sys.reactor5.x[4])

Build the reusable JuMP-based MPC controller for the UConn ASM3 flowsheet.

This is the expensive one-time setup step. The function:
1. allocates state and control trajectories,
2. registers the ODE prediction model,
3. adds move-rate and tracking constraints,
4. builds the objective terms used by the notebook logs.
"""
function build_controller(sys, sol; cfg::MPCConfig = MPCConfig(), y_sym = sys.reactor5.x[4])
    PH, CH, Δt = cfg.PH, cfg.CH, cfg.Δt
    N = PH + 1
    Nc = CH + 1

    decision_vars(sys, [sys.reactor2.KLa, sys.reactor4.KLa])

    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "tol", 1e-6)
    set_optimizer_attribute(model, "max_cpu_time", 300.0)
    set_optimizer_attribute(model, "print_level", 0)
    cfg.verbose || set_silent(model)

    # The notebook controller only keeps the five reactor state blocks. This is
    # the smallest state set that still matches the latest working local ASM3
    # MPC setup.
    reactor_syms = (:reactor1, :reactor2, :reactor3, :reactor4, :reactor5)
    x_lb = fill(0.0, 13)
    x_ub = [20.0; fill(1.0e6, 12)...]
    x_vars = Dict{Num, Vector{JuMP.VariableRef}}()
    c_ic = Dict{Any, JuMP.ConstraintRef}()

    u0_full = Dict(unknowns(sys) .=> sol[1])
    u0_key_by_str = Dict(string(k) => k for k in keys(u0_full))
    diffvars = Any[]
    for eq in ModelingToolkit.diff_equations(sys)
        var, _ = ModelingToolkit.var_from_nested_derivative(eq.lhs)
        k = get(u0_key_by_str, string(var), nothing)
        k === nothing || push!(diffvars, k)
    end
    u0_dict = Dict(k => u0_full[k] for k in unique(diffvars))
    isempty(u0_dict) && (u0_dict = u0_full)

    for rname in reactor_syms
        r = getproperty(sys, rname)
        for i in 1:13
            # Each reactor state becomes one trajectory over the prediction
            # horizon. These are the core decision trajectories of the MPC.
            traj = @variable(
                model,
                [1:N],
                lower_bound = x_lb[i],
                upper_bound = x_ub[i],
                base_name = "$(rname)_x$(i)",
            )
            x_vars[r.x[i]] = traj
            c_ic[r.x[i]] = _mk_ic!(model, traj, r.x[i]; override = u0_dict, rhs0 = 0.0)
        end
    end

    # A few algebraic flow variables still need trajectories because the
    # structurally simplified flowsheet keeps them in the ODE residuals. This
    # smaller set is what avoids the old time-limit behavior.
    aux_flow_rhs = Dict(
        sys.splitter1.In.flow_rate => Ini2vecflow,
        sys.mixer1.Out1.flow_rate => 1.5 * Ini2vecflow,
        sys.mixer3.Out1.flow_rate => 1.75 * Ini2vecflow,
        sys.clarifier.inlet_stream.flow_rate => Ini2vecflow,
    )
    for (var_mtk, rhs0) in aux_flow_rhs
        haskey(u0_key_by_str, string(var_mtk)) || continue
        traj = @variable(
            model,
            [1:N],
            lower_bound = 0.0,
            upper_bound = 1.0e6,
            base_name = sanitize_mpc_name(var_mtk),
        )
        x_vars[var_mtk] = traj
        c_ic[var_mtk] = _mk_ic!(model, traj, var_mtk; override = u0_full, rhs0 = float(rhs0))
    end

    # The manipulated variables are the aeration coefficients in reactors 2 and
    # 4. They are allowed to move freely only over the control horizon.
    @variable(model, cfg.KLa_min <= KLa2[1:N] <= cfg.KLa_max, base_name = "KLa2")
    @variable(model, cfg.KLa_min <= KLa4[1:N] <= cfg.KLa_max, base_name = "KLa4")

    p_disc = [sys.reactor2.KLa, sys.reactor4.KLa]
    p_disc_vars = Dict(sys.reactor2.KLa => KLa2, sys.reactor4.KLa => KLa4)

    iv_JuMP = @variable(model, 0 <= t_stage <= 150)
    JuMP.fix(t_stage, 0.0; force = true)
    iv_MTK = ModelingToolkit.get_iv(sys)
    iv_map = Dict(iv_MTK => iv_JuMP)

    # Register the prediction model once. After this, each MPC trigger only
    # updates ICs, setpoints, and warm-start guesses before re-solving.
    register_odesystem(
        model,
        sys,
        (0.0, PH * Δt),
        Δt,
        "IE";
        x_vars = x_vars,
        p_disc = p_disc,
        p_disc_vars = p_disc_vars,
        t_map = iv_map,
    )

    if Nc < N
        # After the control horizon ends, hold the last move constant.
        for k in (Nc + 1):N
            @constraint(model, KLa2[k] == KLa2[Nc])
            @constraint(model, KLa4[k] == KLa4[Nc])
        end
    end

    @variable(model, d1_2)
    @variable(model, d1_4)
    c_d1_2 = @constraint(model, KLa2[1] - d1_2 == 0.0)
    c_d1_4 = @constraint(model, KLa4[1] - d1_4 == 0.0)

    @variable(model, ΔKLa2[2:N])
    @constraint(model, [k = 2:N], ΔKLa2[k] == KLa2[k] - KLa2[k - 1])
    @constraint(model, [k = 2:N], -cfg.ΔKLa_max <= ΔKLa2[k] <= cfg.ΔKLa_max)

    @variable(model, ΔKLa4[2:N])
    @constraint(model, [k = 2:N], ΔKLa4[k] == KLa4[k] - KLa4[k - 1])
    @constraint(model, [k = 2:N], -cfg.ΔKLa_max <= ΔKLa4[k] <= cfg.ΔKLa_max)

    @variable(model, sp_param)
    JuMP.fix(sp_param, cfg.sp; force = true)

    @variable(model, s_up[1:N] >= 0)
    @variable(model, s_dn[1:N] >= 0)

    # These soft constraints let the controller violate the setpoint band if it
    # must, but every violation is penalized in the objective.
    @constraint(model, [k = 1:N], x_vars[y_sym][k] <= sp_param + s_up[k])
    @constraint(model, [k = 1:N], x_vars[y_sym][k] >= sp_param - s_dn[k])

    term_tr = @expression(
        model,
        sum((x_vars[y_sym][k] - sp_param)^2 for k = 1:N) +
        cfg.wN * (x_vars[y_sym][N] - sp_param)^2,
    )
    term_soft = @expression(model, cfg.ρ * (sum(s_up) + sum(s_dn)))
    term_d = @expression(
        model,
        cfg.Rsm * (sum(ΔKLa2[k]^2 for k = 2:N) + sum(ΔKLa4[k]^2 for k = 2:N)),
    )
    term_d1 = @expression(model, cfg.R1 * (d1_2^2 + d1_4^2))
    term_energy = @expression(
        model,
        (24 * cfg.Δt) * cfg.α * cfg.V * sum(KLa2[k] + KLa4[k] for k = 1:N),
    )

    # The full objective is stored in named pieces so the notebook can report
    # how much each modeling choice contributes at every solve.
    @objective(model, Min, term_tr + term_soft + term_d + term_d1 + cfg.w_energy * term_energy)

    return MPCController(
        cfg,
        sys,
        model,
        N,
        Nc,
        x_vars,
        c_ic,
        KLa2,
        KLa4,
        c_d1_2,
        c_d1_4,
        sp_param,
        y_sym,
        term_tr,
        term_soft,
        term_d,
        term_d1,
        term_energy,
        nothing,
        nothing,
        nothing,
    ), u0_dict
end

"""
    warm_start_trajectories!(ctrl)

Shift the previously solved trajectories forward by one stage.

Warm starts are important in MPC because consecutive optimization problems are
similar. Reusing the last solution often reduces solve time significantly.
"""
function warm_start_trajectories!(ctrl::MPCController)
    if ctrl.last_state_trajs !== nothing
        for (var, traj) in ctrl.x_vars
            old = get(ctrl.last_state_trajs, var, nothing)
            old === nothing && continue
            for k in 1:(ctrl.N - 1)
                JuMP.set_start_value(traj[k], old[k + 1])
            end
            JuMP.set_start_value(traj[ctrl.N], old[end])
        end
    end

    if ctrl.last_KLa2 === nothing || ctrl.last_KLa4 === nothing
        return nothing
    end
    for k in 1:(ctrl.N - 1)
        JuMP.set_start_value(ctrl.KLa2[k], ctrl.last_KLa2[k + 1])
        JuMP.set_start_value(ctrl.KLa4[k], ctrl.last_KLa4[k + 1])
    end
    JuMP.set_start_value(ctrl.KLa2[ctrl.N], ctrl.last_KLa2[end])
    JuMP.set_start_value(ctrl.KLa4[ctrl.N], ctrl.last_KLa4[end])
    return nothing
end

"""
    setpoint!(ctrl, sp_new)

Update the tracked setpoint of the controller in place.

The controller model is reused, so changing the setpoint simply means updating
the fixed JuMP parameter used in the objective and soft constraints.
"""
function setpoint!(ctrl::MPCController, sp_new::Real)
    JuMP.fix(ctrl.sp_param, float(sp_new); force = true)
    return nothing
end

"""
    make_logctx(sys)

Create the log container used by the closed-loop notebook example.

The returned object also stores a fast state-index map so callbacks can extract
the tracked state from the ODE integrator without recomputing symbol lookups.
"""
function make_logctx(sys)
    unk = unknowns(sys)
    return UConnMPCLog(i_state = Dict(var => i for (i, var) in pairs(unk)))
end

"""
    _log_current_state!(logctx, ctrl, integ)

Append the current plant state and applied controls to the log.

This helper is called from the logging callback so the notebook can later plot
the closed-loop trajectories.
"""
function _log_current_state!(logctx::UConnMPCLog, ctrl::MPCController, integ)
    y_idx = logctx.i_state[ctrl.y_sym]
    push!(logctx.ts, float(integ.t))
    push!(logctx.setpoints, float(JuMP.fix_value(ctrl.sp_param)))
    push!(logctx.SNH4, float(integ.u[y_idx]))
    push!(logctx.KLa2, float(integ.ps[ctrl.sys.reactor2.KLa]))
    push!(logctx.KLa4, float(integ.ps[ctrl.sys.reactor4.KLa]))
    return nothing
end

"""
    mpc_solve_step!(ctrl, integ, logctx=nothing)

Run one online MPC update.

The function performs the standard receding-horizon sequence:
1. copy the current plant state into the IC constraints,
2. enforce first-move consistency with the previously applied input,
3. warm start the optimization variables,
4. solve the JuMP model,
5. apply the first move if the solve status is acceptable,
6. record metrics and print a status line.
"""
function mpc_solve_step!(ctrl::MPCController, integ, logctx::Union{Nothing, UConnMPCLog}=nothing)
    cfg = ctrl.cfg
    sys = ctrl.sys
    i_state = logctx === nothing ? Dict(var => i for (i, var) in pairs(unknowns(sys))) : logctx.i_state

    # Refresh all state IC constraints from the latest simulated plant state.
    for (var, cref) in ctrl.c_ic
        JuMP.set_normalized_rhs(cref, float(integ.u[i_state[var]]))
    end

    K2_prev = float(integ.ps[sys.reactor2.KLa])
    K4_prev = float(integ.ps[sys.reactor4.KLa])
    JuMP.set_normalized_rhs(ctrl.c_d1_2, K2_prev)
    JuMP.set_normalized_rhs(ctrl.c_d1_4, K4_prev)

    # Tighten the first move so it cannot jump further than the allowed
    # move-rate limit from the previously applied value.
    JuMP.set_lower_bound(ctrl.KLa2[1], max(cfg.KLa_min, K2_prev - cfg.ΔKLa_max))
    JuMP.set_upper_bound(ctrl.KLa2[1], min(cfg.KLa_max, K2_prev + cfg.ΔKLa_max))
    JuMP.set_lower_bound(ctrl.KLa4[1], max(cfg.KLa_min, K4_prev - cfg.ΔKLa_max))
    JuMP.set_upper_bound(ctrl.KLa4[1], min(cfg.KLa_max, K4_prev + cfg.ΔKLa_max))

    warm_start_trajectories!(ctrl)
    optimize!(ctrl.model)

    st = JuMP.termination_status(ctrl.model)
    accepted = is_accepted_mpc_status(st; accepted_statuses = default_mpc_accepted_statuses())

    info = if accepted && JuMP.has_values(ctrl.model)
        # Apply only the first optimized move. The rest of the optimized
        # trajectory is used only as a prediction and warm start.
        K2 = float(JuMP.value(ctrl.KLa2[1]))
        K4 = float(JuMP.value(ctrl.KLa4[1]))
        abs(K2 - K2_prev) < 1e-3 && (K2 = K2_prev)
        abs(K4 - K4_prev) < 1e-3 && (K4 = K4_prev)
        integ.ps[sys.reactor2.KLa] = K2
        integ.ps[sys.reactor4.KLa] = K4
        ctrl.last_KLa2 = collect(float.(JuMP.value.(ctrl.KLa2)))
        ctrl.last_KLa4 = collect(float.(JuMP.value.(ctrl.KLa4)))
        ctrl.last_state_trajs = Dict(
            var => collect(float.(JuMP.value.(traj))) for (var, traj) in ctrl.x_vars
        )
        (
            status = st,
            obj = float(JuMP.objective_value(ctrl.model)),
            tr = float(JuMP.value(ctrl.term_tr)),
            soft = float(JuMP.value(ctrl.term_soft)),
            d = float(JuMP.value(ctrl.term_d)),
            d1 = float(JuMP.value(ctrl.term_d1)),
            energy = float(JuMP.value(ctrl.term_energy)),
        )
    else
        integ.ps[sys.reactor2.KLa] = K2_prev
        integ.ps[sys.reactor4.KLa] = K4_prev
        (
            status = st,
            obj = NaN,
            tr = NaN,
            soft = NaN,
            d = NaN,
            d1 = NaN,
            energy = NaN,
        )
    end

    if logctx !== nothing
        push!(logctx.solve_times, float(integ.t))
        push!(logctx.statuses, string(info.status))
        push!(logctx.objectives, info.obj)
        push!(logctx.track_terms, info.tr)
        push!(logctx.soft_terms, info.soft)
        push!(logctx.move_terms, info.d)
        push!(logctx.first_move_terms, info.d1)
        push!(logctx.energy_terms, info.energy)
    end

    if cfg.show_status
        print_uconn_mpc_status(
            _uconn_status_io(cfg),
            ctrl,
            info.status,
            info.obj,
            integ.t,
            integ.u[i_state[ctrl.y_sym]],
            JuMP.fix_value(ctrl.sp_param),
            integ.ps[sys.reactor2.KLa],
            integ.ps[sys.reactor4.KLa];
            digits = cfg.status_digits,
            prefix = cfg.status_prefix,
        )
    end

    return integ.ps[sys.reactor2.KLa], integ.ps[sys.reactor4.KLa], info
end

"""
    run_uconn_closed_loop_mpc(sys, sol; ...)

Run the full closed-loop UConn ASM3 MPC demonstration.

This helper wires together:
1. the reusable controller,
2. the ODE plant simulation,
3. callbacks for logging, MPC triggers, and setpoint changes,
4. the final DataFrame used by the notebook plots.
"""
function run_uconn_closed_loop_mpc(
    sys,
    sol;
    y_sym = sys.reactor5.x[4],
    simulation_span = (0.0, 5.0),
    simulation_span_long = (0.0, 10.0),
    cfg_hi::MPCConfig = MPCConfig(sp = 3.0),
    cfg_lo::MPCConfig = MPCConfig(sp = 2.0),
)
    ctrl, u0_dict = build_controller(sys, sol; cfg = cfg_hi, y_sym = y_sym)
    logctx = make_logctx(sys)

    Δt_ctrl = cfg_hi.Δt
    t0, tf = simulation_span_long
    switch_time = last(simulation_span)
    n_ctrl = max(0, Int(floor((tf - t0) / Δt_ctrl + 1e-9)))
    mpc_times = [t0 + k * Δt_ctrl for k in 0:(n_ctrl - 1)]
    log_times = collect(t0:Δt_ctrl:tf)
    tstops = sort(unique(vcat(log_times, mpc_times, [switch_time])))

    function _log_affect!(integ)
        _log_current_state!(logctx, ctrl, integ)
        return nothing
    end

    function _mpc_affect!(integ)
        mpc_solve_step!(ctrl, integ, logctx)
        return nothing
    end

    function _sp_change_affect!(integ)
        setpoint!(ctrl, cfg_lo.sp)
        return nothing
    end

    callbacks = Any[
        PresetTimeCallback(log_times, _log_affect!; save_positions = (false, false)),
        PresetTimeCallback(mpc_times, _mpc_affect!; save_positions = (false, false)),
    ]
    if switch_time > t0 && switch_time < tf
        push!(callbacks, PresetTimeCallback([switch_time], _sp_change_affect!; save_positions = (false, false)))
    end

    cl_guesses = Dict(
        sys.splitter1.In.flow_rate => Ini2vecflow,
        sys.mixer1.Out1.flow_rate => 1.5 * Ini2vecflow,
        sys.mixer3.Out1.flow_rate => 1.75 * Ini2vecflow,
        sys.clarifier.inlet_stream.flow_rate => Ini2vecflow,
        sys.clarifier.recycle_stream.flow_rate => 0.4 * Ini2vecflow,
    )
    # The plant simulation still needs a few flow guesses for stable
    # initialization of the full ASM3 flowsheet.
    prob_cl = ODEProblem(
        sys,
        u0_dict,
        simulation_span_long;
        guesses = cl_guesses,
        warn_initialize_determined = false,
    )
    sol_cl = solve(
        prob_cl,
        FBDF();
        adaptive = true,
        initializealg = OrdinaryDiffEqNonlinearSolve.BrownFullBasicInit(),
        callback = CallbackSet(callbacks...),
        tstops = tstops,
        saveat = log_times,
    )

    log_df = DataFrame(
        time = logctx.ts,
        setpoint = logctx.setpoints,
        SNH4 = logctx.SNH4,
        KLa2 = logctx.KLa2,
        KLa4 = logctx.KLa4,
    )

    return (
        ctrl = ctrl,
        sol = sol_cl,
        log = logctx,
        log_df = log_df,
    )
end
