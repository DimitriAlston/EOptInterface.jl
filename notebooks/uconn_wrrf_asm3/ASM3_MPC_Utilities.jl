using EOptInterface
using JuMP
using ModelingToolkit
using OrdinaryDiffEq
using DiffEqCallbacks
using Ipopt
using DataFrames
using Printf

Base.@kwdef struct MPCConfig
    PH::Int = 20
    CH::Int = 10
    Δt::Float64 = 0.2
    sp::Float64 = 3.0
    ρ::Float64 = 4e3
    Rsm::Float64 = 5.0
    R1::Float64 = 10.0
    wN::Float64 = 2.0
    α::Float64 = 0.03
    V::Float64 = 1000.0
    KLa_horizon_max::Float64 = 4000.0
    KLa_step_max::Float64 = 400.0
    ΔKLa_max::Float64 = 200.0
    show_status::Bool = true
    status_digits::Int = 3
    status_prefix::String = "[MPC]"
    status_stream::Symbol = :stderr
    verbose::Bool = false
end

mutable struct MPCController
    cfg::MPCConfig
    sys
    model::JuMP.Model
    N::Int
    Nc::Int
    x_vars::Dict{Any, Vector{JuMP.VariableRef}}
    c_ic::Dict{Any, JuMP.ConstraintRef}
    KLa2::Vector{JuMP.VariableRef}
    KLa4::Vector{JuMP.VariableRef}
    c_d1_2::JuMP.ConstraintRef
    c_d1_4::JuMP.ConstraintRef
    term_tr
    term_soft
    term_d
    term_d1
    term_energy
    y_sym
    last_KLa2::Union{Nothing, Vector{Float64}}
    last_KLa4::Union{Nothing, Vector{Float64}}
end

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
end

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

function _mk_ic!(
    model::JuMP.Model,
    v::Vector{JuMP.VariableRef},
    x;
    override::AbstractDict = Dict(),
    rhs0::Float64 = 0.0,
)
    rhs = _rhs_from_override(override, x, rhs0)
    return get_ic_constraint!(model, v; idx = 1, rhs_if_missing = rhs)
end

function _var_bounds(x, x0::Real)
    name = string(x)
    mag = max(abs(float(x0)), 1.0)
    if occursin("flow_rate", name)
        return (0.0, max(20.0, 5.0 * mag))
    elseif occursin("T(t)", name)
        return (0.0, max(20.0, mag + 5.0))
    else
        return (0.0, max(20.0, 10.0 * mag))
    end
end

_base_name_for(x) = replace(string(x), r"[^A-Za-z0-9]+" => "_")

function _shift_start_values!(vars::Vector{JuMP.VariableRef}, prev::Union{Nothing, Vector{Float64}}, fallback::Float64)
    if prev === nothing || isempty(prev)
        for v in vars
            JuMP.set_start_value(v, fallback)
        end
        return nothing
    end
    n = length(vars)
    m = min(length(prev), n)
    if m >= 2
        for k in 1:(m - 1)
            JuMP.set_start_value(vars[k], prev[k + 1])
        end
    end
    JuMP.set_start_value(vars[m], prev[m])
    for k in (m + 1):n
        JuMP.set_start_value(vars[k], prev[end])
    end
    return nothing
end

function _fmt_uconn_scalar(x; digits::Int = 3)
    if x isa Real && isfinite(float(x))
        return @sprintf("%.*f", digits, float(x))
    end
    return "NaN"
end

function _uconn_status_io(cfg::MPCConfig)
    if cfg.status_stream === :stdout
        return stdout
    elseif cfg.status_stream === :stderr
        return stderr
    end
    error("Unsupported UConn MPC status stream $(cfg.status_stream). Use :stdout or :stderr.")
end

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

function build_controller(sys, sol; cfg::MPCConfig = MPCConfig(), y_sym = sys.reactor5.x[4])
    PH, CH, Δt = cfg.PH, cfg.CH, cfg.Δt
    N = PH + 1
    Nc = CH + 1
    tspan_pred = (0.0, PH * Δt)

    decision_vars(sys, [sys.reactor2.KLa, sys.reactor4.KLa])

    model = Model(Ipopt.Optimizer)
    cfg.verbose || set_silent(model)

    x_vars = Dict{Any, Vector{JuMP.VariableRef}}()
    c_ic = Dict{Any, JuMP.ConstraintRef}()

    u0_full = Dict(unknowns(sys) .=> sol[1])
    u0_key_by_str = Dict(string(k) => k for k in keys(u0_full))

    diffvars = Any[]
    for eq in ModelingToolkit.diff_equations(sys)
        var, _ = ModelingToolkit.var_from_nested_derivative(eq.lhs)
        k = get(u0_key_by_str, string(var), nothing)
        k === nothing || push!(diffvars, k)
    end
    diffvar_set = Set(unique(diffvars))

    for x in unknowns(sys)
        x0 = get(u0_full, x, 0.0)
        lb, ub = _var_bounds(x, x0)
        xv = @variable(model, [1:N], lower_bound = lb, upper_bound = ub, base_name = _base_name_for(x))
        x_vars[x] = xv
        if x in diffvar_set
            c_ic[x] = _mk_ic!(model, xv, x; override = u0_full, rhs0 = float(x0))
        end
    end

    @variable(model, 0.0 <= KLa2[1:N] <= cfg.KLa_horizon_max, base_name = "KLa2")
    @variable(model, 0.0 <= KLa4[1:N] <= cfg.KLa_horizon_max, base_name = "KLa4")

    p_disc = [sys.reactor2.KLa, sys.reactor4.KLa]
    p_disc_vars = Dict(sys.reactor2.KLa => KLa2, sys.reactor4.KLa => KLa4)

    iv_JuMP = @variable(model, 0 <= t_stage <= 150)
    JuMP.fix(t_stage, tspan_pred[1]; force = true)

    iv_MTK = ModelingToolkit.get_iv(sys)
    iv_map = Dict(iv_MTK => iv_JuMP)

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

    @variable(model, s_up[1:N] >= 0)
    @variable(model, s_dn[1:N] >= 0)

    @constraint(model, [k = 1:N], x_vars[y_sym][k] <= cfg.sp + s_up[k])
    @constraint(model, [k = 1:N], x_vars[y_sym][k] >= cfg.sp - s_dn[k])

    term_energy = @expression(model, cfg.α * cfg.V * sum(KLa2[k] + KLa4[k] for k = 1:N))
    term_tr = @expression(model, sum((x_vars[y_sym][k] - cfg.sp)^2 for k = 1:N) + cfg.wN * (x_vars[y_sym][N] - cfg.sp)^2)
    term_soft = @expression(model, cfg.ρ * (sum(s_up) + sum(s_dn)))
    term_d = @expression(model, cfg.Rsm * sum(ΔKLa2[k]^2 for k = 2:N) + cfg.Rsm * sum(ΔKLa4[k]^2 for k = 2:N))
    term_d1 = @expression(model, cfg.R1 * d1_2^2 + cfg.R1 * d1_4^2)

    @objective(model, Min, term_tr + term_soft + term_d + term_d1)

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
        term_tr,
        term_soft,
        term_d,
        term_d1,
        term_energy,
        y_sym,
        nothing,
        nothing,
    ), u0_full
end

function _mpc_affect!(integ, ctrl::MPCController, i_state::AbstractDict, logctx::UConnMPCLog)
    for (var, cref) in ctrl.c_ic
        idx = get(i_state, var, nothing)
        idx === nothing && continue
        JuMP.set_normalized_rhs(cref, float(integ.u[idx]))
    end

    K2_prev = float(integ.ps[ctrl.sys.reactor2.KLa])
    K4_prev = float(integ.ps[ctrl.sys.reactor4.KLa])
    JuMP.set_normalized_rhs(ctrl.c_d1_2, K2_prev)
    JuMP.set_normalized_rhs(ctrl.c_d1_4, K4_prev)

    JuMP.set_lower_bound(ctrl.KLa2[1], max(0.0, K2_prev - ctrl.cfg.ΔKLa_max))
    JuMP.set_upper_bound(ctrl.KLa2[1], min(ctrl.cfg.KLa_step_max, K2_prev + ctrl.cfg.ΔKLa_max))
    JuMP.set_lower_bound(ctrl.KLa4[1], max(0.0, K4_prev - ctrl.cfg.ΔKLa_max))
    JuMP.set_upper_bound(ctrl.KLa4[1], min(ctrl.cfg.KLa_step_max, K4_prev + ctrl.cfg.ΔKLa_max))

    _shift_start_values!(ctrl.KLa2, ctrl.last_KLa2, K2_prev)
    _shift_start_values!(ctrl.KLa4, ctrl.last_KLa4, K4_prev)

    optimize!(ctrl.model)
    status = JuMP.termination_status(ctrl.model)
    status_io = _uconn_status_io(ctrl.cfg)
    push!(logctx.solve_times, float(integ.t))
    push!(logctx.statuses, string(status))

    if JuMP.has_values(ctrl.model)
        ctrl.last_KLa2 = value.(ctrl.KLa2)
        ctrl.last_KLa4 = value.(ctrl.KLa4)
        push!(logctx.objectives, JuMP.objective_value(ctrl.model))

        KLa2_now = value(ctrl.KLa2[1])
        KLa4_now = value(ctrl.KLa4[1])
        if abs(KLa2_now - K2_prev) < 1e-3
            KLa2_now = K2_prev
        end
        if abs(KLa4_now - K4_prev) < 1e-3
            KLa4_now = K4_prev
        end
        integ.ps[ctrl.sys.reactor2.KLa] = clamp(KLa2_now, 0.0, ctrl.cfg.KLa_step_max)
        integ.ps[ctrl.sys.reactor4.KLa] = clamp(KLa4_now, 0.0, ctrl.cfg.KLa_step_max)
        if ctrl.cfg.show_status
            print_uconn_mpc_status(
                status_io,
                ctrl,
                status,
                JuMP.objective_value(ctrl.model),
                integ.t,
                integ.u[i_state[ctrl.y_sym]],
                ctrl.cfg.sp,
                integ.ps[ctrl.sys.reactor2.KLa],
                integ.ps[ctrl.sys.reactor4.KLa],
                digits = ctrl.cfg.status_digits,
                prefix = ctrl.cfg.status_prefix,
            )
        end
    else
        push!(logctx.objectives, NaN)
        integ.ps[ctrl.sys.reactor2.KLa] = K2_prev
        integ.ps[ctrl.sys.reactor4.KLa] = K4_prev
        if ctrl.cfg.show_status
            print_uconn_mpc_status(
                status_io,
                ctrl,
                status,
                NaN,
                integ.t,
                integ.u[i_state[ctrl.y_sym]],
                ctrl.cfg.sp,
                integ.ps[ctrl.sys.reactor2.KLa],
                integ.ps[ctrl.sys.reactor4.KLa],
                digits = ctrl.cfg.status_digits,
                prefix = ctrl.cfg.status_prefix,
            )
        end
    end
    return nothing
end

function run_uconn_closed_loop_mpc(
    sys,
    sol;
    y_sym = sys.reactor5.x[4],
    simulation_span = (0.0, 5.0),
    simulation_span_long = (0.0, 10.0),
    cfg_hi::MPCConfig = MPCConfig(sp = 3.0),
    cfg_lo::MPCConfig = MPCConfig(sp = 2.0),
)
    ctrl_hi, u0_full = build_controller(sys, sol; cfg = cfg_hi, y_sym = y_sym)
    ctrl_lo, _ = build_controller(sys, sol; cfg = cfg_lo, y_sym = y_sym)

    set_attribute(ctrl_hi.model, "max_cpu_time", 20.0)
    set_attribute(ctrl_hi.model, "acceptable_tol", 1e-4)
    set_attribute(ctrl_lo.model, "max_cpu_time", 20.0)
    set_attribute(ctrl_lo.model, "acceptable_tol", 1e-4)

    Δt_ctrl = cfg_hi.Δt
    times_ctrl = collect(first(simulation_span):Δt_ctrl:last(simulation_span))
    times_ctrl2 = collect(last(simulation_span):Δt_ctrl:last(simulation_span_long))
    savetimes = collect(first(simulation_span_long):Δt_ctrl / 10:last(simulation_span_long))

    unk = unknowns(sys)
    i_state = Dict(var => i for (i, var) in pairs(unk))
    logctx = UConnMPCLog(i_state = i_state)

    setpoint_fn(tnow) = tnow < last(simulation_span) ? cfg_hi.sp : cfg_lo.sp

    function _log_affect!(integ)
        push!(logctx.ts, float(integ.t))
        push!(logctx.setpoints, float(setpoint_fn(integ.t)))
        push!(logctx.SNH4, float(integ.u[i_state[y_sym]]))
        push!(logctx.KLa2, float(integ.ps[sys.reactor2.KLa]))
        push!(logctx.KLa4, float(integ.ps[sys.reactor4.KLa]))
        return nothing
    end

    function _mpc_affect_hi!(integ)
        _mpc_affect!(integ, ctrl_hi, i_state, logctx)
    end

    function _mpc_affect_lo!(integ)
        _mpc_affect!(integ, ctrl_lo, i_state, logctx)
    end

    log_cb = PresetTimeCallback(savetimes, _log_affect!; save_positions = (false, false))
    mpc_cb = PresetTimeCallback(times_ctrl, _mpc_affect_hi!; save_positions = (false, false))
    mpc_cb2 = PresetTimeCallback(times_ctrl2, _mpc_affect_lo!; save_positions = (false, false))
    callbacks = CallbackSet(log_cb, mpc_cb, mpc_cb2)

    prob_cl = ODEProblem(sys, u0_full, simulation_span_long; warn_initialize_determined = false)
    sol_cl = solve(
        prob_cl,
        FBDF();
        adaptive = true,
        initializealg = OrdinaryDiffEqNonlinearSolve.BrownFullBasicInit(),
        callback = callbacks,
        tstops = union(savetimes, times_ctrl, times_ctrl2),
        saveat = savetimes,
    )

    log_df = DataFrame(
        time = logctx.ts,
        setpoint = logctx.setpoints,
        SNH4 = logctx.SNH4,
        KLa2 = logctx.KLa2,
        KLa4 = logctx.KLa4,
    )

    return (
        ctrl_hi = ctrl_hi,
        ctrl_lo = ctrl_lo,
        sol = sol_cl,
        log = logctx,
        log_df = log_df,
    )
end
