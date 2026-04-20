#
# DMC utilities
#
# This file implements the lighter DMC path.
# Unlike the tracking MPC code, it starts from a numeric step response, not a
# full MTK ODE or DAE system.

"""
    DMCRegistration

Container returned by `register_dmcsystem(...)`.

It stores the step-response data, the JuMP variables, and the prediction
constraints.

Think of this object as the DMC counterpart of `TrackingMPCController`.
It is the object that keeps all DMC-related JuMP references together so later
code can update data, inspect predictions, or build an objective on top of the
prediction model.
"""
struct DMCRegistration
    s::Matrix{Float64}
    N::Int
    P::Int
    M::Int
    Ny::Int
    u::Vector{JuMP.VariableRef}
    y_pred::Matrix{JuMP.VariableRef}
    y_meas::Vector{JuMP.VariableRef}
    u_hist::Vector{JuMP.VariableRef}
    du::Vector{JuMP.VariableRef}
    move_constraints::Vector{JuMP.ConstraintRef}
    prediction_constraints::Matrix{JuMP.ConstraintRef}
end

function _step_response_matrix(s; output_count::Union{Nothing, Int}=nothing)
    if s isa AbstractVector
        # A vector means one output. Convert it to a column matrix.
        mat = reshape(Float64.(collect(s)), :, 1)
        if !isnothing(output_count)
            output_count >= 1 || error("`output_count` must be at least 1.")
            return repeat(mat, 1, output_count)
        end
        return mat
    elseif s isa AbstractMatrix
        !isnothing(output_count) && size(s, 2) != output_count &&
            error("`output_count=$(output_count)` does not match size(s, 2) = $(size(s, 2)).")
        return Float64.(Matrix(s))
    else
        error("Step response `s` must be a vector or matrix.")
    end
end

function _float_vector(data, len::Int, label::AbstractString)
    values = Float64.(collect(data))
    length(values) == len || error("`$(label)` must have length $(len), got $(length(values)).")
    return values
end

function _hold_u(ctx::DMCRegistration, idx::Int)
    hist_len = ctx.N - 1
    if hist_len > 0 && idx <= hist_len
        # Early indices use the fixed input history.
        return ctx.u_hist[idx]
    elseif idx <= hist_len + ctx.M
        # Middle indices use optimized future moves.
        return ctx.u[idx - hist_len]
    else
        # After the control horizon, hold the last move.
        return ctx.u[ctx.M]
    end
end

"""
    register_dmcsystem(model, s, P; kwargs...)

Register linear DMC prediction constraints in a JuMP model using a numeric step
response `s`.

Arguments:
- `model`: JuMP model
- `s`: step response coefficients. A vector means one output. A matrix means
  `size(s) = (N, Ny)` where `N` is the DMC model horizon and `Ny` is the number
  of predicted outputs.
- `P`: prediction horizon length in samples

Keyword arguments:
- `M`: control horizon length in samples. Defaults to `P`.
- `y_meas0`: current measured output vector. Defaults to zeros.
- `u_hist0`: past absolute manipulated-input history of length `N - 1`, ordered
  from oldest to newest. Defaults to zeros.
- `u_bounds`: lower and upper bound tuple for the future manipulated input.
- `output_count`: optional output count. If `s` is a vector and `output_count > 1`,
  the same step response is reused for each output column.
- `base_name`: prefix for created JuMP variables.
- `store_ext`: if true, store the registration in `model.ext[:dmc]`.

Return a `DMCRegistration`.

Algorithm:
1. Convert the step response `s` into a dense matrix with one column per output.
2. Check the requested prediction horizon `P` and control horizon `M`.
3. Create JuMP variables for future manipulated inputs and predicted outputs.
4. Fix the current measured outputs and past input history to the supplied
   values.
5. Add move-difference variables when `M > 1`.
6. For each prediction step, build the DMC output equation from the stored step
   response coefficients.
7. Return a registration object that keeps the variables, dimensions, and
   update helpers together.
"""
function register_dmcsystem(model::JuMP.Model,
                            s,
                            P::Integer;
                            M::Union{Nothing, Integer}=nothing,
                            y_meas0=nothing,
                            u_hist0=nothing,
                            u_bounds::Tuple{Real, Real}=(-Inf, Inf),
                            output_count::Union{Nothing, Integer}=nothing,
                            base_name::AbstractString="dmc",
                            store_ext::Bool=true)
    step_mat = _step_response_matrix(s; output_count=output_count)
    N, Ny = size(step_mat)
    P = Int(P)
    P >= 1 || error("`P` must be at least 1.")
    M = isnothing(M) ? P : Int(M)
    1 <= M <= P || error("`M` must satisfy 1 <= M <= P.")

    y_meas_init = isnothing(y_meas0) ? zeros(Float64, Ny) : _float_vector(y_meas0, Ny, "y_meas0")
    u_hist_len = max(N - 1, 0)
    u_hist_init = isnothing(u_hist0) ? zeros(Float64, u_hist_len) : _float_vector(u_hist0, u_hist_len, "u_hist0")

    # Build the future control moves and predicted outputs.
    @variable(model, u_bounds[1] <= u[1:M] <= u_bounds[2], base_name="$(base_name)_u")
    @variable(model, y_pred[1:P, 1:Ny], base_name="$(base_name)_y")
    @variable(model, y_meas[1:Ny], base_name="$(base_name)_y_meas")
    for l in 1:Ny
        JuMP.fix(y_meas[l], y_meas_init[l]; force=true)
    end

    if u_hist_len > 0
        @variable(model, u_hist[1:u_hist_len], base_name="$(base_name)_u_hist")
        for i in 1:u_hist_len
            JuMP.fix(u_hist[i], u_hist_init[i]; force=true)
        end
    else
        u_hist = JuMP.VariableRef[]
    end

    if M > 1
        @variable(model, du[1:(M - 1)], base_name="$(base_name)_du")
        move_constraints = [
            @constraint(model, du[k] == u[k + 1] - u[k]) for k in 1:(M - 1)
        ]
    else
        du = JuMP.VariableRef[]
        move_constraints = JuMP.ConstraintRef[]
    end

    ctx = DMCRegistration(
        step_mat,
        N,
        P,
        M,
        Ny,
        u,
        y_pred,
        y_meas,
        u_hist,
        du,
        move_constraints,
        Matrix{JuMP.ConstraintRef}(undef, P, Ny),
    )

    for l in 1:Ny
        # `dk` is the bias-correction term.
        # It shifts the prediction so the first point matches the current
        # measurement.
        yp0 = sum(
            ctx.s[i, l] * (_hold_u(ctx, ctx.N - i + 1) - _hold_u(ctx, ctx.N - i))
            for i in 1:(ctx.N - 1)
        ) + ctx.s[ctx.N, l] * _hold_u(ctx, 1)

        dk = ctx.y_meas[l] - yp0
        for j in 1:P
            pred_expr = sum(
                ctx.s[i, l] * (_hold_u(ctx, ctx.N - i + j + 1) - _hold_u(ctx, ctx.N - i + j))
                for i in 1:(ctx.N - 1)
            ) + ctx.s[ctx.N, l] * _hold_u(ctx, j + 1) + dk
            ctx.prediction_constraints[j, l] = @constraint(model, ctx.y_pred[j, l] == pred_expr)
        end
    end

    if store_ext
        dmc_store = get!(model.ext, :dmc, Dict{String, Any}())
        dmc_store[String(base_name)] = ctx
    end

    return ctx
end

"""
    update_dmc_state!(ctx; y_meas, u_hist)

Update the fixed measurement and history variables in an existing DMC block.
"""
function update_dmc_state!(ctx::DMCRegistration; y_meas=nothing, u_hist=nothing)
    if y_meas !== nothing
        y_vals = _float_vector(y_meas, ctx.Ny, "y_meas")
        for i in 1:ctx.Ny
            # Re-fix the current measurement.
            JuMP.fix(ctx.y_meas[i], y_vals[i]; force=true)
        end
    end

    if u_hist !== nothing
        hist_len = ctx.N - 1
        hist_vals = _float_vector(u_hist, hist_len, "u_hist")
        for i in 1:hist_len
            # Re-fix the input history.
            JuMP.fix(ctx.u_hist[i], hist_vals[i]; force=true)
        end
    end

    return ctx
end

"""
    summarize_dmc_registration(ctx)

Return a short summary of a registered DMC block.

This is useful when you want to check the horizon lengths and array sizes.
"""
function summarize_dmc_registration(ctx::DMCRegistration)
    return (
        N = ctx.N,
        P = ctx.P,
        M = ctx.M,
        Ny = ctx.Ny,
        u_count = length(ctx.u),
        prediction_rows = size(ctx.y_pred, 1),
        prediction_cols = size(ctx.y_pred, 2),
        du_count = length(ctx.du),
    )
end
