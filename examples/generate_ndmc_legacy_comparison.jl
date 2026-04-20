using CSV
using DataFrames
using Statistics

length(ARGS) == 4 || error("usage: julia --project=examples generate_ndmc_legacy_comparison.jl <closed_loop_csv> <legacy_csv> <comparison_csv> <summary_csv>")

closed_loop_path, legacy_path, comparison_path, legacy_summary_path = ARGS
legacy_script = joinpath(@__DIR__, "MPC_NDMC.jl")

include(legacy_script)

current_df = CSV.read(closed_loop_path, DataFrame)
current_time = current_df.time

legacy_time_raw = collect(time)
legacy_y_raw = yout
legacy_u = collect(uopt)

legacy_state_idx = [clamp(round(Int, t - first(legacy_time_raw)) + 1, 1, length(legacy_time_raw)) for t in current_time]
legacy_q = [legacy_u[min(fld(round(Int, t), 20) + 1, length(legacy_u))] for t in current_time]

legacy_df = DataFrame(
    time = current_time,
    C1 = legacy_y_raw[legacy_state_idx, 1],
    C2 = legacy_y_raw[legacy_state_idx, 2],
    C3 = legacy_y_raw[legacy_state_idx, 3],
    Cmix = legacy_y_raw[legacy_state_idx, 4],
    cO = legacy_y_raw[legacy_state_idx, 5],
    Q_air = legacy_q,
)
CSV.write(legacy_path, legacy_df)

comp_df = DataFrame(
    time = current_df.time,
    current_C1 = current_df.C1,
    legacy_C1 = legacy_df.C1,
    diff_C1 = current_df.C1 .- legacy_df.C1,
    current_C2 = current_df.C2,
    legacy_C2 = legacy_df.C2,
    diff_C2 = current_df.C2 .- legacy_df.C2,
    current_C3 = current_df.C3,
    legacy_C3 = legacy_df.C3,
    diff_C3 = current_df.C3 .- legacy_df.C3,
    current_Q = current_df.Q_air,
    legacy_Q = legacy_df.Q_air,
    diff_Q = current_df.Q_air .- legacy_df.Q_air,
    current_Cmix = current_df.Cmix,
    legacy_Cmix = legacy_df.Cmix,
    diff_Cmix = current_df.Cmix .- legacy_df.Cmix,
    current_cO = current_df.cO,
    legacy_cO = legacy_df.cO,
    diff_cO = current_df.cO .- legacy_df.cO,
)
CSV.write(comparison_path, comp_df)

rmse(x) = sqrt(mean(x .^ 2))
summary_df = DataFrame(
    metric = [
        "RMSE diff C1", "RMSE diff C2", "RMSE diff C3", "RMSE diff Cmix", "RMSE diff cO", "RMSE diff Q",
        "Max abs diff C1", "Max abs diff C2", "Max abs diff C3", "Max abs diff Cmix", "Max abs diff cO", "Max abs diff Q",
        "Final diff C1", "Final diff C2", "Final diff C3", "Final diff Cmix", "Final diff cO", "Final diff Q",
    ],
    value = [
        rmse(comp_df.diff_C1), rmse(comp_df.diff_C2), rmse(comp_df.diff_C3), rmse(comp_df.diff_Cmix), rmse(comp_df.diff_cO), rmse(comp_df.diff_Q),
        maximum(abs.(comp_df.diff_C1)), maximum(abs.(comp_df.diff_C2)), maximum(abs.(comp_df.diff_C3)), maximum(abs.(comp_df.diff_Cmix)), maximum(abs.(comp_df.diff_cO)), maximum(abs.(comp_df.diff_Q)),
        comp_df.diff_C1[end], comp_df.diff_C2[end], comp_df.diff_C3[end], comp_df.diff_Cmix[end], comp_df.diff_cO[end], comp_df.diff_Q[end],
    ],
)
CSV.write(legacy_summary_path, summary_df)
