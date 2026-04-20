using Pkg

ENV["JULIA_PKG_PRECOMPILE_AUTO"] = "0"
ENV["JULIA_NUM_PRECOMPILE_TASKS"] = "1"

function pid_is_alive(pid::Integer)
    pid <= 0 && return false
    rc = ccall(:kill, Cint, (Cint, Cint), pid, 0)
    return rc == 0 || Base.Libc.errno() == Base.Libc.EPERM
end

function clear_stale_precompile_pidfiles!()
    compiled_root = joinpath(homedir(), ".julia", "compiled", "v$(VERSION.major).$(VERSION.minor)")
    isdir(compiled_root) || return nothing
    for (root, _, files) in walkdir(compiled_root)
        for file in files
            endswith(file, ".pidfile") || continue
            pidfile = joinpath(root, file)
            text = try
                strip(read(pidfile, String))
            catch
                ""
            end
            pid = try
                parse(Int, split(text)[1])
            catch
                0
            end
            pid_is_alive(pid) && continue
            rm(pidfile; force = true)
        end
    end
    return nothing
end

repo_root = normpath(joinpath(@__DIR__, ".."))
examples_env = @__DIR__
Pkg.activate(examples_env; io = devnull)
isfile(joinpath(examples_env, "Manifest.toml")) || Pkg.instantiate(; io = devnull)
clear_stale_precompile_pidfiles!()

using EOptInterface
using ModelingToolkit, Unitful, IfElse, NativeNaNMath
using ModelingToolkit: t_nounits as t, D_nounits as D
using OrdinaryDiffEq, DiffEqCallbacks
using Logging

with_logger(Logging.NullLogger()) do
    include(joinpath(repo_root, "notebooks", "uconn_wrrf_asm3", "ASM3_UConn_BNR_Unit_CSTR_Emulation.jl"))
end

include(joinpath(repo_root, "notebooks", "uconn_wrrf_asm3", "ASM3_MPC_Utilities.jl"))

using Printf
using Plots
using Measures
using DataFrames
using JuMP, Ipopt
using MathOptInterface
const MOI = MathOptInterface

include(joinpath(repo_root, "notebooks", "eoi_publication_plots.jl"))
apply_eoi_publication_style!()

println("Starting closed-loop UConn ASM3 MPC...")
flush(stderr)

cfg_hi = MPCConfig(show_status = true, status_stream = :stderr, verbose = false)
cfg_lo = MPCConfig(sp = 2.0, show_status = true, status_stream = :stderr, verbose = false)

uconn_closed_loop = run_uconn_closed_loop_mpc(sys, sol; cfg_hi = cfg_hi, cfg_lo = cfg_lo)

println("final solve status = ", isempty(uconn_closed_loop.log.statuses) ? "n/a" : last(uconn_closed_loop.log.statuses))
println("final objective = ", isempty(uconn_closed_loop.log.objectives) ? NaN : round(last(uconn_closed_loop.log.objectives); digits = 3))
display(last(uconn_closed_loop.log_df, min(10, nrow(uconn_closed_loop.log_df))))

p = asm3_plot_closed_loop_pair(uconn_closed_loop.log_df)
outdir = joinpath(repo_root, "examples", "generated")
mkpath(outdir)
outfile = joinpath(outdir, "uconn_asm3_closed_loop_pair.png")
savefig(p, outfile)
println("saved figure: ", outfile)
