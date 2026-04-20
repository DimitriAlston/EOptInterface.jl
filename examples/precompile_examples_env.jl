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
Pkg.activate(@__DIR__; io = devnull)
isfile(joinpath(@__DIR__, "Manifest.toml")) || Pkg.instantiate(; io = devnull)
clear_stale_precompile_pidfiles!()

println("Precompiling examples environment from: ", Base.active_project())
Pkg.precompile()
println("examples-env-precompile-ok")
