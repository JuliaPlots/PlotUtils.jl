using Pkg, PlotUtils

LibGit2 = Pkg.GitTools.LibGit2
TOML = Pkg.TOML

Plots_jl = joinpath(mkpath(tempname()), "Plots.jl")
Plots_toml = joinpath(Plots_jl, "Project.toml")

# clone and checkout the latest stable version of Plots
depot = joinpath(first(DEPOT_PATH), "registries", "General", "P", "Plots", "Versions.toml")
stable = maximum(VersionNumber.(keys(TOML.parse(read(depot, String)))))
for i ∈ 1:6
    try
        global repo = Pkg.GitTools.ensure_clone(stdout, Plots_jl, "https://github.com/JuliaPlots/Plots.jl")
        break
    catch err
        @warn err
        sleep(20i)
    end
end
@assert isfile(Plots_toml) "spurious network error: clone failed, bailing out"
tag = LibGit2.GitObject(repo, "v$stable")
hash = string(LibGit2.target(tag))
LibGit2.checkout!(repo, hash)

# fake the supported PlotUtils version for testing (for `Pkg.develop`)
plotutils_version = Pkg.Types.read_package(normpath(@__DIR__, "..", "Project.toml")).version
toml = TOML.parse(read(Plots_toml, String))
toml["compat"]["PlotUtils"] = string(plotutils_version)
open(Plots_toml, "w") do io
  TOML.print(io, toml)
end
Pkg.develop(path=Plots_jl)
Pkg.status(["PlotUtils", "Plots"])

# test basic plots creation & display
using Plots, Test

@time for i ∈ 1:length(Plots._examples)
  i ∈ Plots._backend_skips[:gr] && continue  # skip unsupported examples
  Plots._examples[i].imports ≡ nothing || continue  # skip examples requiring optional test deps
  show(devnull, Plots.test_examples(:gr, i; disp = false))  # trigger display logic
end
