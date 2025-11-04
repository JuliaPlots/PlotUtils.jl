using Pkg, PlotUtils, Test

LibGit2 = Pkg.GitTools.LibGit2
TOML = Pkg.TOML

failsafe_clone_checkout(path, toml, url) = begin
    local repo
    for i in 1:6
        try
            repo = Pkg.GitTools.ensure_clone(stdout, path, url)
            break
        catch err
            @warn err
            sleep(20i)
        end
    end

    @assert isfile(toml) "spurious network error: clone failed, bailing out"

    name, _ = splitext(basename(url))
    registries = joinpath(first(DEPOT_PATH), "registries")
    general = joinpath(registries, "General")
    versions = joinpath(general, name[1:1], name, "Versions.toml")
    if !isfile(versions)
        mkpath(general)
        run(setenv(`tar xf $general.tar.gz`; dir = general))
    end
    @assert isfile(versions)

    stable = maximum(VersionNumber.(keys(TOML.parse(read(versions, String)))))
    tag = LibGit2.GitObject(repo, "v$stable")
    hash = string(LibGit2.target(tag))
    LibGit2.checkout!(repo, hash)
    nothing
end

fake_supported_version!(path, toml) = begin
    # fake the supported PlotUtils version for testing (for `Pkg.develop`)
    PlotUtils_version =
        Pkg.Types.read_package(normpath(@__DIR__, "..", "Project.toml")).version
    parsed_toml = TOML.parse(read(toml, String))
    parsed_toml["compat"]["PlotUtils"] = string(PlotUtils_version)
    open(toml, "w") do io
        TOML.print(io, parsed_toml)
    end
    nothing
end

develop_stable_Plots() = begin
    tmpd = mktempdir()
    Plots_jl = joinpath(tmpd, "Plots.jl")
    toml = joinpath(Plots_jl, "Project.toml")

    failsafe_clone_checkout(Plots_jl, toml, "https://github.com/JuliaPlots/Plots.jl")
    fake_supported_version!(Plots_jl, toml)

    Pkg.develop(path = Plots_jl)
    Pkg.status(["PlotUtils", "Plots"])
    nothing
end

develop_stable_Makie(extended = false) = begin
    tmpd = mktempdir()
    Makie_jl = joinpath(tmpd, "Makie.jl")
    toml = joinpath(Makie_jl, "Makie", "Project.toml")

    failsafe_clone_checkout(Makie_jl, toml, "https://github.com/MakieOrg/Makie.jl")
    fake_supported_version!(Makie_jl, toml)

    Pkg.develop(path = joinpath(tmpd, "Makie.jl", "ComputePipeline"))
    Pkg.develop(path = joinpath(tmpd, "Makie.jl", "Makie"))
    extended && Pkg.develop(path = joinpath(tmpd, "Makie.jl", "ReferenceTests"))
    Pkg.develop(path = joinpath(tmpd, "Makie.jl", "CairoMakie"))
    # Pkg.develop(path = joinpath(tmpd, "Makie.jl", "GLMakie"))
    Pkg.status(["PlotUtils", "Makie"])
    nothing
end

develop_stable_Plots()
using Plots

@testset "downstream Plots" begin
    # test basic plots creation & display (Plots tests are too long to run)
    withenv("GKSwstype" => "nul") do
        @time for i in 1:length(Plots._examples)
            i ∈ Plots._backend_skips[:gr] && continue  # skip unsupported examples
            Plots._examples[i].imports ≡ nothing || continue  # skip examples requiring optional test deps
            show(devnull, Plots.test_examples(:gr, i; disp = false))  # trigger display logic
        end
    end
    @test true
end

extended = tryparse(Bool, get(ENV, "CI", "false")) === true  # extended test in CI

develop_stable_Makie(extended)
using CairoMakie

@testset "downstream Makie" begin
    Pkg.test("Makie")

    let f = Figure()  # taken from https://docs.makie.org/dev/reference/blocks/axis#yscale
        for (i, scale) in enumerate([identity, log10, log2, log, sqrt, Makie.logit])
            row, col = fldmod1(i, 3)
            Axis(
                f[row, col], yscale = scale, title = string(scale),
                yminorticksvisible = true, yminorgridvisible = true,
                yminorticks = IntervalsBetween(5)
            )
            lines!(range(0.01, 0.99, length = 200))
        end

        fn = "$(tempname()).png"
        save(fn, f)
        @test isfile(fn)
        rm(fn)
    end

    # extended && Pkg.test("CairoMakie")
end
