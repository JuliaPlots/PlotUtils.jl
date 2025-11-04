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

const DEBUG = true

@testset "downstream Plots" begin
    script = tempname()
    write(
        script,
        """
        include(joinpath("$(@__DIR__)", "downstream.jl"))
        develop_stable_Plots()
        using Plots

        # test basic plots creation & display (Plots tests are too long to run)
        withenv("GKSwstype" => "nul") do
            @time for i in 1:length(Plots._examples)
                i ∈ Plots._backend_skips[:gr] && continue  # skip unsupported examples
                Plots._examples[i].imports ≡ nothing || continue  # skip examples requiring optional test deps
                show(devnull, Plots.test_examples(:gr, i; disp = false))  # trigger display logic
            end
        end
        exit()
        """,
    )
    DEBUG && print(read(script, String))
    @test run(```$(Base.julia_cmd()) $script```) |> success
    rm(script)
end

const EXTENDED = tryparse(Bool, get(ENV, "CI", "false")) === true  # extended test in CI

@testset "downstream Makie" begin
    # EXTENDED && Pkg.test("CairoMakie")

    script = tempname()
    write(
        script,
        """
        include(joinpath("$(@__DIR__)", "downstream.jl"))
        develop_stable_Makie($EXTENDED)
        using CairoMakie

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

            fn = tempname() * ".png"
            save(fn, f)
            @assert isfile(fn)
            rm(fn)
        end

        exit()
        """,
    )
    DEBUG && print(read(script, String))
    @test run(```$(Base.julia_cmd()) $script```) |> success
    rm(script)
end

@testset "adaptative test Plots" begin
    script = tempname()
    write(
        script,
        """
        include(joinpath("$(@__DIR__)", "downstream.jl"))
        develop_stable_Plots()
        using Plots

        let test_funcs = [
            (x -> 2, 0, 1),
            (x -> x, -1, 1),
            (x -> 5x, -1, 1),
            (x -> 1 / x, 0, 1),
            (x -> 1 / x, -1, 1),
            (x -> tan(x), 0, 4),
            (x -> -1 / x, 0, 1),
            (x -> 1 / abs(x), -1, 1),
            (x -> 1.0e6x, -1, 1),
            (x -> 1.0e50x, -1, 1),
            (x -> log(1 + sin(cos(x))), -6, 6),
            (x -> sin(x^3) + cos(x^3), 0, 6.28),
            (x -> sin(x), -5, 200),
            (x -> sin(1 / x), -2, 2),
            (x -> sin(1 / x), 0, 2),
            (x -> sin(x^4), -4, 4),
            (x -> sin(300x), -4, 4),
            (x -> 1 + x^2 + 0.0125log(abs(1 - 3(x - 1))), -2, 2),
            (x -> sin(exp(x)), -6, 6),
            (x -> 1 / sin(x), -10, 10),
            (x -> sin(x) / x, -6, 6),
            (x -> tan(x^3 - x + 1) + 1 / (x + 3exp(x)), -2, 2),
        ]

        plots = [plot(func...) for func in test_funcs]

        pitchfork = plot(x -> sqrt(x), 0, 10)
        plot!(pitchfork, x -> -sqrt(x), 0, 10)
        plot!(pitchfork, x -> 0, -10, 0)
        plot!(pitchfork, x -> 0, 0, 10, linestyle = :dash)
        push!(plots, pitchfork)

        np = length(plots)
        m = ceil(Int, √(np)) + 1
        n, r = divrem(np, m)
        r == 0 || (n += 1)
        append!(plots, [plot() for _ in 1:(m * n - np)])

        @assert length(plots) == m * n

        png(plot(plots...; layout = (m, n), size = (m * 600, n * 400)), "grid")
        """,
    )
    DEBUG && print(read(script, String))
    @test run(```$(Base.julia_cmd()) $script```) |> success
    rm(script)
end
