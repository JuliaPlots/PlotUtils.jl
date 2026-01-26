using Pkg

LibGit2 = Pkg.GitTools.LibGit2
TOML = Pkg.TOML

failsafe_clone_checkout(path, url, pkg = nothing) = begin
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
    obj = LibGit2.GitObject(repo, "v$stable")
    hash = if isa(obj, LibGit2.GitTag)
        LibGit2.target(obj)
    else
        LibGit2.GitHash(obj)
    end |> string
    LibGit2.checkout!(repo, hash)

    toml = if pkg ≢ nothing && (fn = joinpath(path, pkg, "Project.toml")) |> isfile  # monorepo layout
        fn
    elseif (fn = joinpath(path, "Project.toml")) |> isfile  # single package toplevel
        fn
    end
    @assert isfile(toml) "$toml does not exist, bailing out !"
    toml
end

fake_supported_version!(path, toml) = begin
    # fake the supported PlotUtils version for testing (for `Pkg.develop`)
    PlotUtils_version = Pkg.Types.read_package(normpath(@__DIR__, "..", "Project.toml")).version
    parsed_toml = TOML.parse(read(toml, String))
    parsed_toml["compat"]["PlotUtils"] = string(PlotUtils_version)
    open(toml, "w") do io
        TOML.print(io, parsed_toml)
    end
    nothing
end

scratch_env_with_PlotUtils() = begin
    Pkg.activate(; temp = true)
    Pkg.develop(path = normpath(@__DIR__, ".."))  # PlotUtils
    nothing
end

develop_stable_Plots() = begin
    scratch_env_with_PlotUtils()
    tmpd = mktempdir()
    Plots_jl = joinpath(tmpd, "Plots.jl")

    toml = failsafe_clone_checkout(Plots_jl, "https://github.com/JuliaPlots/Plots.jl", "Plots")
    fake_supported_version!(Plots_jl, toml)

    Pkg.develop(path = dirname(toml))
    Pkg.status(["PlotUtils", "Plots"])
    nothing
end

develop_stable_Makie(extended = false) = begin
    scratch_env_with_PlotUtils()
    tmpd = mktempdir()
    Makie_jl = joinpath(tmpd, "Makie.jl")

    toml = failsafe_clone_checkout(Makie_jl, "https://github.com/MakieOrg/Makie.jl", "Makie")
    fake_supported_version!(Makie_jl, toml)

    Pkg.develop(path = joinpath(Makie_jl, "ComputePipeline"))
    Pkg.develop(path = joinpath(Makie_jl, "Makie"))
    extended && Pkg.develop(path = joinpath(Makie_jl, "ReferenceTests"))
    Pkg.develop(path = joinpath(Makie_jl, "CairoMakie"))
    # Pkg.develop(path = joinpath(Makie_jl, "GLMakie"))
    Pkg.status(["PlotUtils", "Makie"])
    nothing
end
