function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    @assert precompile(
        Tuple{
            Core.kwftype(typeof(optimize_ticks)),
            NamedTuple{(:k_min, :k_max),Tuple{Int,Int}},
            typeof(optimize_ticks),
            Float64,
            Float64,
        },
    )
    for C in (RGB, RGBA), T in (Colors.FixedPointNumbers.N0f8, Float32, Float64)
        @assert precompile(cgrad, (ColorScheme{Vector{C{T}},String,String},))
        for V in (Vector{Float64}, typeof(get_range(3)))
            @assert precompile(cgrad, (ColorScheme{Vector{C{T}},String,String}, V))
            @assert precompile(
                CategoricalColorGradient,
                (ColorScheme{Vector{C{T}},String,String}, V),
            )
            @assert precompile(
                ContinuousColorGradient,
                (ColorScheme{Vector{C{T}},String,String}, V),
            )
        end
    end
    for T in (Float64, Int)
        @assert precompile(zscale, (Vector{T},))
    end
    @assert precompile(optimize_ticks, (Int, Int))
    @assert precompile(optimize_datetime_ticks, (Int, Int))
    @assert precompile(adapted_grid, (Function, Tuple{Int,Int}))
end
