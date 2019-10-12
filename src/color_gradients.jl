
mutable struct ColorLibrary
    defaults::Dict{Symbol, Symbol}
    lib::Dict{Symbol, Vector{RGBA{Float64}}}
    ColorLibrary(defaults = Dict(:default => :sequential), lib = Dict{Symbol, Vector{RGBA{Float64}}}()) = new(defaults, lib)
end

ColorLibrary(lib::Dict{Symbol, Vector{RGBA{Float64}}}) =
    ColorLibrary(Dict(:default => first(keys(lib))), lib)

function ColorLibrary(lib::Dict{Symbol, Vector{RGBA{Float64}}}, default::Symbol)
    in(default, keys(lib)) || error("There is no gradient named $default in lib")
    ColorLibrary(Dict(:default => default), lib)
end

default_cgrad(cl::Symbol = _gradients[1]; kwargs...) = default_cgrad(color_libraries[cl]; kwargs...)

function default_cgrad(cl::ColorLibrary; default = nothing, sequential = nothing, diverging = nothing)
    default == nothing || (cl.defaults[:default] = default)
    sequential == nothing || (cl.defaults[:sequential] = sequential)
    diverging == nothing || (cl.defaults[:diverging] = diverging)
end

const color_libraries = Dict{Symbol, ColorLibrary}()

function getgradient(gradient::Symbol = :default, clib::Symbol = _gradients[1])
    haskey(color_libraries, clib) || error("There is no color library named $clib . Use clibraries() to get a list of available color libraries")
    cl = color_libraries[clib]
    getgradient(gradient, cl)
end

function getgradient(gradient::Symbol, cl::ColorLibrary)
    while haskey(cl.defaults, gradient)
        gradient = cl.defaults[gradient]
    end
    haskey(cl.lib, gradient) && return cl.lib[gradient]

    potentials = [name for (name, library) in color_libraries if haskey(library.lib, gradient)]
    length(potentials) == 0 && error("There is no gradient named $gradient . Use cgradients() to get a list of gradients in the current color library, clibraries() to get a list of available color libraries")
    length(potentials) > 1 && @warn("$gradient is found in more than one library: $(join(potentials, ", ")). Choosing $(potentials[1])")
    color_libraries[potentials[1]][gradient]
end

getindex(cl::ColorLibrary, key::Symbol) = getgradient(key, cl)

const _categorical_gradients = Symbol[]
is_categorical(grad) = grad in _categorical_gradients

function register_gradient_colors(name::Symbol, colors::AbstractVector{C}, color_library::Symbol = :default, is_categorical = false) where {C<:Colorant}
    is_categorical && push!(_categorical_gradients, name)
    haskey(color_libraries, color_library) || register_color_library(color_library, ColorLibrary())
    color_libraries[color_library].lib[name] = colors
end

function register_color_library(name::Symbol, color_library::ColorLibrary = ColorLibrary())
    color_libraries[name] = color_library
end

"""
    clibrary(grad::Symbol)

Set the active color library. A list of possible libraries can be printed with `clibraries()`
"""
function clibrary(grad::Symbol)
    haskey(color_libraries, grad) || error("$grad is not a defined color library, valid choices are: "*join([":$(library)"  for library in keys(color_libraries)], ", "))
    _gradients[1] = grad
end

const _gradients = [:Plots]

"""
    clibraries()

List the available color libraries on the system
"""
clibraries() = collect(keys(color_libraries))

"""
    cgradients([color_library::Symbol])

List available color gradients in color_library (defaults to the currently loaded library)
"""
cgradients(color_library::Symbol = _gradients[1]) = collect(keys(color_libraries[color_library].lib))

# --------------------------------------------------------------------------



"Continuous gradient between values.  Wraps a list of bounding colors and the values they represent."
struct ColorGradient
  colors::Vector{RGBA{Float64}}
  values::Vector{Float64}
end

# if the values aren't passed in, pass to the cgrad method for processing
ColorGradient(colors; kw...) = cgrad(colors; kw...)
Base.length(::ColorGradient) = @error "length(::ColorGradient) was called, possibly due to a spuriously broadcast call to a function accepting ColorGradient (e.g. PlotUtils.plot_color). Please open an issue on the library you're using (e.g. Plots)"
Base.:(==)(g1::ColorGradient, g2::ColorGradient) = (g1.colors == g2.colors) && (g1.values == g2.values)
Base.hash(g::ColorGradient) = hash(g.colors) | hash(g.values)

Base.getindex(gradient::ColorGradient, idx::Integer) = gradient.colors[mod1(idx, length(gradient.colors))]

function Base.getindex(gradient::ColorGradient, z::Number)
    isnan(z) && return invisible()
    cs = gradient.colors
    vs = gradient.values
    n = length(cs)
    @assert n > 0 && n == length(vs)

    # can we just return the first color?
    if z <= vs[1] || n == 1
        return cs[1]
    end

    # find the bounding colors and interpolate
    for i in 2:n
        if z <= vs[i]
            return interpolate_rgb(cs[i-1], cs[i], (z - vs[i-1]) / (vs[i] - vs[i-1]))
        end
    end

    # if we get here, return the last color
    cs[end]
end


# --------------------------------------------------------------------------

function cgrad_reverse(s::Symbol)
    str = string(s)
    rev, s = if length(str) > 2 && str[end-1:end] == "_r"
        (true, Symbol(str[1:end-2]))
    else
        (false, s)
    end
end

function iscgrad_symbol(s::Symbol)
    rev, s = cgrad_reverse(s)
    lib = color_libraries[_gradients[1]]
    haskey(lib.lib,s) && return true
    haskey(lib.defaults,s) && return true
    for library in values(color_libraries)
        haskey(library.lib, s) && return true
    end
    return false
end

function cgrad_colors(s::Symbol; color_library = _gradients[1])
    rev, s = cgrad_reverse(s)
    if rev
        reverse(getgradient(s, color_library))
    else
        getgradient(s, color_library)
    end
end

cgrad_colors(grad::ColorGradient) = copy(grad.colors)
cgrad_colors(cs::Vector{RGBA{Float64}}) = cs
cgrad_colors(cs::AbstractVector) = RGBA{Float64}[plot_color(c) for c in cs]

function _color_list(arg, ::Nothing)
    cgrad_colors(arg)
end

_color_list(arg, alpha) = RGBA{Float64}.(convert.(RGB{Float64}, cgrad_colors(arg)), alpha)

cgrad(arg::Symbol, cl::Symbol, values; kw...) = cgrad(cgrad_colors(arg, color_library = cl), values; kw...)

# construct a ColorGradient when given explicit values
function cgrad(arg, values; alpha = nothing)
    colors = _color_list(arg, alpha)
    values = if length(colors) == length(values) && values[1] == 0 && values[end] == 1
        values
    else
        # merge values into the default range, then recompute color list and return vals
        # vals = merge(collect(range(0; stop=1, length=length(colors))), collect(values))
        vals = sort(unique(vcat(range(0; stop=1, length=length(colors)), values)))
        grad = ColorGradient(colors)
        colors = RGBA{Float64}[grad[z] for z in range(0; stop=1, length=length(vals))]
        vals
    end

    # construct and return the gradient
    ColorGradient(colors, values)
end

cgrad(arg::Symbol, cl::Symbol; kw...) = cgrad(cgrad_colors(arg, color_library = cl); kw...)

# construct a ColorGradient automatically
function cgrad(arg; alpha = nothing, scale = :identity)
    colors = _color_list(arg, alpha)
    values = if scale in (:log, :log10)
        log10.(range(1; stop=10, length=length(colors)))
    elseif scale == :log2
        log2.(range(1; stop=2, length=length(colors)))
    elseif scale == :ln
        log.(range(1; stop=pi, length=length(colors)))
    elseif scale in (:exp, :exp10)
        (exp10.(range(0; stop=1, length=length(colors))) .- 1) ./ 9
    elseif isa(arg, ColorGradient)
        arg.values
    else
        range(0; stop=1, length=length(colors))
    end

    # construct and return the gradient
    ColorGradient(colors, values)
end



# the default gradient
cgrad(; kw...) = cgrad(:default; kw...)


cvec(s::Symbol, n::Integer = 10; kw...) = cvec(cgrad(s; kw...), n)
cvec(grad::ColorGradient, n::Integer = 10; kw...) = RGBA{Float64}[grad[z] for z in range(0; stop=1, length=n)]

function sample_evenly(v::AbstractVector, n::Integer = length(v))
    idx = Int[round(Int, x) for x in range(1; stop=length(v), length=n)]
    v[idx]
end


# allows passing a ColorGradient to rgba_string and get a useful response by picking the first color - introduced because the plotly backend to Plots uses this functionality
rgba_string(cg::ColorGradient) = rgba_string(cg[1])

include("gradients/matplotlib.jl")
include("gradients/cmocean.jl")
include("gradients/colorbrewer.jl")
include("gradients/colorcet.jl")
include("gradients/misc.jl")
