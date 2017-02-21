type ColorLibrary
    defaults::Dict{Symbol, Symbol}
    lib::Dict{Symbol, Vector{RGBA{Float64}}}
    ColorLibrary(defaults = Dict(:default => :sequential), lib = Dict{Symbol, Vector{RGBA{Float64}}}()) = new(defaults, lib)
end

ColorLibrary(lib::Dict{Symbol, Vector{RGBA{Float64}}}) =
    ColorLibrary(Dict(:default => keys(lib)[1]), lib)

function ColorLibrary(lib::Dict{Symbol, Vector{RGBA{Float64}}}, default::Symbol)
    in(default, keys(lib)) || error("There is no gradient named $default in lib")
    ColorLibrary(Dict(:default => default), lib)
end

function cgraddefaults(cl::ColorLibrary; default = nothing, sequential = nothing, diverging = nothing)
    default == nothing || (cl.defaults[:default] = default)
    sequential == nothing || (cl.defaults[:sequential] = sequential)
    diverging == nothing || (cl.defaults[:diverging] = diverging)
end

const color_libraries = Dict{Symbol, ColorLibrary}()

function getgradient(gradient::Symbol = :default, clibrary::Symbol = _gradients[1])
    haskey(color_libraries, clibrary) || error("There is no color library named $clibrary . Use clibraries() to get a list of available color libraries")
    cl = color_libraries[clibrary]
    getgradient(gradient, cl)
end

function getgradient(gradient::Symbol, cl::ColorLibrary)
    while haskey(cl.defaults, gradient)
        gradient = cl.defaults[gradient]
    end
    haskey(cl.lib, gradient) && return cl.lib[gradient]

    potentials = [name for (name, library) in color_libraries if haskey(library.lib, gradient)]
    length(potentials) == 0 && error("There is no gradient named $gradient . Use cgradients() to get a list of gradients in the current color library, clibraries() to get a list of available color libraries")
    length(potentials) > 1 && warn("$gradient is found in more than one library: $(join(potentials, ", ")). Choosing $(potentials[1])")
    color_libraries[potentials[1]][gradient]
end

getindex(cl::ColorLibrary, key::Symbol) = getgradient(key, cl)

function register_gradient_colors{C<:Colorant}(name::Symbol, colors::AbstractVector{C}, color_library::Symbol = :default)
    haskey(color_libraries, color_library) || register_color_library(color_library, ColorLibrary())
    color_libraries[color_library].lib[name] = colors
end

function register_color_library(name::Symbol, color_library::ColorLibrary)
    color_libraries[name] = color_library
end

"""
    clibrary(grad::Symbol)

Set the active color library. A list of possible libraries can be printed with `clibraries()`
"""
function set_clibrary(grad::Symbol)
    haskey(color_libraries, grad) || error("$grad is not a defined color library, valid choices are: "*join([":$(library)"  for library in keys(color_libraries)], ", "))
    _gradients[1] = grad
end

function clibrary(grad::Symbol)
    haskey(color_libraries, grad) || error("$grad is not a defined color library, valid choices are: "*join([":$(library)"  for library in keys(color_libraries)], ", "))
    color_libraries[grad]
end

const _rainbowColors = [colorant"purple", colorant"blue", colorant"green", colorant"orange", colorant"red"]
const _testColors = [colorant"darkblue", colorant"blueviolet",  colorant"darkcyan",colorant"green",
                     darken(colorant"yellow",0.3), colorant"orange", darken(colorant"red",0.2)]

const Plots_internal = ColorLibrary(Dict(:default => :heat), Dict(
    :blues        => [colorant"lightblue", colorant"darkblue"],
    :reds         => [colorant"lightpink", colorant"darkred"],
    :greens       => [colorant"lightgreen", colorant"darkgreen"],
    :redsblues    => [colorant"darkred", RGB(0.8,0.85,0.8), colorant"darkblue"],
    :bluesreds    => [colorant"darkblue", RGB(0.8,0.85,0.8), colorant"darkred"],
    :heat         => [colorant"lightyellow", colorant"orange", colorant"darkred"],
    :grays        => [RGB(.05,.05,.05),RGB(.95,.95,.95)],
    :rainbow      => _rainbowColors,
    :lightrainbow => map(lighten, _rainbowColors),
    :darkrainbow  => map(darken, _rainbowColors),
    :darktest     => _testColors,
    :lighttest    => map(c -> lighten(c, 0.3), _testColors),
  ))


register_color_library(:Plots_internal, Plots_internal)
const _gradients = [:matplotlib]

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
immutable ColorGradient
  colors::Vector{RGBA{Float64}}
  values::Vector{Float64}
end

# if the values aren't passed in, pass to the cgrad method for processing
ColorGradient(colors; kw...) = cgrad(colors; kw...)


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

function iscgrad_symbol(s::Symbol; color_library = _gradients[1])
    rev, s = cgrad_reverse(s)
    lib = isa(color_library, Symbol) ? color_libraries[color_library] : color_library
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

function _color_list(arg, ::Void; color_library = _gradients[1])
    cgrad_colors(arg; color_library = color_library)
end

function _color_list(arg, alpha; color_library = _gradients[1])
    colors = cgrad_colors(arg; color_library = color_library)
    for i in eachindex(colors)
        colors[i] = RGBA{Float64}(convert(RGB{Float64}, colors[i]), alpha)
    end
    colors
end

# construct a ColorGradient when given explicit values
function cgrad(arg, values; alpha = nothing, color_library = _gradients[1])
    colors = _color_list(arg, alpha; color_library = color_library)
    values = if length(colors) == length(values) && values[1] == 0 && values[end] == 1
        values
    else
        # merge values into the default range, then recompute color list and return vals
        # vals = merge(collect(linspace(0, 1, length(colors))), collect(values))
        vals = sort(unique(vcat(linspace(0, 1, length(colors)), values)))
        grad = ColorGradient(colors)
        colors = RGBA{Float64}[grad[z] for z in linspace(0, 1, length(vals))]
        vals
    end

    # construct and return the gradient
    ColorGradient(colors, values)
end

# construct a ColorGradient automatically
function cgrad(arg; alpha = nothing, scale = :identity, color_library = _gradients[1])
    colors = _color_list(arg, alpha, color_library = color_library)
    values = if scale in (:log, :log10)
        log10(linspace(1,10,30))
    elseif scale == :log2
        log2(linspace(1,2,30))
    elseif scale == :ln
        log(linspace(1,pi,30))
    elseif scale in (:exp, :exp10)
        (exp10(linspace(0,1,30)) - 1) / 9
    elseif isa(arg, ColorGradient)
        arg.values
    else
        linspace(0, 1, length(colors))
    end

    # construct and return the gradient
    ColorGradient(colors, values)
end

# the default gradient
cgrad(; kw...) = cgrad(:default; kw...)


cvec(s::Symbol, n::Integer = 10; kw...) = cvec(cgrad(s; kw...), n)
cvec(grad::ColorGradient, n::Integer = 10; kw...) = RGBA{Float64}[grad[z] for z in linspace(0,1,n)]

function sample_evenly(v::AbstractVector, n::Integer = length(v))
    idx = Int[round(Int, x) for x in linspace(1, length(v), n)]
    v[idx]
end

include("gradients/matplotlib.jl")
include("gradients/cmocean.jl")
include("gradients/colorbrewer.jl")
