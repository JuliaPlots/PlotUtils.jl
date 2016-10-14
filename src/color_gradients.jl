
const _rainbowColors = [colorant"purple", colorant"blue", colorant"green", colorant"orange", colorant"red"]
const _testColors = [colorant"darkblue", colorant"blueviolet",  colorant"darkcyan",colorant"green",
                     darken(colorant"yellow",0.3), colorant"orange", darken(colorant"red",0.2)]

const _gradients = Dict{Symbol,Vector{RGBA{Float64}}}(
    :blues        => [colorant"lightblue", colorant"darkblue"],
    :reds         => [colorant"lightpink", colorant"darkred"],
    :greens       => [colorant"lightgreen", colorant"darkgreen"],
    :redsblues    => [colorant"darkred", RGB(0.8,0.85,0.8), colorant"darkblue"],
    :bluesreds    => [colorant"darkblue", RGB(0.8,0.85,0.8), colorant"darkred"],
    :heat         => [colorant"lightyellow", colorant"orange", colorant"darkred"],
    :grays        => [RGB(.95,.95,.95),RGB(.05,.05,.05)],
    :rainbow      => _rainbowColors,
    :lightrainbow => map(lighten, _rainbowColors),
    :darkrainbow  => map(darken, _rainbowColors),
    :darktest     => _testColors,
    :lighttest    => map(c -> lighten(c, 0.3), _testColors),
  )

function register_gradient_colors{C<:Colorant}(name::Symbol, colors::AbstractVector{C})
    _gradients[name] = colors
end

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

cgrad_colors(s::Symbol) = copy(_gradients[s])
cgrad_colors(grad::ColorGradient) = copy(grad.colors)
cgrad_colors(cs::Vector{RGBA{Float64}}) = cs
cgrad_colors(cs::AbstractVector) = RGBA{Float64}[plot_color(c) for c in cs]

function _color_list(arg, ::Void)
    cgrad_colors(arg)
end
function _color_list(arg, alpha)
    colors = cgrad_colors(arg)
    for i in eachindex(colors)
        colors[i] = RGBA{Float64}(convert(RGB{Float64}, colors[i]), alpha)
    end
    colors
end

# construct a ColorGradient when given explicit values
function cgrad(arg, values; alpha = nothing)
    colors = _color_list(arg, alpha)
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
function cgrad(arg; alpha = nothing, scale = :identity)
    colors = _color_list(arg, alpha)
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

const _default_gradient = Ref(:inferno)

# the default gradient
cgrad(; kw...) = cgrad(_default_gradient[]; kw...)


cvec(s::Symbol, n::Integer = 10; kw...) = cvec(cgrad(s; kw...), n)
cvec(grad::ColorGradient, n::Integer = 10; kw...) = RGBA{Float64}[grad[z] for z in linspace(0,1,n)]

function sample_evenly(v::AbstractVector, n::Integer = length(v))
    idx = Int[round(Int, x) for x in linspace(1, length(v), n)]
    v[idx]
end

include("gradients/matplotlib.jl")
include("gradients/cmocean.jl")
