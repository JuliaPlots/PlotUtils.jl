
"""
colors can be single RGBA colors, or lists of colors (Vector{RGBA}), or gradients
"""

# abstract AbstractPlotColor

# Base.getindex(scheme::AbstractPlotColor, i::Integer) = getColor(scheme, i)

# the one-arg cases, meant for single colors
plot_color(x) = x  # pass through
plot_color(s::AbstractString) = parse(RGBA, s)
plot_color(s::Symbol) = haskey(_gradients,s) ? cgrad(s) : parse(RGBA, s)
plot_color(b::Bool) = b ? RGBA(0,0,0,1) : RGBA(0,0,0,0)
plot_color(::Void) = RGBA(0,0,0,0)
plot_color(c::Colorant) = convert(RGBA, c)
plot_color(cs::AbstractVector) = RGBA[plot_color(c) for c in cs]
plot_color(cs::AbstractArray) = map(plot_color, cs)

# no alpha override
plot_color(x, ::Void) = plot_color(x)

# alpha override
plot_color(x, α::Number) = RGBA(convert(RGB, plot_color(x)), α)
plot_color(c::Colorant, α::Number) = RGBA(red(c), green(c), blue(c), α)
plot_color(s::Symbol, α::Number) = (haskey(_gradients,s) ? cgrad(s, alpha=α) : RGBA(convert(RGB, plot_color(s)), α))
plot_color(cs::AbstractVector, α::Number) = RGBA[plot_color(c,α) for c in cs]
plot_color(cs::AbstractArray, α::Number) = (a = Array(RGBA, size(cs)); map!(c -> plot_color(c,α), a, cs))

# --------------------------------------------------------------

# getColor(scheme::AbstractPlotColor) = getColor(scheme, 1)
# getColorVector(scheme::AbstractPlotColor) = [getColor(scheme)]

# colorscheme(scheme::AbstractPlotColor) = scheme
# colorscheme(s::AbstractString; kw...) = colorscheme(Symbol(s); kw...)
# colorscheme(s::Symbol; kw...) = haskey(_gradients, s) ? ColorGradient(s; kw...) : ColorWrapper(convertColor(s); kw...)
# colorscheme{T<:Real}(s::Symbol, vals::AVec{T}; kw...) = ColorGradient(s, vals; kw...)
# colorscheme(cs::AVec, vs::AVec; kw...) = ColorGradient(cs, vs; kw...)
# colorscheme{T<:Colorant}(cs::AVec{T}; kw...) = ColorGradient(cs; kw...)
# # colorscheme(f::Function; kw...) = ColorFunction(f; kw...)
# colorscheme(v::AVec; kw...) = ColorVector(v; kw...)
# colorscheme(m::AMat; kw...) = size(m,1) == 1 ? map(c->colorscheme(c; kw...), m) : [colorscheme(m[:,i]; kw...) for i in 1:size(m,2)]'
# colorscheme(c::Colorant; kw...) = ColorWrapper(c; kw...)


# --------------------------------------------------------------


# convertColor(c::AbstractString) = parse(Colorant, c)
# convertColor(c::Symbol) = parse(Colorant, string(c))
# convertColor(c::Colorant) = c
# convertColor(cvec::AbstractVector) = map(convertColor, cvec)
# convertColor(c::AbstractPlotColor) = c
# convertColor(v::Void) = RGBA(0,0,0,0)
# convertColor(b::Bool) = b ? RGBA(0,0,0,1) : RGBA(0,0,0,0)

# function convertColor(c, α::Real)
#   c = convertColor(c)
#   RGBA(RGB(getColor(c)), α)
# end
# convertColor(cs::AVec, α::Real) = map(c -> convertColor(c, α), cs)
# convertColor(c, α::Void) = convertColor(c)

# # backup... try to convert
# getColor(c) = convertColor(c)


# --------------------------------------------------------------

# "Continuous gradient between values.  Wraps a list of bounding colors and the values they represent."
# immutable ColorGradient <: AbstractPlotColor
#   colors::Vector
#   values::Vector

#   function ColorGradient{S<:Real}(cs::AVec, vals::AVec{S} = linspace(0, 1, length(cs)); alpha = nothing)
#     if length(cs) == length(vals)
#       return new(convertColor(cs,alpha), collect(vals))
#     end

#     # interpolate the colors for each value
#     vals = merge(linspace(0, 1, length(cs)), vals)
#     grad = ColorGradient(cs)
#     cs = [getColorZ(grad, z) for z in linspace(0, 1, length(vals))]
#     new(convertColor(cs, alpha), vals)
#   end
# end



# Base.getindex(cs::ColorGradient, i::Integer) = getColor(cs, i)
# Base.getindex(cs::ColorGradient, z::Number) = getColorZ(cs, z)


# # create a gradient from a symbol (blues, reds, etc) and vector of boundary values
# function ColorGradient{T<:Real}(s::Symbol, vals::AVec{T} = 0:0; kw...)
#   haskey(_gradients, s) || error("Invalid gradient symbol.  Choose from: ", sort(collect(keys(_gradients))))
#   cs = _gradients[s]
#   if vals == 0:0
#     vals = linspace(0, 1, length(cs))
#   end
#   ColorGradient(cs, vals; kw...)
# end

# function ColorGradient(grad::ColorGradient; alpha = nothing)
#   ColorGradient(convertColor(grad.colors, alpha), grad.values)
# end

# # anything else just gets the default gradient
# function ColorGradient(cw; alpha=nothing)
#     ColorGradient(default_gradient(), alpha=alpha)
# end

# getColor(gradient::ColorGradient, idx::Int) = gradient.colors[mod1(idx, length(gradient.colors))]

# function getColorZ(gradient::ColorGradient, z::Real)
#   cs = gradient.colors
#   vs = gradient.values
#   n = length(cs)
#   @assert n > 0 && n == length(vs)

#   # can we just return the first color?
#   if z <= vs[1] || n == 1
#     return cs[1]
#   end

#   # find the bounding colors and interpolate
#   for i in 2:n
#     if z <= vs[i]
#       return interpolate_rgb(cs[i-1], cs[i], (z - vs[i-1]) / (vs[i] - vs[i-1]))
#     end
#   end

#   # if we get here, return the last color
#   cs[end]
# end

# getColorVector(gradient::ColorGradient) = gradient.colors

# # for 0.3
# Colors.RGBA(c::Colorant) = RGBA(red(c), green(c), blue(c), alpha(c))
# Colors.RGB(c::Colorant) = RGB(red(c), green(c), blue(c))


# --------------------------------------------------------------

# "Wraps a function, taking an index and returning a Colorant"
# immutable ColorFunction <: AbstractPlotColor
#   f::Function
# end

# getColor(scheme::ColorFunction, idx::Int) = scheme.f(idx)

# # --------------------------------------------------------------

# "Wraps a function, taking an z-value and returning a Colorant"
# immutable ColorZFunction <: AbstractPlotColor
#   f::Function
# end

# getColorZ(scheme::ColorZFunction, z::Real) = scheme.f(z)

# --------------------------------------------------------------

# "Wraps a vector of colors... may be vector of Symbol/String/Colorant"
# immutable ColorVector <: AbstractPlotColor
#   v::Vector{Colorant}
#   ColorVector(v::AVec; alpha = nothing) = new(convertColor(v,alpha))
# end

# getColor(scheme::ColorVector, idx::Int) = convertColor(scheme.v[mod1(idx, length(scheme.v))])
# getColorVector(scheme::ColorVector) = scheme.v


# # --------------------------------------------------------------

# "Wraps a single color"
# immutable ColorWrapper <: AbstractPlotColor
#   c::RGBA
#   ColorWrapper(c::Colorant; alpha = nothing) = new(convertColor(c, alpha))
# end

# ColorWrapper(s::Symbol; alpha = nothing) = ColorWrapper(convertColor(parse(Colorant, s), alpha))

# getColor(scheme::ColorWrapper, idx::Int) = scheme.c
# getColorZ(scheme::ColorWrapper, z::Real) = scheme.c
# convertColor(c::ColorWrapper, α::Void) = c.c

# --------------------------------------------------------------




# ----------------------------------------------------------------------------------



# ----------------------------------------------------------------------------------

# TODO: this needs to be added to Plots, maybe in args.jl?

# # converts a symbol or string into a colorant (Colors.RGB), and assigns a color automatically
# function getSeriesRGBColor(c, sp::Subplot, n::Int)

#   if c == :auto
#     c = autopick(sp[:color_palette], n)
#   end

#   # c should now be a subtype of AbstractPlotColor
#   colorscheme(c)
# end
