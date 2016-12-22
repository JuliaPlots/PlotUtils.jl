
invisible() = RGBA{Float64}(0.,0.,0.,0.)

# the one-arg cases, meant for single colors
plot_color(s::AbstractString) = parse(RGBA{Float64}, s)

plot_color(s::Symbol; color_library::Symbol = _gradients[1]) = iscgrad_symbol(s, color_library = color_library) ? cgrad(s) : parse(RGBA{Float64}, s)
plot_color(b::Bool) = b ? error("plot_color(true) not allowed.") : invisible()
plot_color(::Void) = invisible()
plot_color(c::Colorant) = convert(RGBA{Float64}, c)
# plot_color(cs::AbstractVector) = RGBA{Float64}[plot_color(c) for c in cs]
# plot_color(cs::AbstractArray) = map(plot_color, cs)
plot_color(grad::ColorGradient) = grad

# no alpha override
plot_color(x, ::Void) = plot_color(x)

# alpha override
plot_color(x, α::Number) = RGBA{Float64}(convert(RGB, plot_color(x)), α)
plot_color(c::Colorant, α::Number) = RGBA{Float64}(red(c), green(c), blue(c), α)

plot_color(s::Symbol, α::Number) = (iscgrad_symbol(s) ? cgrad(s, alpha=α) : RGBA{Float64}(convert(RGB, plot_color(s)), α))
plot_color(grad::ColorGradient, α::Number) = cgrad(grad, alpha=α)

function plot_color(cs::AbstractArray)
    a = Array(RGBA{Float64}, size(cs))
    for i in eachindex(cs)
        a[i] = plot_color(cs[i])
    end
    a
end

# plot_color(cs::AbstractVector, α::Number) = RGBA{Float64}[plot_color(c,α) for c in cs]
function plot_color(cs::AbstractArray, α::Number)
    a = Array(RGBA{Float64}, size(cs))
    for i in eachindex(cs)
        a[i] = plot_color(cs[i], α)
    end
    a
end
    # map!(c -> plot_color(c,α), a, cs))


# convenience conversions from numeric arrays to gradient values
# note: we need the first version because of dispatch
# function plot_color{T<:Number}(zs::AbstractVector{T})
#     grad = cgrad()
#     zmin, zmax = extrema(zs)
#     RGBA{Float64}[grad[(z-zmin)/(zmax-zmin)] for z in zs]
# end
function plot_color{T<:Number}(zs::AbstractArray{T})
    grad = cgrad()
    zmin, zmax = extrema(zs)
    a = Array(RGBA{Float64}, size(zs))
    for i in eachindex(zs)
        a[i] = grad[(zs[i]-zmin)/(zmax-zmin)]
    end
    a
end

# function plot_color{T<:Number}(zs::AbstractVector{T}, α::Number)
#     cs = plot_color(zs)
#     RGBA{Float64}[RGBA{Float64}(convert(RGB, c), α) for c in cs]
# end
function plot_color{T<:Number}(zs::AbstractArray{T}, α::Number)
    cs = plot_color(zs)
    a = Array(RGBA{Float64}, size(zs))
    for i in eachindex(zs)
        a[i] = RGBA{Float64}(convert(RGB, cs[i]), α)
    end
    a
end
