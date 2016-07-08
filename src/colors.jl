
# the one-arg cases, meant for single colors
plot_color(s::AbstractString) = parse(RGBA{Float64}, s)
plot_color(s::Symbol) = haskey(_gradients,s) ? cgrad(s) : parse(RGBA{Float64}, s)
plot_color(b::Bool) = b ? error("plot_color(true) not allowed.") : RGBA{Float64}(0,0,0,0)
plot_color(::Void) = RGBA{Float64}(0,0,0,0)
plot_color(c::Colorant) = convert(RGBA{Float64}, c)
plot_color(cs::AbstractVector) = RGBA{Float64}[plot_color(c) for c in cs]
plot_color(cs::AbstractArray) = map(plot_color, cs)
plot_color(grad::ColorGradient) = grad

# no alpha override
plot_color(x, ::Void) = plot_color(x)

# alpha override
plot_color(x, α::Number) = RGBA{Float64}(convert(RGB, plot_color(x)), α)
plot_color(c::Colorant, α::Number) = RGBA{Float64}(red(c), green(c), blue(c), α)
plot_color(s::Symbol, α::Number) = (haskey(_gradients,s) ? cgrad(s, alpha=α) : RGBA{Float64}(convert(RGB, plot_color(s)), α))
plot_color(cs::AbstractVector, α::Number) = RGBA{Float64}[plot_color(c,α) for c in cs]
plot_color(cs::AbstractArray, α::Number) = (a = Array(RGBA{Float64}, size(cs)); map!(c -> plot_color(c,α), a, cs))
plot_color(grad::ColorGradient, α::Number) = cgrad(grad, alpha=α)
