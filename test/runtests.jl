using PlotUtils
using Base.Test

# TODO: real tests

# ----------------------
# colors

const C = RGBA{Float64}

@test plot_color(nothing) == C(0,0,0,0)
@test plot_color(false) == C(0,0,0,0)
@test_throws ErrorException plot_color(true)

@test plot_color(:red) == parse(C, :red)
@test plot_color("red") == parse(C, "red")
@test_throws ErrorException plot_color("notacolor")

@test plot_color(colorant"red") == C(1,0,0,1)

grad = cgrad()
@test typeof(grad) == ColorGradient
@test plot_color(grad) === grad

cs = plot_color(rand(10))
@test typeof(cs) == Vector{C}
@test length(cs) == 10

cs = plot_color(rand(4,4))
@test typeof(cs) == Matrix{C}
@test length(cs) == 16
@test size(cs) == (4,4)


cs = plot_color(rand(10), 0.5)
@test typeof(cs) == Vector{C}
@test length(cs) == 10
for c in cs
    @test alpha(c) == 0.5
end

cs = plot_color(rand(4,4), 0.5)
@test typeof(cs) == Matrix{C}
@test length(cs) == 16
@test size(cs) == (4,4)
for c in cs
    @test alpha(c) == 0.5
end

# ----------------------
# ticks

@test optimize_ticks(-1,2) == ([-1.0,0.0,1.0,2.0],-1.0,2.0)
