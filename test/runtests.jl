using PlotUtils
using Test

# TODO: real tests

# ----------------------
# colors

const C = RGBA{Float64}

@testset "colors" begin

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

    grad = cgrad([:red, "blue"])
    @test grad.colors == C[colorant"red", colorant"blue"]
    @test grad.values == collect(range(0,stop=1,length=2))

    grad = cgrad([:red, "blue"], alpha = 0.5)
    @test grad.colors == C[C(1,0,0,0.5), C(0,0,1,0.5)]
    @test grad.values == collect(range(0,stop=1,length=2))

    grad = cgrad([:red,:blue], [0,0.1,1])
    @test grad.colors == C[C(1,0,0), C(0.5,0,0.5), C(0,0,1)]
    @test grad.values == [0,0.1,1]

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
end

# ----------------------
# gradients

@testset "gradients" begin
    grad = cgrad(:inferno)
    @test length(grad.colors) == 30
    @test length(grad.values) == 30
    @test RGB(grad.colors[1]) == RGB(0.001462, 0.000466, 0.013866)
    @test RGB(grad.colors[end]) == RGB(0.988362, 0.998364, 0.644924)
end

# ----------------------
# ticks

@testset "ticks" begin
    @test optimize_ticks(-1,2) == ([-1.0,0.0,1.0,2.0],-1.0,2.0)
end
