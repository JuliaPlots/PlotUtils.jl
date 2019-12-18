using PlotUtils
using Test
using Statistics: mean

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

# Copied from Plots.is_uniformly_spaced to avoid dependency on recent version
# on Plots which is not used on Travis.
function is_uniformly_spaced(v; tol=1e-6)
  dv = diff(v)
  maximum(dv) - minimum(dv) < tol * mean(abs.(dv))
end

@testset "ticks" begin
    @test optimize_ticks(-1,2) == ([-1.0,0.0,1.0,2.0],-1.0,2.0)
    
    @testset "small range $x, $(i)ϵ" for x in exp10.(-12:12), i in -5:5
        y = x + i*eps(x)
        x,y = minmax(x,y)
        ticks = PlotUtils.optimize_ticks(x, y)[1]
        @test issorted(ticks)
        @test all(x .<= ticks .<= y)
        # Fails:
        # @test allunique(ticks)
    end

    @testset "digits $((10^n)-1)*10^$i" for n in 1:9, i in -9:9
        y0 = 10^n
        x0 = y0-1
        x, y = (x0,y0) .* 10.0^i
        ticks = optimize_ticks(x, y)[1]
        @test length(ticks) >= 2
        @test issorted(ticks)
        @test all(x .<= ticks .<= y)
        @test is_uniformly_spaced(ticks)
    end
end
