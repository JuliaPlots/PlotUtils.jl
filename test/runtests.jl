using PlotUtils
using Base.Test

# TODO: real tests

@test plot_color(nothing) == RGBA{Float64}(0,0,0,0)

@test optimize_ticks(-1,2) == ([-1.0,0.0,1.0,2.0],-1.0,2.0)
