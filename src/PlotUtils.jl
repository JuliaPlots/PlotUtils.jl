
__precompile__()

module PlotUtils

using Reexport
@reexport using Colors

export
    ColorGradient,
    cgrad,
    cvec,
    rgb_string,
    rgba_string,
    invisible,
    get_color_palette,
    isdark,
    plot_color

include("color_utils.jl")
include("color_gradients.jl")
include("colors.jl")

export
    optimize_ticks

include("ticks.jl")

end # module
