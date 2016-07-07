module PlotUtils

using Reexport
@reexport using Colors

export
    # colorscheme,
    # ColorScheme,
    ColorGradient,
    # ColorVector,
    # ColorWrapper,
    # ColorFunction,
    # ColorZFunction,
    # getColor,
    # getColorZ,
    cgrad,
    plot_color

include("color_gradients.jl")
include("color_utils.jl")
include("colors.jl")

export
    optimize_ticks

include("ticks.jl")

end # module
