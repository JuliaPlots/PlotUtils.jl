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
    rgba_string,
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
