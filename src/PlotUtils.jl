
__precompile__()

module PlotUtils

using Reexport
@reexport using Colors

import Base: getindex, setindex!

export
    ColorGradient,
    cgrad,
    cvec,
    rgb_string,
    rgba_string,
    invisible,
    get_color_palette,
    isdark,
    plot_color,
    set_color_library,
    list_color_libraries,
    list_gradients

include("color_utils.jl")
include("color_gradients.jl")
include("colors.jl")

export
    optimize_ticks

include("ticks.jl")

end # module
