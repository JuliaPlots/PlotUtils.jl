
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
    clibrary,
    clibraries,
    cgradients,
    setdefaults,
    adapted_grid


include("color_utils.jl")
include("color_gradients.jl")
include("colors.jl")
include("adapted_grid.jl")

export
    optimize_ticks

include("ticks.jl")

end # module
