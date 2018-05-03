
__precompile__()

module PlotUtils

using Reexport
@reexport using Colors
using Dates
import Base: getindex

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
    adapted_grid,
    clibrary,
    clibraries,
    cgradients,
    default_cgrad

include("color_utils.jl")
include("color_gradients.jl")
include("colors.jl")
include("adapted_grid.jl")

export
    optimize_ticks,
    optimize_datetime_ticks

include("ticks.jl")

end # module
