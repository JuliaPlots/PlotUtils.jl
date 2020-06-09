
__precompile__()

module PlotUtils

using ColorSchemes
using Dates
using Reexport
using Printf
@reexport using Colors
import Base: getindex
import Random: MersenneTwister

export
    ColorGradient,
    ColorPalette,
    cgrad,
    palette,
    color_list,
    cvec,
    rgb_string,
    rgba_string,
    invisible,
    get_color_palette,
    isdark,
    plot_color,
    adapted_grid,
    default_cgrad,
    zscale

include("color_utils.jl")
include("colors.jl")
include("colorschemes.jl")
include("adapted_grid.jl")
include("intervals.jl")

export
    optimize_ticks,
    optimize_datetime_ticks

include("ticks.jl")

const _default_colorscheme = generate_colorscheme()

end # module
