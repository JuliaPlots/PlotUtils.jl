# An abstract color list is either a Palette or a gradient
abstract type AbstractColorList end

get_colorscheme(acl::AbstractColorList) = acl.colors
color_list(acl::AbstractColorList) = color_list(acl.colors)
plot_color(cl::AbstractColorList) = cl

# Interface
Base.show(io::IO, m::MIME"image/svg+xml", acl::AbstractColorList) = show(io, m, color_list(acl))
Base.length(acl::AbstractColorList) = length(color_list(acl))
Base.getindex(acl::AbstractColorList, x) = getindex(color_list(acl), x)
Base.size(acl::AbstractColorList) = size(color_list(acl))
Base.IndexStyle(::Type{<:AbstractColorList}) = IndexLinear()
Base.iterate(acl::AbstractColorList) = iterate(color_list(acl))
Base.iterate(acl::AbstractColorList, s) = iterate(color_list(acl), s)
Base.lastindex(acl::AbstractColorList) = lastindex(color_list(acl))
Base.firstindex(acl::AbstractColorList) = firstindex(color_list(acl))
Base.get(acl::AbstractColorList, args...) = get(get_colorscheme(acl), args...)


## ColorGradient

abstract type ColorGradient <: AbstractColorList end

Base.getindex(cg::ColorGradient, x::Union{AbstractFloat, AbstractVector{<:AbstractFloat}}) =
    get(cg, x)
function Base.get(cg::ColorGradient, v::AbstractArray, rangescale = (0.0, 1.0))
    rangescale == :extrema && (rangescale = extrema(v))
    map(x -> get(cg, x, rangescale), v)
end

Base.:(==)(cg1::ColorGradient, cg2::ColorGradient) = color_list(cg1) == color_list(cg2) && cg1.values == cg2.values
Base.hash(cg::ColorGradient) = hash(color_list(cg)) ⊻ hash(cg.values)

## Continuous Color Gradient

struct ContinuousColorGradient <: ColorGradient
    colors::ColorScheme
    values::Vector{Float64}

    function ContinuousColorGradient(colors, values = get_range(colors))
        c, v = prepare_continuous_cgrad_colors(colors, values)
        new(c, v)
    end
end

plot_color(cg::ContinuousColorGradient, α::Number) =
    ContinuousColorGradient(plot_color(cg.colors, α), cg.values)

Base.show(io::IO, m::MIME"image/svg+xml", cg::ContinuousColorGradient) =
    show(io, m, cg[get_range(100)])

function Base.get(cg::ContinuousColorGradient, x::AbstractFloat, rangescale = (0.0, 1.0))
    isfinite(x) || return invisible()
    rangescale = get_rangescale(rangescale)
    allunique(rangescale) || return cg.colors[1]
    x = clamp(x, rangescale...)
    if rangescale != (0.0, 1.0)
        x = ColorSchemes.remap(x, rangescale..., 0, 1)
    end
    c, v = cg.colors, cg.values
    if x in v
        index = findfirst(==(x), v)
    else
        xx = x # avoid boxing
        i = findlast(t -> t < xx, v)
        r = (x - v[i]) / (v[i + 1] - v[i])
        index = (i + r - 1) / (length(v) - 1)
    end
    return c[index]
end

Base.reverse(cg::ContinuousColorGradient) =
    ContinuousColorGradient(reverse(cg.colors), reverse(1 .- cg.values))

# required for GR's nonuniformcellarray
function ColorSchemes.getinverse(cg::ContinuousColorGradient, c)
    alpha(c) == 0 && return NaN
    z = getinverse(to_rgb(get_colorscheme(cg)), to_rgb(c))
    cr = get_range(cg.colors)
    if z in cr
        i = findfirst(==(z), cr)
        return cg.values[i]
    else
        i = findlast(x -> x < z, cr)
        ColorSchemes.remap(z, cr[i], cr[i + 1], cg.values[i], cg.values[i + 1])
    end
end

function prepare_continuous_cgrad_colors(c, v)
    v = sort(unique(clamp.(v, 0, 1)))
    0 in v || pushfirst!(v, 0)
    1 in v || push!(v, 1)
    nv = length(v)
    nc = length(c)
    c, v = if nc != nv
        value_range = get_range(nv)
        color_range = get_range(nc)
        values = [0.0]
        for i in 2:nv
            inds = findall(x -> value_range[i - 1] < x < value_range[i], color_range)
            if !isempty(inds)
                append!(
                    values,
                    ColorSchemes.remap(
                        color_range[inds],
                        value_range[i - 1],
                        value_range[i],
                        v[i - 1],
                        v[i],
                    ),
                )
            end
            push!(values, v[i])
        end
        colors = c[sort(unique([value_range; color_range]))]
        colors, values
    else
        color_list(c), v
    end
    return ColorScheme(plot_color(c)), v
end


## Categorical Color Gradient

struct CategoricalColorGradient <: ColorGradient
    colors::ColorScheme
    values::Vector{Float64}

    function CategoricalColorGradient(colors, values = get_range(colors))
        c, v = prepare_categorical_cgrad_colors(colors, values)
        new(c, v)
    end
end

plot_color(cg::CategoricalColorGradient, α::Number) =
    CategoricalColorGradient(plot_color(cg.colors, α), cg.values)

Base.show(io::IO, m::MIME"image/svg+xml", cg::CategoricalColorGradient) =
    show(io, m, cg[get_range(100)])

function Base.get(cg::CategoricalColorGradient, x::AbstractFloat, rangescale = (0.0, 1.0))
    isfinite(x) || return invisible()
    rangescale = get_rangescale(rangescale)
    x = clamp.(x, rangescale...)
    if rangescale != (0.0, 1.0)
        x = ColorSchemes.remap(x, rangescale..., 0, 1)
    end
    return cg.colors[x == 0 ? 1 : findlast(z -> z < x, cg.values)]
end

Base.reverse(cg::CategoricalColorGradient) =
    CategoricalColorGradient(reverse(cg.colors), reverse(1 .- cg.values))

# required for GR's nonuniformcellarray
function ColorSchemes.getinverse(cg::CategoricalColorGradient, c)
    alpha(c) == 0 && return NaN
    i = findfirst(==(c), color_list(cg))
    return (cg.values[i] + cg.values[i + 1]) / 2
end

function prepare_categorical_cgrad_colors(c, v)
    v = sort(unique(clamp.(v, 0, 1)))
    0 in v || pushfirst!(v, 0)
    1 in v || push!(v, 1)
    c = ColorScheme(plot_color(c[get_range(length(v) - 1)]))
    return c, v
end


## cgrad

"""
    cgrad(colors, [values]; categorical = nothing, scale = nothing, rev = false, alpha = nothing)

Construct a Colorgradient from `colors` and `values`.

`colors` can be a symbol for Colorschemes.jl `ColorScheme`s, a `ColorScheme`, a vector of colors, a `ColorGradient` or a `ColorPalette`.
If `values` is an integer, it specifies the numbers of colors chosen equidistantly from the colorscheme specified by colors.
Otherwise vectors are accepted.
For continuous color gradients `values` indicate where between 0 and 1 the colors are positioned.
For categorical color gradients `values` indicate where a color ends and where a new one begins between 0 and 1.
0 and 1 are added to `values` if not already present.

If `rev` is `true` colors are reversed.
`scale` accepts the symbols `:log`, `:log10`, `:log2`, `:ln`, `:exp`, `:exp10` or functions.
If `alpha` is set, it is applied to all colors.
"""
function cgrad(
    colors::ColorScheme,
    values;
    categorical = nothing,
    scale = nothing,
    rev = false,
    alpha = nothing,
)
    if categorical !== nothing
        colors, values = prepare_categorical_cgrad_colors(colors, values)
    end

    if alpha !== nothing
        rgbs = convert.(RGB, colors.colors)
        colors = ColorScheme(RGBA.(rgbs, alpha))
    end
    if rev
        colors = reverse(colors)
    end
    values = if scale in (:log, :log10) || scale isa typeof(log10)
        log10.(ColorSchemes.remap(values, 0, 1, 1, 10))
    elseif scale == :log2 || scale isa typeof(log2)
        log2.(ColorSchemes.remap(values, 0, 1, 1, 2))
    elseif scale == :ln || scale isa typeof(log)
        log.(ColorSchemes.remap(values, 0, 1, 1, ℯ))
    elseif scale in (:exp, :exp10) || scale isa typeof(exp10) || scale isa typeof(exp)
        ColorSchemes.remap(exp10.(values), 1, 10, 0, 1)
    elseif scale isa Function
        v = scale.(values)
        ColorSchemes.remap(v, extrema(v)..., 0, 1)
    else
        values
    end

    if categorical !== nothing
        return CategoricalColorGradient(colors, values)
    else
        return ContinuousColorGradient(colors, values)
    end
end

function cgrad(colors::ColorScheme, n::Int=length(colors); categorical = nothing, kwargs...)
    values = get_range(n + (categorical !== nothing))
    return cgrad(colors, values; categorical = categorical, kwargs...)
end

function cgrad(colors, args...; kwargs...)
    if colors === :default
        colors = :inferno
    end
    return cgrad(get_colorscheme(colors), args...; kwargs...)
end

cgrad(; kw...) = cgrad(DEFAULT_COLOR_GRADIENT[]; kw...)

function default_cgrad(cg; kw...)
    DEFAULT_COLOR_GRADIENT[] = cgrad(cg; kw...)
end

function get_rangescale(rangescale)
    rangescale == :clamp && return (0.0, 1.0)
    rangescale == :extrema && return extrema(x)
    (rangescale isa NTuple{2, Number}) ||
    error("rangescale ($rangescale) not supported, should be :clamp, :extrema or tuple (minVal, maxVal).  Got $(rangescale).")
    return rangescale
end


## ColorPalette

struct ColorPalette <: AbstractColorList
    colors::ColorScheme
end

plot_color(cp::ColorPalette, α::Number) = palette(cp, alpha = α)

Base.reverse(cp::ColorPalette) = ColorPalette(reverse(cp.colors))

"""
    palette(colors, [n]; rev = false, alpha = nothing)

Construct a `ColorPalette`.
Accepts symbols for Colorschemes.jl `ColorScheme`s, `ColorScheme`s, vectors of colors, `ColorGradient`s and `ColorPalettes`.
If 'n' is provided, `n` colors are chosen equidistantly from the colorscheme.
If `rev` is `true` colors are reversed.
"""
function palette(cs; rev = false, alpha = nothing)
    cs = get_colorscheme(cs)
    if alpha !== nothing
        rgbs = convert.(RGB, cs.colors)
        cs = ColorScheme(RGBA.(rgbs, alpha))
    end
    if rev
        cs = reverse(cs)
    end
    return ColorPalette(cs)
end

function palette(cs, n; kwargs...)
    cs = get_colorscheme(cs)
    z = get_range(n)
    return palette(cs[z]; kwargs...)
end

## Utils

get_range(n::Int) = range(0, stop = 1, length = n)
get_range(cs) = get_range(length(cs))

get_colorscheme(v::AbstractVector{<:Colorant}) = ColorScheme(v)
get_colorscheme(v::AbstractVector) = ColorScheme(parse.(Colorant, v))
function get_colorscheme(sym::Symbol)
    if haskey(MISC_COLORSCHEMES, sym)
        return MISC_COLORSCHEMES[sym]
    end
    sym = get(COLORSCHEME_ALIASES, sym, sym)
    if sym === :default || sym === :auto
        _default_colorscheme
    elseif haskey(ColorSchemes.colorschemes, sym)
        ColorSchemes.colorschemes[sym]
    else
        error("Unknown ColorScheme `:$sym`. Check https://juliagraphics.github.io/ColorSchemes.jl/stable/ for available ColorSchemes.")
    end
end
get_colorscheme(cs::ColorScheme) = cs


function cvec(cs, n = 10; kw...)
    cg = cgrad(cs; kw...)
    return RGBA{Float64}[cg[z] for z in get_range(n)]
end


color_list(c) = get_colorscheme(c).colors
color_list(v::AbstractVector) = plot_color.(v)

get_color_palette(v, n) = palette(v)
get_color_palette(cg::ColorGradient, n) = palette(cg[get_zvalues(n)])

to_rgb(c::Colorant) = RGB(c)
to_rgb(cs::ColorScheme) = ColorScheme(to_rgb.(cs.colors))
to_rgb(cg::ColorGradient) = ColorGradient(to_rgb(cg.colors))
to_rgb(cp::ColorPalette) = ColorPalette(to_rgb(cp.colors))


# allows passing a ColorGradient to rgba_string and get a useful response by picking the first color - introduced because the plotly backend to Plots uses this functionality
rgba_string(cg::T) where {T <: Union{ColorScheme, ColorGradient, ColorPalette}} =
    rgba_string(cg[1])


is_colorscheme(sym) =
    sym in keys(ColorSchemes.colorschemes) ||
    sym in keys(COLORSCHEME_ALIASES) || sym in keys(MISC_COLORSCHEMES)


const DEFAULT_COLOR_GRADIENT = Ref(cgrad(ColorSchemes.colorschemes[:inferno]))


## Compat

const COLORSCHEME_ALIASES = Dict{Symbol, Symbol}(
    # colorbrewer
    :Spectral => :Spectral_11,
    :RdYlGn => :RdYlGn_11,
    :RdBu => :RdBu_11,
    :PiYG => :PiYG_11,
    :PRGn => :PRGn_11,
    :RdYlBu => :RdYlBu_11,
    :BrBG => :BrBG_11,
    :RdGy => :RdGy_11,
    :Set2 => :Set2_8,
    :Accent => :Accent_8,
    :Set1 => :Set1_9,
    :Set3 => :Set3_12,
    :Dark2 => :Dark2_8,
    :Paired => :Paired_12,
    :Pastel2 => :Pastel2_8,
    :Pastel1 => :Pastel1_9,
    :PuOr => :PuOr_11,
    :OrRd => :OrRd_9,
    :PuBu => :PuBu_9,
    :BuPu => :BuPu_9,
    :Oranges => :Oranges_9,
    :BuGn => :BuGn_9,
    :YlOrBr => :YlOrBr_9,
    :YlGn => :YlGn_9,
    :Reds => :Reds_9,
    :RdPu => :RdPu_9,
    :Greens => :Greens_9,
    :YlGnBu => :YlGnBu_9,
    :Purples => :Purples_9,
    :GnBu => :GnBu_9,
    :Greys => :Greys_9,
    :YlOrRd => :YlOrRd_8,
    :PuRd => :PuRd_9,
    :Blues => :Blues_9,
    :PuBuGn => :PuBuGn_9,
    # colorcet
    :bgy => :linear_bgy_10_95_c74_n256,
    :bgyw => :linear_bgyw_15_100_c68_n256,
    :bjy => :diverging_linear_bjy_30_90_c45_n256,
    :bkr => :diverging_bkr_55_10_c35_n256,
    :bky => :diverging_bky_60_10_c30_n256,
    :blues => :linear_blue_5_95_c73_n256,
    :bmw => :linear_bmw_5_95_c86_n256,
    :colorwheel => :cyclic_mygbm_30_95_c78_n256_s25,
    :dimgreys => :linear_grey_10_95_c0_n256,
    :fire => :linear_kryw_0_100_c71_n256,
    :greys => :linear_grey_0_100_c0_n256,
    :gwv => :diverging_gwv_55_95_c39_n256,
    :cinferno => :linear_bmy_10_95_c78_n256,
    :isolum => :isoluminant_cgo_80_c38_n256,
    :kb => :linear_ternary_blue_0_44_c57_n256,
    :kdc => :linear_blue_5_95_c73_n256,
    :kg => :linear_ternary_green_0_46_c42_n256,
    :kgy => :linear_green_5_95_c69_n256,
    :kr => :linear_ternary_red_0_50_c52_n256,
    :diverging => :diverging_bwr_40_95_c42_n256,
    :cyclic1 => :cyclic_mrybm_35_75_c68_n256,
    :cyclic2 => :cyclic_mygbm_30_95_c78_n256,
    :cyclic3 => :cyclic_wrwbw_40_90_c42_n256,
)

## Plots misc colorschemes

const RAINBOW_COLORS = RGBA{Float64}[
    colorant"purple",
    colorant"blue",
    colorant"green",
    colorant"orange",
    colorant"red",
]
const TEST_COLORS = RGBA{Float64}[
    colorant"darkblue",
    colorant"blueviolet",
    colorant"darkcyan",
    colorant"green",
    darken(colorant"yellow", 0.3),
    colorant"orange",
    darken(colorant"red", 0.2),
]

const MISC_COLORSCHEMES = Dict{Symbol, ColorScheme}(
    :blues => ColorScheme(RGBA{Float64}[colorant"lightblue", colorant"darkblue"]),
    :reds => ColorScheme(RGBA{Float64}[colorant"lightpink", colorant"darkred"]),
    :greens => ColorScheme(RGBA{Float64}[colorant"lightgreen", colorant"darkgreen"]),
    :redsblues => ColorScheme(RGBA{Float64}[
        colorant"darkred",
        RGB(0.8, 0.85, 0.8),
        colorant"darkblue",
    ]),
    :bluesreds => ColorScheme(RGBA{Float64}[
        colorant"darkblue",
        RGB(0.8, 0.85, 0.8),
        colorant"darkred",
    ]),
    :heat => ColorScheme(RGBA{Float64}[
        colorant"lightyellow",
        colorant"orange",
        colorant"darkred",
    ]),
    :grays => ColorScheme(RGBA{Float64}[RGB(0.05, 0.05, 0.05), RGB(0.95, 0.95, 0.95)]),
    :rainbow => ColorScheme(RAINBOW_COLORS),
    :lightrainbow => ColorScheme(map(lighten, RAINBOW_COLORS)),
    :darkrainbow => ColorScheme(map(darken, RAINBOW_COLORS)),
    :darktest => ColorScheme(TEST_COLORS),
    :lighttest => ColorScheme(map(c -> lighten(c, 0.3), TEST_COLORS)),
    :cividis => ColorScheme([
        RGB(0.0000, 0.1262, 0.3015),
        RGB(0.0000, 0.1292, 0.3077),
        RGB(0.0000, 0.1321, 0.3142),
        RGB(0.0000, 0.1350, 0.3205),
        RGB(0.0000, 0.1379, 0.3269),
        RGB(0.0000, 0.1408, 0.3334),
        RGB(0.0000, 0.1437, 0.3400),
        RGB(0.0000, 0.1465, 0.3467),
        RGB(0.0000, 0.1492, 0.3537),
        RGB(0.0000, 0.1519, 0.3606),
        RGB(0.0000, 0.1546, 0.3676),
        RGB(0.0000, 0.1574, 0.3746),
        RGB(0.0000, 0.1601, 0.3817),
        RGB(0.0000, 0.1629, 0.3888),
        RGB(0.0000, 0.1657, 0.3960),
        RGB(0.0000, 0.1685, 0.4031),
        RGB(0.0000, 0.1714, 0.4102),
        RGB(0.0000, 0.1743, 0.4172),
        RGB(0.0000, 0.1773, 0.4241),
        RGB(0.0000, 0.1798, 0.4307),
        RGB(0.0000, 0.1817, 0.4347),
        RGB(0.0000, 0.1834, 0.4363),
        RGB(0.0000, 0.1852, 0.4368),
        RGB(0.0000, 0.1872, 0.4368),
        RGB(0.0000, 0.1901, 0.4365),
        RGB(0.0000, 0.1930, 0.4361),
        RGB(0.0000, 0.1958, 0.4356),
        RGB(0.0000, 0.1987, 0.4349),
        RGB(0.0000, 0.2015, 0.4343),
        RGB(0.0000, 0.2044, 0.4336),
        RGB(0.0000, 0.2073, 0.4329),
        RGB(0.0055, 0.2101, 0.4322),
        RGB(0.0236, 0.2130, 0.4314),
        RGB(0.0416, 0.2158, 0.4308),
        RGB(0.0576, 0.2187, 0.4301),
        RGB(0.0710, 0.2215, 0.4293),
        RGB(0.0827, 0.2244, 0.4287),
        RGB(0.0932, 0.2272, 0.4280),
        RGB(0.1030, 0.2300, 0.4274),
        RGB(0.1120, 0.2329, 0.4268),
        RGB(0.1204, 0.2357, 0.4262),
        RGB(0.1283, 0.2385, 0.4256),
        RGB(0.1359, 0.2414, 0.4251),
        RGB(0.1431, 0.2442, 0.4245),
        RGB(0.1500, 0.2470, 0.4241),
        RGB(0.1566, 0.2498, 0.4236),
        RGB(0.1630, 0.2526, 0.4232),
        RGB(0.1692, 0.2555, 0.4228),
        RGB(0.1752, 0.2583, 0.4224),
        RGB(0.1811, 0.2611, 0.4220),
        RGB(0.1868, 0.2639, 0.4217),
        RGB(0.1923, 0.2667, 0.4214),
        RGB(0.1977, 0.2695, 0.4212),
        RGB(0.2030, 0.2723, 0.4209),
        RGB(0.2082, 0.2751, 0.4207),
        RGB(0.2133, 0.2780, 0.4205),
        RGB(0.2183, 0.2808, 0.4204),
        RGB(0.2232, 0.2836, 0.4203),
        RGB(0.2281, 0.2864, 0.4202),
        RGB(0.2328, 0.2892, 0.4201),
        RGB(0.2375, 0.2920, 0.4200),
        RGB(0.2421, 0.2948, 0.4200),
        RGB(0.2466, 0.2976, 0.4200),
        RGB(0.2511, 0.3004, 0.4201),
        RGB(0.2556, 0.3032, 0.4201),
        RGB(0.2599, 0.3060, 0.4202),
        RGB(0.2643, 0.3088, 0.4203),
        RGB(0.2686, 0.3116, 0.4205),
        RGB(0.2728, 0.3144, 0.4206),
        RGB(0.2770, 0.3172, 0.4208),
        RGB(0.2811, 0.3200, 0.4210),
        RGB(0.2853, 0.3228, 0.4212),
        RGB(0.2894, 0.3256, 0.4215),
        RGB(0.2934, 0.3284, 0.4218),
        RGB(0.2974, 0.3312, 0.4221),
        RGB(0.3014, 0.3340, 0.4224),
        RGB(0.3054, 0.3368, 0.4227),
        RGB(0.3093, 0.3396, 0.4231),
        RGB(0.3132, 0.3424, 0.4236),
        RGB(0.3170, 0.3453, 0.4240),
        RGB(0.3209, 0.3481, 0.4244),
        RGB(0.3247, 0.3509, 0.4249),
        RGB(0.3285, 0.3537, 0.4254),
        RGB(0.3323, 0.3565, 0.4259),
        RGB(0.3361, 0.3593, 0.4264),
        RGB(0.3398, 0.3622, 0.4270),
        RGB(0.3435, 0.3650, 0.4276),
        RGB(0.3472, 0.3678, 0.4282),
        RGB(0.3509, 0.3706, 0.4288),
        RGB(0.3546, 0.3734, 0.4294),
        RGB(0.3582, 0.3763, 0.4302),
        RGB(0.3619, 0.3791, 0.4308),
        RGB(0.3655, 0.3819, 0.4316),
        RGB(0.3691, 0.3848, 0.4322),
        RGB(0.3727, 0.3876, 0.4331),
        RGB(0.3763, 0.3904, 0.4338),
        RGB(0.3798, 0.3933, 0.4346),
        RGB(0.3834, 0.3961, 0.4355),
        RGB(0.3869, 0.3990, 0.4364),
        RGB(0.3905, 0.4018, 0.4372),
        RGB(0.3940, 0.4047, 0.4381),
        RGB(0.3975, 0.4075, 0.4390),
        RGB(0.4010, 0.4104, 0.4400),
        RGB(0.4045, 0.4132, 0.4409),
        RGB(0.4080, 0.4161, 0.4419),
        RGB(0.4114, 0.4189, 0.4430),
        RGB(0.4149, 0.4218, 0.4440),
        RGB(0.4183, 0.4247, 0.4450),
        RGB(0.4218, 0.4275, 0.4462),
        RGB(0.4252, 0.4304, 0.4473),
        RGB(0.4286, 0.4333, 0.4485),
        RGB(0.4320, 0.4362, 0.4496),
        RGB(0.4354, 0.4390, 0.4508),
        RGB(0.4388, 0.4419, 0.4521),
        RGB(0.4422, 0.4448, 0.4534),
        RGB(0.4456, 0.4477, 0.4547),
        RGB(0.4489, 0.4506, 0.4561),
        RGB(0.4523, 0.4535, 0.4575),
        RGB(0.4556, 0.4564, 0.4589),
        RGB(0.4589, 0.4593, 0.4604),
        RGB(0.4622, 0.4622, 0.4620),
        RGB(0.4656, 0.4651, 0.4635),
        RGB(0.4689, 0.4680, 0.4650),
        RGB(0.4722, 0.4709, 0.4665),
        RGB(0.4756, 0.4738, 0.4679),
        RGB(0.4790, 0.4767, 0.4691),
        RGB(0.4825, 0.4797, 0.4701),
        RGB(0.4861, 0.4826, 0.4707),
        RGB(0.4897, 0.4856, 0.4714),
        RGB(0.4934, 0.4886, 0.4719),
        RGB(0.4971, 0.4915, 0.4723),
        RGB(0.5008, 0.4945, 0.4727),
        RGB(0.5045, 0.4975, 0.4730),
        RGB(0.5083, 0.5005, 0.4732),
        RGB(0.5121, 0.5035, 0.4734),
        RGB(0.5158, 0.5065, 0.4736),
        RGB(0.5196, 0.5095, 0.4737),
        RGB(0.5234, 0.5125, 0.4738),
        RGB(0.5272, 0.5155, 0.4739),
        RGB(0.5310, 0.5186, 0.4739),
        RGB(0.5349, 0.5216, 0.4738),
        RGB(0.5387, 0.5246, 0.4739),
        RGB(0.5425, 0.5277, 0.4738),
        RGB(0.5464, 0.5307, 0.4736),
        RGB(0.5502, 0.5338, 0.4735),
        RGB(0.5541, 0.5368, 0.4733),
        RGB(0.5579, 0.5399, 0.4732),
        RGB(0.5618, 0.5430, 0.4729),
        RGB(0.5657, 0.5461, 0.4727),
        RGB(0.5696, 0.5491, 0.4723),
        RGB(0.5735, 0.5522, 0.4720),
        RGB(0.5774, 0.5553, 0.4717),
        RGB(0.5813, 0.5584, 0.4714),
        RGB(0.5852, 0.5615, 0.4709),
        RGB(0.5892, 0.5646, 0.4705),
        RGB(0.5931, 0.5678, 0.4701),
        RGB(0.5970, 0.5709, 0.4696),
        RGB(0.6010, 0.5740, 0.4691),
        RGB(0.6050, 0.5772, 0.4685),
        RGB(0.6089, 0.5803, 0.4680),
        RGB(0.6129, 0.5835, 0.4673),
        RGB(0.6168, 0.5866, 0.4668),
        RGB(0.6208, 0.5898, 0.4662),
        RGB(0.6248, 0.5929, 0.4655),
        RGB(0.6288, 0.5961, 0.4649),
        RGB(0.6328, 0.5993, 0.4641),
        RGB(0.6368, 0.6025, 0.4632),
        RGB(0.6408, 0.6057, 0.4625),
        RGB(0.6449, 0.6089, 0.4617),
        RGB(0.6489, 0.6121, 0.4609),
        RGB(0.6529, 0.6153, 0.4600),
        RGB(0.6570, 0.6185, 0.4591),
        RGB(0.6610, 0.6217, 0.4583),
        RGB(0.6651, 0.6250, 0.4573),
        RGB(0.6691, 0.6282, 0.4562),
        RGB(0.6732, 0.6315, 0.4553),
        RGB(0.6773, 0.6347, 0.4543),
        RGB(0.6813, 0.6380, 0.4532),
        RGB(0.6854, 0.6412, 0.4521),
        RGB(0.6895, 0.6445, 0.4511),
        RGB(0.6936, 0.6478, 0.4499),
        RGB(0.6977, 0.6511, 0.4487),
        RGB(0.7018, 0.6544, 0.4475),
        RGB(0.7060, 0.6577, 0.4463),
        RGB(0.7101, 0.6610, 0.4450),
        RGB(0.7142, 0.6643, 0.4437),
        RGB(0.7184, 0.6676, 0.4424),
        RGB(0.7225, 0.6710, 0.4409),
        RGB(0.7267, 0.6743, 0.4396),
        RGB(0.7308, 0.6776, 0.4382),
        RGB(0.7350, 0.6810, 0.4368),
        RGB(0.7392, 0.6844, 0.4352),
        RGB(0.7434, 0.6877, 0.4338),
        RGB(0.7476, 0.6911, 0.4322),
        RGB(0.7518, 0.6945, 0.4307),
        RGB(0.7560, 0.6979, 0.4290),
        RGB(0.7602, 0.7013, 0.4273),
        RGB(0.7644, 0.7047, 0.4258),
        RGB(0.7686, 0.7081, 0.4241),
        RGB(0.7729, 0.7115, 0.4223),
        RGB(0.7771, 0.7150, 0.4205),
        RGB(0.7814, 0.7184, 0.4188),
        RGB(0.7856, 0.7218, 0.4168),
        RGB(0.7899, 0.7253, 0.4150),
        RGB(0.7942, 0.7288, 0.4129),
        RGB(0.7985, 0.7322, 0.4111),
        RGB(0.8027, 0.7357, 0.4090),
        RGB(0.8070, 0.7392, 0.4070),
        RGB(0.8114, 0.7427, 0.4049),
        RGB(0.8157, 0.7462, 0.4028),
        RGB(0.8200, 0.7497, 0.4007),
        RGB(0.8243, 0.7532, 0.3984),
        RGB(0.8287, 0.7568, 0.3961),
        RGB(0.8330, 0.7603, 0.3938),
        RGB(0.8374, 0.7639, 0.3915),
        RGB(0.8417, 0.7674, 0.3892),
        RGB(0.8461, 0.7710, 0.3869),
        RGB(0.8505, 0.7745, 0.3843),
        RGB(0.8548, 0.7781, 0.3818),
        RGB(0.8592, 0.7817, 0.3793),
        RGB(0.8636, 0.7853, 0.3766),
        RGB(0.8681, 0.7889, 0.3739),
        RGB(0.8725, 0.7926, 0.3712),
        RGB(0.8769, 0.7962, 0.3684),
        RGB(0.8813, 0.7998, 0.3657),
        RGB(0.8858, 0.8035, 0.3627),
        RGB(0.8902, 0.8071, 0.3599),
        RGB(0.8947, 0.8108, 0.3569),
        RGB(0.8992, 0.8145, 0.3538),
        RGB(0.9037, 0.8182, 0.3507),
        RGB(0.9082, 0.8219, 0.3474),
        RGB(0.9127, 0.8256, 0.3442),
        RGB(0.9172, 0.8293, 0.3409),
        RGB(0.9217, 0.8330, 0.3374),
        RGB(0.9262, 0.8367, 0.3340),
        RGB(0.9308, 0.8405, 0.3306),
        RGB(0.9353, 0.8442, 0.3268),
        RGB(0.9399, 0.8480, 0.3232),
        RGB(0.9444, 0.8518, 0.3195),
        RGB(0.9490, 0.8556, 0.3155),
        RGB(0.9536, 0.8593, 0.3116),
        RGB(0.9582, 0.8632, 0.3076),
        RGB(0.9628, 0.8670, 0.3034),
        RGB(0.9674, 0.8708, 0.2990),
        RGB(0.9721, 0.8746, 0.2947),
        RGB(0.9767, 0.8785, 0.2901),
        RGB(0.9814, 0.8823, 0.2856),
        RGB(0.9860, 0.8862, 0.2807),
        RGB(0.9907, 0.8901, 0.2759),
        RGB(0.9954, 0.8940, 0.2708),
        RGB(1.0000, 0.8979, 0.2655),
        RGB(1.0000, 0.9018, 0.2600),
        RGB(1.0000, 0.9057, 0.2593),
        RGB(1.0000, 0.9094, 0.2634),
        RGB(1.0000, 0.9131, 0.2680),
        RGB(1.0000, 0.9169, 0.2731),
    ]),
)
