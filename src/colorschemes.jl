## Type Piracy
# Until ColorSchemes 3.7 is released is released
Base.getindex(cs::ColorScheme, i::AbstractFloat) = get(cs, i)
Base.getindex(cs::ColorScheme, i::AbstractVector{<: AbstractFloat}) = get(cs, i)
Base.lastindex(cs::ColorScheme) = lastindex(cs.colors)


## ColorGradient

struct ColorGradient
    colors::ColorScheme
end

Base.show(io::IO, m::MIME"image/svg+xml", cg::ColorGradient) = show(io, m, cg.colors)
Base.length(cg::ColorGradient) = length(cg.colors)
Base.getindex(cg::ColorGradient, x) = getindex(cg.colors, x)
Base.size(cg::ColorGradient) = size(cg.colors)
Base.IndexStyle(::Type{<:ColorGradient}) = IndexStyle(cg.colors)
Base.iterate(cg::ColorGradient) = iterate(cg.colors)
Base.iterate(cg::ColorGradient, s) = iterate(cg.colors, s)
Base.reverse(cg::ColorGradient) = ColorGradient(reverse(cg.colors))
Base.get(cg::ColorGradient, args...) = get(cg.colors, args...)
Base.lastindex(cg::ColorGradient) = lastindex(cg.colors)

"""
    cgrad(cs; scale = identity, rev = false, alpha = nothing)

Construct a `ColorGradient`. Accepts symbols for Colorschemes.jl `ColorScheme`s,
`ColorScheme`s, vectors of colors, `ColorGradient`s and `ColorPalettes`.
If `rev` is `true` colors are reversed. `scale` accepts symbols the `:log`, `:log10`,
`:log2`, `:ln`, `:exp`, `:exp10` or functions. If `alpha` is set, it is applied to all
colors.
"""
function cgrad(cs; scale = identity, rev = false, alpha = nothing)
    if cs === :default
        cs = :inferno
    end
    cs = get_colorscheme(cs)
    if alpha !== nothing
        rgbs = convert.(RGB, cs.colors)
        cs = ColorScheme(RGBA.(rgbs, alpha))
    end
    if rev
        cs = reverse(cs)
    end
    sf = scale_function(scale)
    if scale !== identity
        n = length(cs)
        cs = ColorScheme(get(cs, sf.(range(1/n, stop = 1, length = n)), :extrema))
    end
    return ColorGradient(cs)
end

const SCALE_FUNCTIONS = Dict{Symbol, Function}(
    :log => log10,
    :log10 => log10,
    :log2 => log2,
    :ln => log,
    :exp => exp10,
    :exp10 => exp10,
)
scale_function(f) = f
scale_function(sym::Symbol) = get(SCALE_FUNCTIONS, sym, identity)

cgrad(; kw...) = cgrad(DEFAULT_COLOR_GRADIENT[]; kw...)

function default_cgrad(cg; kw...)
    DEFAULT_COLOR_GRADIENT[] = cgrad(cg; kw...)
end


## ColorPalette

struct ColorPalette
    colors::ColorScheme
end

Base.show(io::IO, m::MIME"image/svg+xml", cp::ColorPalette) = show(io, m, cp.colors)
Base.length(cg::ColorPalette) = length(cg.colors)
Base.getindex(cg::ColorPalette, x) = getindex(cg.colors, x)
Base.size(cg::ColorPalette) = size(cg.colors)
Base.IndexStyle(::Type{<:ColorPalette}) = IndexStyle(cg.colors)
Base.iterate(cg::ColorPalette) = iterate(cg.colors)
Base.iterate(cg::ColorPalette, s) = iterate(cg.colors, s)
Base.reverse(cg::ColorPalette) = ColorGradient(reverse(cg.colors))
function Base.get(cg::ColorPalette, x, rangescale)
    rangescale == :clamp && (rangescale = (0.0, 1.0))
    rangescale == :extrema && (rangescale = extrema(x))
    (rangescale isa NTuple{2, Number}) || error("rangescale ($rangescale) not supported, should be :clamp, :extrema or tuple (minVal, maxVal).  Got $(rangescale).")
   x isa AbstractRange && (x = collect(x))
   x = clamp.(x, rangescale...)
   index = Int.(round.(ColorSchemes.remap(x, rangescale..., 1, length(cscheme))))
   return cg[index]
end
Base.lastindex(cp::ColorPalette) = lastindex(cp.colors)

"""
    palette(cs; scale = identity, rev = false, alpha = nothing)

Construct a `ColorGradient`. Accepts symbols for Colorschemes.jl `ColorScheme`s,
`ColorScheme`s, vectors of colors, `ColorGradient`s and `ColorPalettes`.
If `rev` is `true` colors are reversed. `scale` accepts symbols the `:log`, `:log10`,
`:log2`, `:ln`, `:exp`, `:exp10` or functions. If `alpha` is set, it is applied to all
colors.
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


## Utils

get_colorscheme(v) = ColorScheme(v)
function get_colorscheme(sym::Symbol)
    if sym === :default
        generate_colorscheme()
    elseif haskey(ColorSchemes.colorschemes, sym)
        ColorSchemes.colorschemes[sym]
    else
        error("Unknown ColorScheme `:$sym`. Check https://juliagraphics.github.io/ColorSchemes.jl/stable/ for available ColorSchemes.")
    end
end
get_colorscheme(cs::ColorScheme) = cs
get_colorscheme(cg::ColorGradient) = cg.colors
get_colorscheme(cp::ColorPalette) = cp.colors


function cvec(cs, n = 10; kw...)
    cg = cgrad(cg; kw...)
    return RGBA{Float64}[cg[z] for z in range(0; stop=1, length=n)]
end


get_color_palette(v, n) = palette(v)
get_color_palette(cg::ColorGradient, n) = palette(cg[get_zvalues(n)])


# allows passing a ColorGradient to rgba_string and get a useful response by picking the first color - introduced because the plotly backend to Plots uses this functionality
rgba_string(cg::T) where T <: Union{ColorScheme, ColorGradient, ColorPalette} = rgba_string(cg[1])


is_colorscheme(sym) = sym in keys(ColorSchemes.colorschemes)


const DEFAULT_COLOR_GRADIENT = Ref(cgrad(ColorSchemes.colorschemes[:inferno]))
