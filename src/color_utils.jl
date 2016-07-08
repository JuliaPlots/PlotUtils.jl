
function interpolate_rgb(c1::Colorant, c2::Colorant, w::Real)
  rgb1 = RGBA(c1)
  rgb2 = RGBA(c2)
  r = interpolate(rgb1.r, rgb2.r, w)
  g = interpolate(rgb1.g, rgb2.g, w)
  b = interpolate(rgb1.b, rgb2.b, w)
  a = interpolate(rgb1.alpha, rgb2.alpha, w)
  RGBA(r, g, b, a)
end


function interpolate(v1::Real, v2::Real, w::Real)
  (1-w) * v1 + w * v2
end


# --------------------------------------------------------------

# Methods to automatically generate gradients for color selection based on
# background color and a short list of seed colors

# here are some magic constants that could be changed if you really want
const _lightness_darkbg = [80.0]
const _lightness_lightbg = [60.0]
const _lch_c_const = [60]

function adjust_lch(color, l, c)
    lch = convert(LCHab, color)
    convert(RGBA{Float64}, LCHab(l, c, lch.h))
end

function lightness_from_background(bgcolor)
	bglight = convert(LCHab, bgcolor).l
	bglight < 50.0 ? _lightness_darkbg[1] : _lightness_lightbg[1]
end

function gradient_from_list(cs)
    zvalues = get_zvalues(length(cs))
    indices = sortperm(zvalues)
    sorted_colors = map(RGBA, cs[indices])
    sorted_zvalues = zvalues[indices]
    ColorGradient(sorted_colors, sorted_zvalues)
end

function generate_colorgradient(bgcolor = plot_color(:white);
                               color_bases = plot_color([colorant"steelblue",colorant"orangered"]),
                               lightness = lightness_from_background(bgcolor),
                               chroma = _lch_c_const[1],
                               n = 17)
    seed_colors = vcat(bgcolor, map(c -> adjust_lch(c, lightness, chroma), color_bases))
	seed_colors = convert(Vector{RGB{Float64}}, seed_colors)
	colors = distinguishable_colors(
		n,
		seed_colors,
		lchoices=Float64[lightness],
		cchoices=Float64[chroma],
		hchoices=linspace(0, 340, 20)
	)[2:end]
	gradient_from_list(colors)
end

function get_color_palette(palette, bgcolor::Colorant, numcolors::Integer)
	grad = if palette == :auto
		generate_colorgradient(bgcolor)
	else
		cgrad(palette)
	end
	zrng = get_zvalues(numcolors)
	RGBA{Float64}[grad[z] for z in zrng]
end

get_color_palette{C<:Colorant}(palette::Vector{C}, bgcolor::Colorant, numcolors::Integer) = palette


# ----------------------------------------------------------------------------------


function getpctrange(n::Int)
    n > 0 || error()
    n == 1 && return zeros(1)
    zs = [0.0, 1.0]
    for i in 3:n
        sorted = sort(zs)
        diffs = diff(sorted)
        widestj = 0
        widest = 0.0
        for (j,d) in enumerate(diffs)
            if d > widest
                widest = d
                widestj = j
            end
        end
        push!(zs, sorted[widestj] + 0.5 * diffs[widestj])
    end
    zs
end

function get_zvalues(n::Int)
    offsets = getpctrange(ceil(Int,n/4)+1)/4
    offsets = vcat(offsets[1], offsets[3:end])
    zvalues = Float64[]
    for offset in offsets
        append!(zvalues, offset + [0.0, 0.5, 0.25, 0.75])
    end
    vcat(zvalues[1], 1.0, zvalues[2:n-1])
end


# ----------------------------------------------------------------------------------

function darken(c, v=0.1)
    rgba = convert(RGBA, c)
    r = max(0, min(rgba.r - v, 1))
    g = max(0, min(rgba.g - v, 1))
    b = max(0, min(rgba.b - v, 1))
    RGBA(r,g,b,rgba.alpha)
end
function lighten(c, v=0.3)
    darken(c, -v)
end


# isbackgrounddark(bgcolor::Color) = Lab(bgcolor).l < 0.5

# move closer to lighter/darker depending on background value
function adjust_away(val, bgval, vmin=0., vmax=100.)
  if bgval < 0.5 * (vmax+vmin)
    tmp = max(val, bgval)
    return 0.5 * (tmp + max(tmp, vmax))
  else
    tmp = min(val, bgval)
    return 0.5 * (tmp + min(tmp, vmin))
  end
end

# borrowed from http://stackoverflow.com/a/1855903:
lightness_level(c::Colorant) = 0.299 * red(c) + 0.587 * green(c) + 0.114 * blue(c)

isdark(c::Colorant) = lightness_level(c) < 0.5
islight(c::Colorant) = !isdark(c)


function Base.convert(::Type{RGB}, h::Unsigned)
  mask = 0x0000FF
  RGB([(x & mask) / 0xFF for x in  (h >> 16, h >> 8, h)]...)
end


# # note: I found this list of hex values in a comment by Tatarize here: http://stackoverflow.com/a/12224359
# const _masterColorList = [
#     0xFFFFFF, 0x000000, 0x0000FF, 0x00FF00, 0xFF0000, 0x01FFFE, 0xFFA6FE, 0xFFDB66, 0x006401, 0x010067,
#     0x95003A, 0x007DB5, 0xFF00F6, 0xFFEEE8, 0x774D00, 0x90FB92, 0x0076FF, 0xD5FF00, 0xFF937E, 0x6A826C,
#     0xFF029D, 0xFE8900, 0x7A4782, 0x7E2DD2, 0x85A900, 0xFF0056, 0xA42400, 0x00AE7E, 0x683D3B, 0xBDC6FF,
#     0x263400, 0xBDD393, 0x00B917, 0x9E008E, 0x001544, 0xC28C9F, 0xFF74A3, 0x01D0FF, 0x004754, 0xE56FFE,
#     0x788231, 0x0E4CA1, 0x91D0CB, 0xBE9970, 0x968AE8, 0xBB8800, 0x43002C, 0xDEFF74, 0x00FFC6, 0xFFE502,
#     0x620E00, 0x008F9C, 0x98FF52, 0x7544B1, 0xB500FF, 0x00FF78, 0xFF6E41, 0x005F39, 0x6B6882, 0x5FAD4E,
#     0xA75740, 0xA5FFD2, 0xFFB167, 0x009BFF, 0xE85EBE
#   ]
# const _allColors = map(RGB, _masterColorList)
# const _darkColors = filter(isdark, _allColors)
# const _lightColors = filter(islight, _allColors)
# const _sortedColorsForDarkBackground = vcat(_lightColors, reverse(_darkColors[2:end]))
# const _sortedColorsForLightBackground = vcat(_darkColors, reverse(_lightColors[2:end]))


make255(x) = round(Int, 255 * x)

function rgba_string(c::Colorant)
	@sprintf("rgba(%d, %d, %d, %1.3f)", make255(red(c)), make255(green(c)), make255(blue(c)), alpha(c))
end

# function webcolor(c::Color)
#     @sprintf("rgb(%d, %d, %d)", [make255(f(c)) for f in [red,green,blue]]...)
# end
# function webcolor(c::TransparentColor)
#     @sprintf("rgba(%d, %d, %d, %1.3f)", [make255(f(c)) for f in [red,green,blue]]..., alpha(c))
# end
# # webcolor(cs::ColorScheme) = webcolor(getColor(cs))
# # webcolor(c) = webcolor(convertColor(c))
# webcolor(c) = webcolor(plot_color(c))
# webcolor(c, α) = webcolor(convertColor(getColor(c), α))
