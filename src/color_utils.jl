
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
		hchoices=range(0; stop=340, length=20)
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

get_color_palette(palette::Vector{C}, bgcolor::Colorant, numcolors::Integer) where {C<:Colorant} = palette


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
        append!(zvalues, offset .+ [0.0, 0.5, 0.25, 0.75])
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

make255(x) = round(Int, 255 * x)

function rgb_string(c::Colorant)
  @sprintf("rgb(%d, %d, %d)", make255(red(c)), make255(green(c)), make255(blue(c)))
end

function rgba_string(c::Colorant)
	@sprintf("rgba(%d, %d, %d, %1.3f)", make255(red(c)), make255(green(c)), make255(blue(c)), alpha(c))
end
