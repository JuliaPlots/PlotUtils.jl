const _logScales = [:ln, :log2, :log10]
const _logScaleBases = Dict(:ln => ℯ, :log2 => 2.0, :log10 => 10.0)

# NOTE: This file was moved from Gadfly.jl, and the original author is Daniel Jones (@dcjones)

# Find the smallest order of magnitude that is larger than xspan This is a
# little opaque because I want to avoid assuming the log function is defined
# over typeof(xspan)
function bounding_order_of_magnitude(xspan::T, base) where T
    a = 1
    step = 1
    while xspan < T(base^a)
        a -= step
    end

    b = 1
    step = 1
    while xspan > T(base^b)
        b += step
    end

    while a + 1 < b
        c = div(a + b, 2)
        if xspan < T(base^c)
            b = c
        else
            a = c
        end
    end

    return b
end

const float_digit_range = floor(Int,log10(floatmin())):ceil(Int,log10(floatmax()))
postdecimal_digits(x) = first(i for i in float_digit_range if x==floor(x; digits=i))

fallback_ticks(x_min::T, x_max::T, k_min, k_max) where T = collect(T, range(x_min, x_max; length=k_min)), x_min, x_max

# Empty catchall
optimize_ticks() = Any[]

"""
    optimize_ticks(xmin, xmax; extend_ticks::Bool=false,
                           Q=[(1.0,1.0), (5.0, 0.9), (2.0, 0.7), (2.5, 0.5), (3.0, 0.2)],
                           k_min::Int=2, k_max::Int=10, k_ideal::Int=5,
                           granularity_weight::Float64=1/4, simplicity_weight::Float64=1/6,
                           coverage_weight::Float64=1/3, niceness_weight::Float64=1/4,
                           strict_span=true, span_buffer = nothing
        )

Find some reasonable values for tick marks.

This is basically Wilkinson's ad-hoc scoring method that tries to balance
tight fit around the data, optimal number of ticks, and simple numbers.

## Arguments:

*    `xmax`:

    The maximum value occurring in the data.

*    `xmin`:

    The minimum value occurring in the data.

*    `extend_ticks`:

    Determines whether to extend tick computation. Defaults to false.

*    `strict_span`:

    True if no ticks should be outside [x_min, x_max]. Defaults to true.

*    `Q`:

    A distribution of nice numbers from which labellings are sampled. Stored in the form (number, score).

*    `k_min`:

    The minimum number of ticks.

*    `k_max`:

    The maximum number of ticks.

*    `k_ideal`:

    The ideal number of ticks.

*    `granularity_weight`:

    Encourages returning roughly the number of labels requested.

*    `simplicity_weight`:

    Encourages nicer labeling sequences by preferring step sizes that appear earlier in Q.
    Also rewards labelings that include 0 as a way to ground the sequence.

*    `coverage_weight`:

    Encourages labelings that do not extend far beyond the range of the data, penalizing unnecessary whitespace.

*    `niceness_weight`:

    Encourages labellings to produce nice ranges.

## Returns:
  `(ticklocations::Vector{Float64}, x_min, x_max)`

## Mathematical details

Wilkinson’s optimization function is defined as the sum of three
components. If the user requests m labels and a possible labeling has
k labels, then the components are `simplicity`, `coverage` and `granularity`.

These components are defined as follows:
```math
\\begin{aligned}
  &\\text{simplicity} = 1 - \\frac{i}{|Q|} + \\frac{v}{|Q|}\\\\
  &\\text{coverage}   = \\frac{x_{max} - x_{min}}{\\mathrm{label}_{max} - \\mathrm{label}_{min}}\\\\
  &\\text{granularity}= 1 - \\frac{\\left|k - m\\right|}{m}
\\end{aligned}
```

and the variables here are:

*  `q`: element of `Q`.
*  `i`: index of `q` in `Q`.
*  `v`: 1 if label range includes 0, 0 otherwise.
"""
function optimize_ticks(x_min::T, x_max::T; extend_ticks::Bool=false,
                           Q=[(1., 1.), (5., .9), (2., .7), (2.5, .5), (3., .2)],
                           k_min::Int=2, k_max::Int=10, k_ideal::Int=5,
                           granularity_weight::Float64=1/4, simplicity_weight::Float64=1/6,
                           coverage_weight::Float64=1/3, niceness_weight::Float64=1/4,
                           strict_span=true, span_buffer=nothing, scale=nothing) where T

    F = float(T)
    Qv = [(F(q[1]), F(q[2])) for q in Q]
    optimize_ticks_typed(x_min, x_max, extend_ticks, Qv, k_min, k_max, k_ideal,
                         granularity_weight, simplicity_weight,
                         coverage_weight, niceness_weight, strict_span, span_buffer, scale)
end

function optimize_ticks_typed(x_min::T, x_max::T, extend_ticks,
                           Q::Vector{Tuple{Float64,Float64}}, k_min, k_max, k_ideal,
                           granularity_weight::Float64, simplicity_weight::Float64,
                           coverage_weight::Float64, niceness_weight::Float64,
                           strict_span, span_buffer, scale) where T
    xspan = x_max - x_min
    F = float(T)
    xspan < eps(F) && return fallback_ticks(x_min, x_max, k_min, k_max)

    n = length(Q)
    is_log_scale = scale ∈ _logScales
    base = F(get(_logScaleBases, scale, 10.))

    # generalizing "order of magnitude"
    z = bounding_order_of_magnitude(xspan, base)

    # find required significant digits for ticks with q * base^z spacing,
    # for q values specified in Q
    x_digits = bounding_order_of_magnitude(max(abs(x_min), abs(x_max)), base)
    q_extra_digits = maximum(postdecimal_digits(q[1]) for q in Q)
    sigdigits(z) = max(1, x_digits - z + q_extra_digits)

    ib = Int(base)
    round_base = (
        isinteger(base) ?
        v -> round(v, sigdigits=sigdigits(z), base=ib) :
        v -> round(v, sigdigits=sigdigits(z))
    )

    high_score = -Inf
    S_best = Array{F}(undef, 1)
    viewmin_best, viewmax_best = x_min, x_max

    # we preallocate arrays that hold all required S arrays for every given
    # the k parameter, so we don't have to create them again and again, which
    # saves many allocations
    prealloc_Ss = [Array{F}(undef, extend_ticks ? Int(3k) : k) for k in k_min:2k_max]

    while 2k_max * base^(z + 1) > xspan
        for (ik, k) in enumerate(k_min:2k_max)
            for (q, qscore) in Q
                tickspan = q * base^z
                span = (k - 1) * tickspan
                span < xspan && continue

                stp = tickspan
                stp < eps(F) && continue

                has_nice_scale = true
                if is_log_scale && !isinteger(stp)
                    qscore = 0  # prefer integer exponents for log scales
                    has_nice_scale = false
                end
                r = (x_max - span) / stp
                isfinite(r) || continue
                r = ceil(Int, r)

                while r * stp <= x_min
                    # Filter or expand ticks
                    if extend_ticks
                        S = prealloc_Ss[ik]
                        for i in 0:(3k - 1)
                            S[i+1] = (r + i - k) * tickspan
                        end
                        # round only those values that end up as viewmin and viewmax
                        # to save computation time
                        S[k + 1] = round_base(S[k + 1])
                        S[2k] = round_base(S[2k])
                        viewmin, viewmax = S[k + 1], S[2k]
                    else
                        S = prealloc_Ss[ik]
                        for i in 0:(k - 1)
                            S[i+1] = (r + i) * tickspan
                        end
                        # round only those values that end up as viewmin and viewmax
                        # to save computation time
                        S[1] = round_base(S[1])
                        S[k] = round_base(S[k])
                        viewmin, viewmax = S[1], S[k]
                    end
                    if strict_span
                        viewmin = max(viewmin, x_min)
                        viewmax = min(viewmax, x_max)
                        buf = something(span_buffer, 0) * (viewmax - viewmin)

                        # filter the S array while reusing its own memory to do so
                        # this works because S is sorted, and we will only overwrite
                        # values that are not needed anymore going forward in the loop

                        # we do this because it saves allocations and leaves S type stable
                        counter = 0
                        @inbounds for i in 1:length(S)
                            if (viewmin - buf) <= S[i] <= (viewmax + buf)
                                counter += 1
                                S[counter] = S[i]
                            end
                        end
                        S = view(S, 1:counter)
                    end

                    # evaluate quality of ticks
                    has_zero = r <= 0 && abs(r) < k

                    # simplicity
                    s = has_zero && has_nice_scale ? 1 : 0

                    # granularity
                    g = 0 < length(S) < 2k_ideal ? 1 - abs(length(S) - k_ideal) / k_ideal : 0

                    # coverage
                    effective_span = (length(S) - 1) * tickspan
                    c = abs(effective_span) > eps(F) ? 1.5xspan / effective_span : 0

                    score = granularity_weight * g +
                            simplicity_weight * s +
                            coverage_weight * c +
                            niceness_weight * qscore

                    # strict limits on coverage
                    if strict_span && span > xspan
                        score -= 10000
                    end
                    if span >= 2xspan
                        score -= 1000
                    end

                    if score > high_score && (k_min <= length(S) <= k_max)
                        if strict_span
                            # make S a copy because it is a view and
                            # could otherwise be mutated in the next runs
                            S = collect(S)
                        end
                        S_best, viewmin_best, viewmax_best = S, viewmin, viewmax
                        high_score = score
                    end
                    r += 1
                end
            end
        end
        z -= 1
    end

    if isinf(high_score)
        if strict_span
            @warn "No strict ticks found"
            return optimize_ticks_typed(x_min, x_max, extend_ticks, Q, k_min, k_max, k_ideal,
                                       granularity_weight, simplicity_weight,
                                       coverage_weight, niceness_weight,
                                       false, span_buffer, scale)
        else
            return fallback_ticks(x_min, x_max, k_min, k_max)
        end
    end

    return S_best, viewmin_best, viewmax_best
end


function optimize_ticks(x_min::Date, x_max::Date; extend_ticks::Bool=false,
                        k_min=nothing, k_max=nothing, scale=:auto,
                        granularity_weight=nothing, simplicity_weight=nothing,
                        coverage_weight=nothing, niceness_weight=nothing,
                        strict_span=true, span_buffer = nothing)
    return optimize_ticks(convert(DateTime, x_min), convert(DateTime, x_max),
                          extend_ticks=extend_ticks, scale=scale)
end


function optimize_ticks(x_min::DateTime, x_max::DateTime; extend_ticks::Bool=false,
                        k_min=nothing, k_max=nothing, scale=:auto,
                        granularity_weight=nothing, simplicity_weight=nothing,
                        coverage_weight=nothing, niceness_weight=nothing,
                        strict_span=true, span_buffer = nothing)
    if x_min == x_max
        x_max += Second(1)
    end

    if year(x_max) - year(x_min) <= 1 && scale != :year
        if year(x_max) == year(x_min) && month(x_max) - month(x_min) <= 1 && scale != :month
            ticks = DateTime[]

            scales = [
                Day(1), Hour(1), Minute(1), Second(1), Millisecond(100),
                Millisecond(10), Millisecond(1)
            ]

            # ticks on week boundries
            if x_min + Day(7) < x_max || scale == :week
                push!(ticks, x_min)
                while true
                    next_month = Date(year(ticks[end]), month(ticks[end])) + Month(1)
                    while ticks[end] + Week(1) < next_month - Day(2)
                        push!(ticks, ticks[end] + Week(1))
                    end
                    push!(ticks, next_month)
                    if next_month >= x_max
                        break
                    end
                end
            else
                scale = nothing
                if scale != :auto
                    # TODO: manually setting scale with :day, :minute, etc
                end

                if scale === nothing
                    for proposed_scale in [Day(1), Hour(1), Minute(1),
                                           Second(1), Millisecond(100),
                                           Millisecond(10), Millisecond(1)]
                        if x_min + proposed_scale < x_max
                            scale = proposed_scale
                            break
                        end
                    end
                end

                if scale === nothing
                    scale = Millisecond(1)
                end

                # round x_min down
                if scale === Day(1)
                    first_tick = DateTime(year(x_min), month(x_min), day(x_min))
                elseif scale === Hour(1)
                    first_tick = DateTime(year(x_min), month(x_min), day(x_min),
                                          hour(x_min))
                elseif scale === Minute(1)
                    first_tick = DateTime(year(x_min), month(x_min), day(x_min),
                                          hour(x_min), minute(x_min))
                elseif scale === Second(1)
                    first_tick = DateTime(year(x_min), month(x_min), day(x_min),
                                          hour(x_min), minute(x_min), second(x_min))
                elseif scale === Millisecond(100)
                    first_tick = DateTime(year(x_min), month(x_min), day(x_min),
                                          hour(x_min), minute(x_min),
                                          second(x_min), millisecond(x_min) % 100)
                elseif scale === Millisecond(10)
                    first_tick = DateTime(year(x_min), month(x_min), day(x_min),
                                          hour(x_min), minute(x_min),
                                          second(x_min), millisecond(x_min) % 10)
                else
                    first_tick = x_min
                end
                push!(ticks, first_tick)

                while ticks[end] < x_max
                    push!(ticks, ticks[end] + scale)
                end
            end

            viewmin, viewmax = ticks[1], ticks[end]
            return ticks, viewmin, viewmax
        else
            ticks = DateTime[]
            push!(ticks, Date(year(x_min), month(x_min)))
            while ticks[end] < x_max
                push!(ticks, ticks[end] + Month(1))
            end
            viewmin, viewmax = ticks[1], ticks[end]

            return ticks, x_min, x_max
        end
    else
        ticks, viewmin, viewmax =
            optimize_ticks(year(x_min), year(x_max + Year(1) - Day(1)), extend_ticks=extend_ticks)

        return DateTime[DateTime(round(y)) for y in ticks],
                        DateTime(round(viewmin)), DateTime(round(viewmax))
    end
end



# Generate ticks suitable for multiple scales.
function multilevel_ticks(viewmin::T, viewmax::T;
                          scales=[0.5, 5.0, 10.0]) where T

    ticks = Dict()
    for scale in scales
        ticks[scale] = optimize_ticks(viewmin, viewmax,
                                      k_min=max(1, round(Int, 2scale)),
                                      k_max=max(3, round(Int, 10scale)),
                                      k_ideal=max(2, round(Int, 15scale)))[1]
    end

    return ticks
end


function multilevel_ticks(viewmin::Date, viewmax::Date;
                          scales=[:year, :month, :day])
    return multilevel_ticks(convert(DateTime, viewmin),
                            convert(DateTime, viewmax),
                            scales=scales)
end


function multilevel_ticks(viewmin::DateTime, viewmax::DateTime;
                          scales=[:year, :month, :day])
    # TODO: This needs to be improved for DateTime
    span = convert(Float64, Dates.toms(viewmax - viewmin))
    ticks = Dict()
    for scale in scales
        if scale == :year
            s = span / Dates.toms(Day(360))
        elseif scale == :month
            s = span / Dates.toms(Day(90))
        else
            s = span / Dates.toms(Day(1))
        end

        ticks[s/20] = optimize_ticks(viewmin, viewmax, scale=scale)[1]
    end

    return ticks
end


# Choose "round" (full seconds/minutes/hours/days/months/years) DateTime ticks
# between x_min and x_max:
function optimize_datetime_ticks(a_min, a_max; k_min = 2, k_max = 4)
    x_min = DateTime(Dates.UTM(Int(round(a_min))))
    x_max = DateTime(Dates.UTM(Int(round(a_max))))

    Δt = x_max - x_min
    if Δt > Dates.Day(365k_min)
        P = Dates.Year
        steplength = Δt / (k_max * Dates.Millisecond(Dates.Day(365)))
    elseif Δt > Dates.Day(30k_min)
        P = Dates.Month
        steplength = Δt / (k_max * Dates.Millisecond(Dates.Day(30)))
    elseif Δt > Dates.Day(k_min)
        P = Dates.Day
        steplength = Δt / (k_max * Dates.Millisecond(Dates.Day(1)))
    elseif Δt > Dates.Hour(k_min)
        P = Dates.Hour
        steplength = Δt / (k_max * Dates.Millisecond(Dates.Hour(1)))
    elseif Δt > Dates.Minute(k_min)
        P = Dates.Minute
        steplength = Δt / (k_max * Dates.Millisecond(Dates.Minute(1)))
    elseif Δt > Dates.Second(k_min)
        P = Dates.Second
        steplength = Δt / (k_max * Dates.Millisecond(Dates.Second(1)))
    else
        P = Dates.Millisecond
        steplength = Δt / (k_max * Dates.Millisecond(1))
    end
    steplength = P(max(1, Int(round(steplength))))

    period_hierarchy = [Dates.Millisecond, Dates.Second, Dates.Minute,
        Dates.Hour, Dates.Day, Dates.Month, Dates.Year]
    i = findfirst(period_hierarchy .== P)
    showtype = i >= 5 ? Date : DateTime
    start = DateTime((PH(x_min) for PH in period_hierarchy[end:-1:i])...) + P(1)
    ticks = collect(start:steplength:x_max)
    labels = string.(showtype.(ticks))

    return Dates.value.(ticks), labels
end
