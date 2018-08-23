
# NOTE: This file was moved from Gadfly.jl, and the original author is Daniel Jones (@dcjones)

# Find the smallest order of magnitude that is larger than xspan This is a
# little opaque because I want to avoid assuming the log function is defined
# over typeof(xspan)
function bounding_order_of_magnitude(xspan::DT) where DT
    one_dt = convert(DT, one(DT))

    a = 1
    step = 1
    while xspan < 10.0^a * one_dt
        a -= step
    end

    b = 1
    step = 1
    while xspan > 10.0^b * one_dt
        b += step
    end

    while a + 1 < b
        c = div(a + b, 2)
        if xspan < 10.0^c * one_dt
            b = c
        else
            a = c
        end
    end

    return b
end


# Empty catchall
optimize_ticks() = Any[]


# Find some reasonable values for tick marks.
#
# This is basically Wilkinson's ad-hoc scoring method that tries to balance
# tight fit around the data, optimal number of ticks, and simple numbers.
#
# Args:
#   x_min: minimum value occuring in the data.
#   x_max: maximum value occuring in the data.
#   Q: tick intervals and scores
#   k_min: minimum number of ticks
#   k_max: maximum number of ticks
#   k_ideal: ideal number of ticks
#   strict_span: true if no ticks should be outside [x_min, x_max].
#
# Returns:
#   A Float64 vector containing tick marks.
#
function optimize_ticks(x_min::T, x_max::T; extend_ticks::Bool=false,
                           Q=[(1.0,1.0), (5.0, 0.9), (2.0, 0.7), (2.5, 0.5), (3.0, 0.2)],
                           k_min::Int=2, k_max::Int=10, k_ideal::Int=5,
                           granularity_weight::Float64=1/4, simplicity_weight::Float64=1/6,
                           coverage_weight::Float64=1/3, niceness_weight::Float64=1/4,
                           strict_span=true, span_buffer = nothing) where T

    Qv = [(Float64(q[1]), Float64(q[2])) for q in Q]
    optimize_ticks_typed(x_min, x_max, extend_ticks, Qv, k_min, k_max, k_ideal,
                         granularity_weight, simplicity_weight,
                         coverage_weight, niceness_weight, strict_span, span_buffer)
end

function optimize_ticks_typed(x_min::T, x_max::T, extend_ticks,
                           Q::Vector{Tuple{Float64,Float64}}, k_min,
                           k_max, k_ideal,
                           granularity_weight::Float64, simplicity_weight::Float64,
                           coverage_weight::Float64, niceness_weight::Float64,
                           strict_span, span_buffer) where T
    one_t = convert(T, one(T))
    if x_max - x_min < eps()*one_t
        R = typeof(1.0 * one_t)
        return R[x_min], x_min - one_t, x_min + one_t
    end

    n = length(Q)

    # generalizing "order of magnitude"
    xspan = x_max - x_min
    z = bounding_order_of_magnitude(xspan)

    high_score = -Inf
    S_best = Array{typeof(1.0 * one_t)}(undef, )
    viewmin_best, viewmax_best = x_min, x_max

    while 2k_max * 10.0^(z+1) * one_t > xspan
        for k in k_min:2k_max
            for (q, qscore) in Q
                tickspan = q * 10.0^z * one_t
                span = (k - 1) * tickspan
                if span < xspan
                    continue
                end

                stp = q*10.0^z
                if stp < eps()
                    continue
                end
                r = ceil((x_max - span) / (stp * one_t))

                while r*stp * one_t <= x_min
                    # Filter or expand ticks
                    if extend_ticks
                        S = Array{typeof(1.0 * one_t)}(undef, Int(3 * k))
                        for i in 0:(3*k - 1)
                            S[i+1] = (r + i - k) * tickspan
                        end
                        viewmin, viewmax = S[k + 1], S[2 * k]
                    else
                        S = Array{typeof(1.0 * one_t)}(undef, k)
                        for i in 0:(k - 1)
                            S[i+1] = (r + i) * tickspan
                        end
                        viewmin, viewmax = S[1], S[end]
                    end
                    if strict_span
                        viewmin = max(viewmin, x_min)
                        viewmax = min(viewmax, x_max)
                        if span_buffer == nothing
                            S = filter(si -> viewmin <= si <= viewmax, S)
                        else
                            buf = span_buffer * (viewmax - viewmin)
                            # @show buf
                            S = filter(si -> viewmin-buf <= si <= viewmax+buf, S)
                        end
                    end

                    # evaluate quality of ticks

                    has_zero = r <= 0 && abs(r) < k

                    # simplicity
                    s = has_zero ? 1.0 : 0.0

                    # granularity
                    g = 0 < length(S) < 2k_ideal ? 1 - abs(length(S) - k_ideal) / k_ideal : 0.0

                    # coverage
                    effective_span = (length(S)-1) * tickspan
                    c = 1.5 * xspan/effective_span

                    score = granularity_weight * g +
                            simplicity_weight * s +
                            coverage_weight * c +
                            niceness_weight * qscore

                    # strict limits on coverage
                    if strict_span && span > xspan
                        score -= 10000
                    end
                    if  span >= 2.0*xspan
                        score -= 1000
                    end

                    if score > high_score && (k_min <= length(S) <= k_max)
                        (S_best, viewmin_best, viewmax_best) = (S, viewmin, viewmax)
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
            @warn("No strict ticks found")
            return optimize_ticks_typed(x_min, x_max, extend_ticks,
                                       Q, k_min,
                                       k_max, k_ideal,
                                       granularity_weight, simplicity_weight,
                                       coverage_weight, niceness_weight,
                                       false, span_buffer)
        else
            R = typeof(1.0 * one_t)
            return R[x_min], x_min - one_t, x_min + one_t
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
                                      k_min=max(1, round(Int, 2*scale)),
                                      k_max=max(3, round(Int, 10*scale)),
                                      k_ideal=max(2, round(Int, 15*scale)))[1]
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
function optimize_datetime_ticks(a_min::Real, a_max::Real;
    k_min::Integer = 2,
    k_max::Integer = 4)

    x_min = DateTime(Dates.UTM(Int(round(a_min))))
    x_max = DateTime(Dates.UTM(Int(round(a_max))))

    Δt = x_max - x_min
    if Δt > Dates.Day(365 * k_min)
        P = Dates.Year
        steplength = Δt / (k_max * Dates.Millisecond(Dates.Day(365)))
    elseif Δt > Dates.Day(30 * k_min)
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
    steplength = P(max(1, Integer(round(steplength))))

    period_hierarchy = [Dates.Millisecond, Dates.Second, Dates.Minute,
        Dates.Hour, Dates.Day, Dates.Month, Dates.Year]
    i = findfirst(period_hierarchy .== P)
    showtype = i >= 5 ? Date : DateTime
    start = DateTime((PH(x_min) for PH in period_hierarchy[end:-1:i])...) + P(1)
    ticks = collect(start:steplength:x_max)
    labels = string.(showtype.(ticks))

    return Dates.value.(ticks), labels
end
