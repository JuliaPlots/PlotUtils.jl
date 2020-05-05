#= 
This file contains code for calculating intervals for plotting purposes.
These functions should, at minimum, take in some form of data input and
return a tuple (min, max) of the limits corresponding to the function 
=#

import Polynomials
using Statistics: std, median
using ImageFiltering

"""
    zscale(input::AbstractArray, 
        nsamples=600;
        contrast=0.25,
        max_reject=0.5,
        min_npixels=5,
        k_rej=2.5,
        max_iterations=5)

Implementation of the `zscale` IRAF function for finding colorbar limits of `input`.

## Arguments
* `nsamples` - The number of samples to use from `input`. If fewer than `nsamples` are present, will use the full input
* `contrast` - The desired contrast
* `k_rej` - The number of standard deviations above which data is rejected
* `max_iteration` - The number of iterations used for fitting samples
* `max_reject` - The maximum number of pixels to reject during the iterative fitting
* `min_npixels` - The minimum number of pixels to calculate the limits after the iterative fitting

# Examples
```jldoctest
julia> img = 0:9999

julia> zscale(img)
(0, 9990)
```
"""
function zscale(input::AbstractArray,
    nsamples::Int = 1000;
    contrast = 0.25,
    max_reject = 0.5,
    min_npixels = 5,
    k_rej = 2.5,
    max_iterations = 5)

    # get samples from finite values of input
    values = float(filter(isfinite, input))
    stride = max(1, round(Int, length(values) / nsamples))
    samples = values[1:stride:end][1:min(nsamples, end)]
    sort!(samples)

    N = length(samples)
    vmin = first(samples)
    vmax = last(samples)

    # fit a line to the sorted samples
    min_pix = max(min_npixels, Int(N * max_reject))
    x = 0:N - 1

    ngood = N
    last_good = N + 1

    # bad pixel mask (array bool faster than bitarray)
    badmask = zeros(Bool, N)

    # get number of neighbors to mask if a pixel is bad
    ngrow = max(1, round(Int, N / 100))
    kernel = centered(ones(Bool, ngrow))
    nleft = ngrow ÷ 2
    nright = ngrow - nleft - 1

    local fit
    # iteratively fit samples and reject sigma-clipped outliers
    for _ in 1:max_iterations
        (ngood ≥ last_good || ngood < min_pix) && break
        
        # linear fit using mask for weighting
        fit = Polynomials.fit(x, samples, 1; weights = .!badmask)
        flat = @. samples - fit(x)

        # k-sigma rejection threshold
        threshold = k_rej * std(flat[.!badmask])

        # detect and reject outliers based on threshold
        @. badmask[!(-threshold ≤ flat ≤ threshold)] = true

        # dialate mask
        badmask .= imfilter(badmask, kernel, Fill(false))[axes(badmask)...] .> 0

        last_good = ngood
        ngood = sum(.!badmask)
    end

    if ngood ≥ min_pix
        slope = contrast > 0 ? fit[1] / contrast : fit[1]
        center = (N - 1) ÷ 2
        m = median(samples)
        vmin = max(vmin, m - (center - 1) * slope)
        vmax = min(vmax, m + (N - center) * slope)
    end

    return vmin, vmax
end
