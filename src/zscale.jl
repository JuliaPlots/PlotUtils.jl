# ----------------------------------------------------------------------------------

# implementation of the "zscale" IRAF function for finding appropriate
# color limits in an iterative manner.

import Polynomials
using Statistics: std, median

function zscale(input,
    nsamples::Int=1000;
    contrast=0.25,
    max_reject=0.5,
    min_npixels=5,
    kref=2.5,
    max_iterations=5)

    values = filter(isfinite, input)
    stride = Int(max(1, length(values) / nsamples))
    samples = values[1:stride:end][1:nsamples]
    sort!(samples)

    N = length(samples)
    vmin = first(samples)
    vmax = last(samples)

    # fit a line to the sorted samples
    min_pix = max(min_npixels, Int(N * max_reject))
    x = 1:N-1

    ngood = N
    last_good = N + 1

    # bad pixel mask (array bool faster than bitarray)
    badmask = zeros(Bool, N)

    # kernel to dilate bad pixel mask
    ngrow = max(1, Int(N / 100))
    kernel = ones(Bool, ngrow)

    # iteratively fit samples and reject sigma-clipped outliers
    for _ in 1:max_iterations
        (ngood ≥ last_good || ngood < min_pix) && break
        
        # linear fit using mask for weighting
        fit = Polynomials.fit(x, samples, 1; weights=.!badmask)
        flat = @. samples - fit(x)

        # k-sigma rejection threshold
        threshold = krej * std(flat[.!badmask])

        # detect and reject outliers based on threshold
        @. badmask[!(-threshold ≤ flat ≤ threshold)] = true

        # enlarge mask
        # TODO

        last_good = ngood
        ngood = sum(.!badmask)
    end

    if ngood ≥ min_pix
        slope = contrast > 0 ? fit[1] / contrast : fit[1]
        center = N ÷ 2
        m = median(samples)
        vmin = max(vmin, m - (center - 1) * slope)
        vmax = min(vmax, m + (N - center) * slope)
    end

    return vmin, vmax
end
