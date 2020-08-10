const RIDLER_CALVARD_MAX_ITER = 100
const RIDLER_CALVARD_TOL = 1e-6

abstract type ThresholdMethod end
struct RidlerCalvard    <: ThresholdMethod end
struct EdgeDetection    <: ThresholdMethod end
struct MinimumThreshold <: ThresholdMethod end

threshold(polim::PolarImage) = threshold(polim, RidlerCalvard())
threshold(polim::PolarImage, ::RidlerCalvard) = threshold(pixels(polim), RidlerCalvard())

# Ridler Calvard method
# ---------------------

# adopted with permission from Jason Merrill for MIT license
function threshold(gray::AbstractArray, ::RidlerCalvard)    
    
    # Single pass over input for both high and low mean is
    # almost 5 times faster because no temporary array allocation.
    function highlowmean(gray, thresh)
        lowmean  = zero(thresh)
        highmean = zero(thresh)
        lowcount = 0
        highcount= 0
        for pixel in gray            
            if pixel < thresh
                lowmean += pixel
                lowcount += 1
            else
                highmean += pixel
                highcount += 1
            end
        end
        return (highmean / highcount, lowmean / lowcount)
    end

    thresh = mean(gray)
    count = 1

    while count < RIDLER_CALVARD_MAX_ITER
        high, low  = highlowmean(gray, thresh)
        thresh_old = thresh
        thresh = (high + low) / 2
        abs(thresh - thresh_old) < RIDLER_CALVARD_TOL && break
        count += 1
    end
    return(thresh)
end


# Edge detection method
# ---------------------

"Specialized type for a fast circular queue."
mutable struct circqueue{T}
    array::Vector{T}
    len::Int
    index::Int #to be replaced
end
circqueue(A::AbstractVector) = circqueue{eltype(A)}(A, length(A), 1)
function pushshift!(f::circqueue{T}, input::T) where T
    ind = f.index
    len = f.len
    output = f.array[ind]
    f.array[ind] = input
    f.index = ifelse(ind == f.len, 1, ind + 1) #ifelse faster than ternary
    output
end

"Fitness function to be optimized for Edge Detection method."
function edgefoptim(im, th)
    t = convert(eltype(im), th)
    s = StreamMean()
    s = update(s, 0) # in case no edges

    # We cache the first row and then
    # work with circular dequeues for speed, using pushshift!
    imq = circqueue(im[:,1])
    thq = circqueue(im[:,1] .< t)

    for j = 2:size(im,2)
        # p stands for pixel, b for boolean;
        #`prev` is one up, `up` is upper left, `left` is one left.
        prevp = im[1,j]
        prevb = prevp < t
        upp = pushshift!(imq, prevp)
        upb = pushshift!(thq, prevb)
        for i = 2:size(im,1)
            p = im[i,j]
            b = p < t
            leftp = pushshift!(imq, p)
            leftb = pushshift!(thq, b)

            if b != upb; s = update(s, abs(p - upp)); end
            if upb != prevb; s = update(s, abs(upp - prevp)); end
            if upb != leftb; s = update(s, abs(upp - leftp)); end
            if leftb != prevb; s = update(s, abs(leftp - prevp)); end

            prevp, prevb = p, b
            upp, upb = leftp, leftb
        end
    end
    # The published algorithm uses:
    #    out = mean(s)/2 
    # but we multiply the contrast mean with the sqrt (because 2D) 
    # of edges count, to avoid spurious results with high threshold
    # values and very little edges.
    out = mean(s) * sqrt(length(s))
end

function threshold(gray::AbstractArray, ::EdgeDetection)
    # maximization so use negative sign
    res = optimize(x -> -edgefoptim(gray, x), 0.01, 0.99)
    minimizer(res)
end

# Cut out box around fθρ(π/2) to reduce for polarimage argument
"""
    threshold(polar_image, EdgeDetection())

Edge Detection method for automatic thresholding after Nobis & Hunziker, 2005.

It optimizes the threshold to find the maximum value of the contrast mean at edges.
The method has been slightly adapted from the original by normalizing the contrast mean with the
sqrt (because 2D) of edges count, to avoid spurious results with high threshold
values and very little edges.
"""
function threshold(polim::PolarImage, ::EdgeDetection)
    rows, cols = size(polim)
    Rmax   = ceil(Int, polim.cl.fθρ(π/2))
    ci     = polim.cl.ci
    cj     = polim.cl.cj
    rowmin = max(1, ci - Rmax)
    rowmax = min(rows, ci + Rmax)
    colmin = max(1, cj - Rmax)
    colmax = min(cols, cj + Rmax)
    threshold(polim.img[rowmin:rowmax, colmin:colmax], EdgeDetection())
end


# Minimum method
# --------------

# Check if histogram is bimodal by detecting change in direction.
# TODO rewrite with diff
function isbimodal(hc) #hc=histcounts
    prevup = hc[2] > hc[1]
    prevc = hc[2]
    mode = 0
    for i = 3:length(hc)
        c = hc[i]
        down = c < prevc

        if prevup & down
            mode += 1
        end

        if mode > 2
            return false
        end

        if !down
            prevup = ifelse(c == prevc, prevup, !down)
        else
            prevup = !down
        end
        prevc = c
    end
    return true
end
# find the minimim between modes if bimodal
function bimodalmin(hc)
    @assert isbimodal(hc)
    prevup = hc[2] > hc[1]
    prevc = hc[2]
    mode = 0
    for i = 3:length(hc)
        c = hc[i]
        down = c < prevc

        if prevup && down
            mode += 1
        end

        if (mode == 1) && (c > prevc)
            return i - 1
        end

        if !down
            prevup = ifelse(c == prevc, prevup, !down)
        else
            prevup = !down
        end
        prevc = c
    end
    return true
end

"""
See paper Glasbey 1993 Analysis of Histogram-Based Thresholding Algorithms.
It does not work for integer, because isbimodal can get stuck in instability.
Smooths the vector according to ``y_i = (y_{i-1} + y_i + y_{i+1})/3```.
"""
function smooth_hist!(hc::AbstractArray{<:AbstractFloat})
    hc[end] = (hc[end-1] + hc[end]) / 3
    prev = hc[1]
    c = hc[2]
    #assume hc[0] = 0 as in paper
    hc[1] = (prev + c)/3
    @inbounds for i = 2:(length(hc) - 1)
        nxt = hc[i+1]
        hc[i] = (prev + c + nxt) / 3
        prev = c
        c = nxt
    end
end
# precision does not seem to increase by increasing binsize
function threshold(img, ::MinimumThreshold; bins=256, maxiter=10_000)
    hist_range = -1 / (bins-1) : 1 / (bins-1) : 1
    counts = fasthist(reshape(img, length(img)), hist_range)
    counts = float(counts) # float required for convergence!
    cnt = 0
    while cnt < maxiter
        cnt += 1
        smooth_hist!(counts)
        isbimodal(counts) && break
    end
    cnt == maxiter && @warn("maximum iteration reached at $cnt in minimum_threshold")
    th = bimodalmin(counts)/bins
end
threshold(polim::PolarImage, ::MinimumThreshold) = threshold(pixels(polim), MinimumThreshold())
