
threshold(gray) = RidlerCalvard(gray)
threshold(polim::PolarImage) = threshold(pixels(polim))

##
## Ridler Calvard method
# adopted with permission from Jason Merrill for MIT license
function RidlerCalvard(gray)
    max_iter = 100
    tol = 1e-6

    thresh = mean(gray)
    count = 1

    # Single pass over input for both high and low mean is
    # almost 5 times faster because no temporary array allocation.
    function highlowmean(gray, thresh)
        lowmean = highmean = zero(thresh)
        lowcount = highcount = 0
        for pixel in gray            
            if pixel < thresh
                lowmean += pixel
                lowcount += 1
            else
                highmean += pixel
                highcount += 1
            end
        end
        highmean/highcount, lowmean/lowcount
    end

    while count < max_iter        
        high, low  = highlowmean(gray, thresh)
        thresh_old = thresh
        thresh = (high + low)/2
        abs(thresh - thresh_old) < tol && break
        count += 1
    end
    return(thresh)
end

# Specialized function for Ufixed for fast threshold comparison
# FixedPointNumbers.Ufixed has value field i with integer value
function edgefoptim{T<:FixedPointNumbers.Ufixed}(im::Array{T}, th)    
    t = convert(T, th)
    s = StreamStats.Mean()
    
    # work with dequeues for speed, using push! and shift!
    imq = im[:,1]
    thq = Bool[k.i < t.i for k in imq]
        
    for j = 2:size(im,2)
        # p stands for pixel, b for boolean
        prevp = im[1,j]
        prevb = prevp.i < t.i
        upp = shift!(imq)        
        upb = shift!(thq)
        for i = 2:size(im,1)
            p = im[i,j]
            b = p.i < t.i
            leftp = shift!(imq)
            leftb = shift!(thq)
            
            b != upb && StreamStats.update!(s, abs(p - upp))
            upb != prevb && StreamStats.update!(s, abs(upp - prevp))
            upb != leftb && StreamStats.update!(s, abs(upp - leftp))
            leftb != prevb && StreamStats.update!(s, abs(leftp - prevp))
            
            push!(imq, prevp)
            push!(thq, prevb)
            prevp, prevb = p, b
            upp, upb = leftp, leftb
        end
        push!(imq, prevp)
        push!(thq, prevb)
    end
    StreamStats.state(s)/2
end

##
## Edge detection method
# General version of the optimization function for edge detection
function edgefoptim(im, th)    
    t = convert(eltype(im), th)
    s = StreamStats.Mean()
    
    # work with dequeues for speed, using push! and shift!
    imq = im[:,1]
    thq = im[:,1] .< t
        
    for j = 2:size(im,2)
        # p stands for pixel, b for boolean; 
        #`prev` is one up, `up` is upper left, `left` is one left.
        prevp = im[1,j]
        prevb = prevp < t
        upp = shift!(imq)        
        upb = shift!(thq)
        for i = 2:size(im,1)
            p = im[i,j]
            b = p < t
            leftp = shift!(imq)
            leftb = shift!(thq)
            
            b != upb && StreamStats.update!(s, abs(p - upp))
            upb != prevb && StreamStats.update!(s, abs(upp - prevp))
            upb != leftb && StreamStats.update!(s, abs(upp - leftp))
            leftb != prevb && StreamStats.update!(s, abs(leftp - prevp))
            
            push!(imq, prevp)
            push!(thq, prevb)
            prevp, prevb = p, b
            upp, upb = leftp, leftb
        end
        push!(imq, prevp)
        push!(thq, prevb)
    end
    StreamStats.state(s)/2
end

function edge_threshold(gray)
    # maximization so use negative sign
    res = Optim.optimize(x->-edgefoptim(gray, x), 0., 1.)
    res.minimum
end

# TODO cut out box around fθρ(π/2) for polarimage argument
#function edge_threshold(polim::PolarImage)

##
## Minimum method
# a fast histogram method derived from julialang PR #8952 
function fasthist(img::AbstractVector, edg::Range)    
    n = length(edg) - 1
    histcount = zeros(Int, n)
    step(edg) <= 0 && error("step(edg) must be positive")
    for pixel in img
        f = (pixel - first(edg))/step(edg)
        if 0 < f <= n
            histcount[iceil(f)] += 1
        end
    end
    edg, histcount
end
# check if histogram is bimodal by detection change in direction
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
    isbimodal(hc) || error("input is not yet bimodal")
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
# See paper Glasbey 1993 Analysis of Histogram-Based Thresholding Algorithms.
# It does not work for integer, because isbimodal can get stuck in instability.
# Smooths the vector according to $y_i = (y_{i-1} + y_i + y_{i+1})/3$.
function smooth_hist!{T<:FloatingPoint}(hc::AbstractArray{T})
    hc[end] = (hc[end-1] + hc[end]) / 3
    prev = hc[1]
    c = hc[2]
    #assume hc[0] = 0 as in paper
    hc[1] = (prev + c)/3
    @inbounds for i = 2:length(hc)-1        
        nxt = hc[i+1]
        hc[i] = (prev + c + nxt)/3
        prev = c 
        c = nxt
    end
end
# precision does not seem to increase by increasing binsize
function minimum_threshold(img; bins=256, maxiter=10_000)       
    hrange, counts = fasthist(reshape(img,length(img)), -1/(bins-1):1/(bins-1):1)
    counts = float(counts)
    cnt = 0
    while cnt < maxiter
        cnt += 1
        smooth_hist!(counts)
        isbimodal(counts) && break
    end     
    cnt == maxiter && warn("maximum iteration reached at $cnt in minimum_threshold")
    th = bimodalmin(counts)/bins
end
minimum_threshold(polim::PolarImage) = minimum_threshold(pixels(polim))