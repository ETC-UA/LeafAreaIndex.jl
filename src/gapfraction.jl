
Ufixed = FixedPointNumbers.Ufixed

#fallback
gapfraction(pixs, thresh) = mean(pixs .> thresh)

#in general specialize on type of input array
gapfraction(pixs::AbstractArray, thresh)
    threshT = convert(eltype(pixs), thresh)
    gfs = zero(Int)        
    for pix in pixs
        gfs += pix > thresh
    end
    return(gfs / length(pixs))
end
# specialized gapfraction method for Ufixed pixels, with threshold already 
# converted to Ufixed type for type stability.
function gapfraction{T<:Ufixed}(pixs::AbstractArray{T}, thresh)
    threshT = convert(T, thresh)
    gfs = zero(Int)        
    for pix in pixs
        gfs += pix.i > thresh.i
    end
    return(gfs / length(pixs))
end

function contactfreq(polim::PolarImage, θ, thresh; ringwidth=RING_WIDTH)
    θ < 0   && throw(DomainError())
    θ > π/2 && throw(DomainError())
    pixs = pixels(polim, max(0, θ-ringwidth), min(pi/2, θ+ringwidth))
    P = gapfraction(pixs, thresh)
    return(-log(P) * cos(θ))
end