
Ufixed = FixedPointNumbers.Ufixed

gapfraction(pixs::AbstractArray, thresh) = mean(pixs .> thresh)

# specialized gapfraction method for Ufixed pixels, with threshold already converted to Ufixed type
function gapfraction{T<:Ufixed}(px::AbstractArray{T}, thresh::T)    
    gfs = 0.        
    for i = 1:length(px)
        gfs += px[i].i > thresh.i
    end
    return(gfs / length(px))
end

function contactfreq(polim::PolarImage, θ, thresh; ringwidth=RING_WIDTH)
    θ < 0   && throw(DomainError())
    θ > π/2 && throw(DomainError())
    pixs = pixels(polim, max(0, θ-ringwidth), min(pi/2, θ+ringwidth))
    P = gapfraction(pixs, convert(eltype(pixs),thresh))
    return(-log(P) * cos(θ))
end