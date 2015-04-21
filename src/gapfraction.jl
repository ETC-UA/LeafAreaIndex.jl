
Ufixed = FixedPointNumbers.Ufixed

#fallback
gapfraction(pixs, thresh) = mean(pixs .> thresh)

#in general specialize on type of input array
function gapfraction(pixs::AbstractArray, thresh)
    threshT = convert(eltype(pixs), thresh)
    gfs = zero(Int)        
    for pix in pixs
        gfs += pix > threshT
    end
    return(gfs / length(pixs))
end
# specialized gapfraction method for Ufixed pixels, with threshold already 
# converted to Ufixed type for type stability.
function gapfraction{T<:Ufixed}(pixs::AbstractArray{T}, thresh)
    threshT = convert(T, thresh)
    gfs = zero(Int)        
    for pix in pixs
        gfs += pix.i > threshT.i
    end
    return(gfs / length(pixs))
end

function loggapfraction(pixs, thresh)
    gf = gapfraction(pixs, thresh)
    gf == zero(gf) && return(log(1/length(pixs)))
    log(gf)
end
