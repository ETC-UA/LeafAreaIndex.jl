
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


function weightedrings(polim::PolarImage, θ1::Real, θ2::Real, N::Integer)
    
    # create edges for θ rings with similar number of pixels each
    θedges = map(polim.cl.fρθ, sqrt(linspace(polim.cl.fθρ(θ1)^2, (polim.cl.fθρ(θ2))^2, N+1)))
    #fix possible floating point roundoff
    θedges[1] = max(θedges[1], θ1) 
    θedges[end] = min(θedges[end], θ2)
    
    # weighted average midpoints
    θdouble = map(polim.cl.fρθ, sqrt(linspace(polim.cl.fθρ(θ1)^2, (polim.cl.fθρ(θ2))^2, 2N+1)))
    θmid = θdouble[2:2:2N]
    
    return θedges, θmid
end
weightedrings(polim::PolarImage, N::Integer) = weightedrings(polim, 0, pi/2, N)

function contactfreqs(polim::PolarImage, θ1::Real, θ2::Real, N::Integer, thresh)
    @checkθ1θ2
    θedges, θmid = weightedrings(polim, θ1, θ2, N)    
    K = zeros(N)
    for i = 1:N        
        logT = loggapfraction(pixels(polim, θedges[i], θedges[i+1]), thresh)
        K[i] = -logT * cos(θmid[i])
    end
    θedges, θmid, K
end
