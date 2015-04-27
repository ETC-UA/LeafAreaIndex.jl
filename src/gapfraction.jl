# For slope correction on contact frequencies, divide each ring in N segments
const AZIMUTH_GROUPS = 8
# for ease of notation:
Ufixed = FixedPointNumbers.Ufixed

#fallback
gapfraction(pixs, thresh) = mean(pixs .> thresh)

#in general specialize on type of input array
function gapfraction(pixs::AbstractArray, thresh)
    threshT = convert(eltype(pixs), thresh)
    gapfrsum = zero(Int)        
    for pix in pixs
        gapfrsum += pix > threshT
    end
    return(gapfrsum / length(pixs))
end
# specialized gapfraction method for Ufixed pixels, with threshold already 
# converted to Ufixed type for type stability.
function gapfraction{T<:Ufixed}(pixs::AbstractArray{T}, thresh)
    threshT = convert(T, thresh)
    gapfrsum = zero(Int)        
    for pix in pixs
        gapfrsum += pix.i > threshT.i
    end
    return(gapfrsum / length(pixs))
end

function loggapfraction(pixs, thresh)
    gf = gapfraction(pixs, thresh)
    gf == zero(gf) && return(log(1/length(pixs)))
    log(gf)
end


function weightedrings(polim::PolarImage, θ1::Real, θ2::Real, N::Integer)
    
    # create edges for θ rings with similar number of pixels each
    θedges = map(polim.cl.fρθ, sqrt(linspace(polim.cl.fθρ(θ1)^2, 
                                             polim.cl.fθρ(θ2)^2, N+1)))
    #fix possible floating point roundoff errors
    θedges[1] = max(θedges[1], θ1) 
    θedges[end] = min(θedges[end], θ2)
    
    # weighted average midpoints
    θdouble = map(polim.cl.fρθ, sqrt(linspace(polim.cl.fθρ(θ1)^2, (polim.cl.fθρ(θ2))^2, 2N+1)))
    θmid = θdouble[2:2:2N]
    
    return θedges, θmid
end
weightedrings(polim::PolarImage, N::Integer) = weightedrings(polim, 0, pi/2, N)


function contactfreqs(polim::PolarImage, θ1::Real, θ2::Real, N::Integer, thresh)
    contactfreqs(polim, polim.slope, θ1, θ2, N, thresh)
end

function contactfreqs(polim::PolarImage, sl::NoSlope, θ1::Real, θ2::Real, 
                      N::Integer, thresh)
    @checkθ1θ2
    θedges, θmid = weightedrings(polim, θ1, θ2, N)    
    K = zeros(N)
    for i = 1:N        
        logT = loggapfraction(pixels(polim, θedges[i], θedges[i+1]), thresh)
        K[i] = -logT * cos(θmid[i])
    end
    θedges, θmid, K
end

function contactfreqs(polim::PolarImage, sl::Slope, θ1::Real, θ2::Real, N::Integer, 
                              thresh; Nϕ=AZIMUTH_GROUPS)
    @checkθ1θ2
    θedges, θmid = weightedrings(polim, θ1, θ2, N)    
    K = zeros(N)
    ϕv = midpoints(linspace(0, 2π, Nϕ+1))
    for i = 1:N
        adj = slope_adj(polim.slope, θmid[i], ϕv)
        # we divide each ring in Nϕ azimuth  segments, calculate the slope 
        # adjustment and loggapfraction per segment, then take average weighted 
        # by segment length.
        segm = segments(polim, θedges[i], θedges[i+1], Nϕ)
        lengths = Int[length(seg) for seg in segm]
        
        T = Float64[gapfraction(seg, thresh) for seg in segm]
        nz = find(T) # avoid 0.^(negative float)
        if isempty(nz)
            Tadj = 0.
        else
            Tadj = sum(T[nz].^(1./adj[nz]) .* lengths[nz])/sum(lengths)
        end        
        if Tadj == 0. #to avoid infinity with log, assume at least 1 sky pixel
            Tadj = 1 / sum(lengths)
        end        
        K[i] = -log(Tadj) * cos(θmid[i])

        # Doesnt work well...
        #logTs = Float64[loggapfraction(seg, thresh) for seg in segm]

        # Following gives same result as NoSlope for Slope(0,0)
        #logTs = loggapfraction(vcat(segm...), thresh)
        #K[i] = - cos(θmid[i]) * sum((logTs .* adj) .* lengths) / sum(lengths)
    end
    θedges, θmid, K
end
