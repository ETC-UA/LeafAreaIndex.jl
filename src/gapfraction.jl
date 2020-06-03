# For slope correction on contact frequencies:
const INCIDENCE_GROUPS = 360 #number of τ groups over [0, pi/2]
const MAX_ITER_τ = 5 # (Schleppi, 2007) says "after a few cycles"
const SLOPE_TOL = 1e-3

#fallback
# WHY??
#gapfraction(pixs, thresh) = mean(pixs .> thresh)

#in general specialize on type of input array
function gapfraction(pixs::AbstractArray, thresh)
    threshT = convert(eltype(pixs), thresh)
    gapfrsum = 0
    for pix in pixs
        gapfrsum += pix > threshT
    end
    return(gapfrsum / length(pixs))
end

function loggapfraction(pixs, thresh)
    gf = gapfraction(pixs, thresh)
    gf == zero(gf) && return(log(1/length(pixs)))
    log(gf)
end

"Create rings with similar amount of pixels per ring. Outputs the edges and weigted midpoints of the rings."
function weightedrings(polim::PolarImage, θ1::Real, θ2::Real, N::Integer)
    
    fθρ = polim.cl.fθρ
    # create edges for θ rings with similar number of pixels each
    θedges = map(polim.cl.fρθ, sqrt.(range(fθρ(θ1)^2, stop=fθρ(θ2)^2, length=N+1)))
                                             
    #fix possible floating point roundoff errors
    θedges[1] = max(θedges[1], θ1) 
    θedges[end] = min(θedges[end], θ2)
    
    # weighted average midpoints
    θdouble = map(polim.cl.fρθ, sqrt.(range(fθρ(θ1)^2, stop=(fθρ(θ2))^2, length=2N+1)))
    θmid = θdouble[2:2:2N]
    
    return θedges, θmid
end

weightedrings(polim::PolarImage, N::Integer) = weightedrings(polim, 0, pi/2, N)

function contactfreqs(polim::PolarImage, θ1::Real, θ2::Real, N::Integer, thresh; 
                      Nτ=INCIDENCE_GROUPS, max_iter=MAX_ITER_τ, tol=SLOPE_TOL)
    checkθ1θ2(θ1,θ2)
    θedges, θmid = weightedrings(polim, θ1, θ2, N)    
    K = zeros(N)
    
    ## WITHOUT SLOPE ##
    if !hasslope(polim)
        for i = 1:N        
            logT = loggapfraction(pixels(polim, θedges[i], θedges[i+1]), thresh)
            K[i] = -logT * cos(θmid[i])
        end
        return θedges, θmid, K
    end

    ## WITH  SLOPE ##
    Nτ = max(Nτ, N) #at least as many incidence groups as zenith groups required
    for i = 1:N
        pixs = pixels(polim, θedges[i], θedges[i+1])
  
        ind_first, ind_last = firstlastind(polim, θedges[i], θedges[i+1])
        cosτ = view(polim.slope.cosτsort, ind_first:ind_last)
        
        K[i] = contactfreqs_iterate(pixs, cosτ, thresh, θmid[i]; 
                                    Nτ=Nτ, max_iter=max_iter, tol=tol)
  
        # Method España et al 2007. Nϕ different here!
        # adj = slope_adj(polim.slope, θmid[i], ϕv)
        # # we divide each ring in Nϕ azimuth  segments, calculate the slope 
        # # adjustment and loggapfraction per segment, then take average weighted 
        # # by segment length.
        # segm = segments(polim, θedges[i], θedges[i+1], Nϕ)
        # lengths = Int[length(seg) for seg in segm]
        
        # T = Float64[gapfraction(seg, thresh) for seg in segm]
        # nz = find(T) # avoid 0.^(negative float)
        # if isempty(nz) #to avoid infinity with log, assume at least 1 sky pixel
        #     Tadj = 1 / sum(lengths)
        # else
        #     Tadj = sum(T[nz].^(1./adj[nz]) .* lengths[nz])/sum(lengths)
        # end      
        # K[i] = -log(Tadj) * cos(θmid[i])
  
    end
    θedges, θmid, K
end
# Method Schleppi et al 2007
function contactfreqs_iterate(pixs::AbstractArray, cosτ::AbstractArray, thresh, θ::Float64;
        Nτ=INCIDENCE_GROUPS, max_iter=MAX_ITER_τ, tol=SLOPE_TOL)

    τmax = π/2
    τ = StatsBase.midpoints(range(0, stop=τmax, length=Nτ+1))    
    Aθτ = fasthist(acos.(cosτ), -1/Nτ : τmax/Nτ : τmax)

    iter = 0
    # initially start with contact frequency K from whole θ ring 
    logT = loggapfraction(pixs, thresh)
    K = - logT * cos(θ)
    while iter < max_iter
        iter += 1 
        sum(Aθτ) == 0 && break
        Tnew = sum(Aθτ .* exp.(-K ./ cos.(τ))) / sum(Aθτ)
        logTnew = log(Tnew)
        abs(logTnew / logT - 1) < tol && break
        K *= logTnew / logT
        logT = logTnew
    end
    K
end


