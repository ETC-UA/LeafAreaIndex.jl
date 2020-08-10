const WP_EST_N = 3
const WP_EST_θ1 = 10 / 180 * pi
const WP_EST_θ2 = 80 / 180 * pi
const F_CUTOFF = 1e-3
const CC_ITER_MAX = 10

function gaplengths(polim::PolarImage, thresh, θ1::Real, θ2::Real)
    checkθ1θ2(θ1,θ2)
    θ1ind = searchsortedfirst(polim.cl.ρ²sort, polim.cl.fθρ(θ1)^2) 
    θ2ind = searchsortedlast(polim.cl.ρ²sort, polim.cl.fθρ(θ2)^2) 

    threshT = convert(eltype(PolarImage), thresh)
    λ = 0 #gaplength
    gaps = Int[]
    for pixel in polim.imgspiral[θ1ind:θ2ind]
        if pixel < threshT
            λ != 0 && push!(gaps, λ)
            λ = 0
        else
            λ += 1
        end
    end
    return(gaps)
end

# methods for P approach
# Probe probability from formulae 17-18, Chen & Cihlar 1995, IEEE
function probeprob(gaps::Vector{Int})
    length(gaps) == 0 && return([1])
    sortedgaps = sort(gaps)
    len_max = sortedgaps[end]
    P = zeros(Int, len_max)
    for len = 1:len_max
        firstind = searchsortedfirst(sortedgaps, len)
        s = 0
        for λ in sortedgaps[firstind:end]
            s += λ - len
        end 
        P[len] = s
    end
    return(P/sum(P)) #normalize probability
end


#Methods for F approach
function reducegaps(gaps::Vector{Int}, Lp, Wp; cutoff=F_CUTOFF)
    # F = theoretical accumulated gap fraction for random foliage (no clumping)
    F = Float64[(1 + Lp*λ/Wp) * exp(-Lp*(1+λ/Wp)) for λ = 0:maximum(gaps)]
    trunc_gaplength = searchsortedfirst(F, cutoff, rev=true)
    sortedgaps = sort(gaps)    
    gapsr = sortedgaps[1:searchsortedlast(sortedgaps, trunc_gaplength)]
end
function Fmeasured(gaps::Vector{Int}, pixs)
    gaprange, gapcount = hist(gaps, 0:maximum(gaps))
    blacks = length(pixs) - sum(gaps)
    Fm = 1 - cumsum([blacks, gapcount .* gaprange[2:end]])/length(pixs)
end



# We estimate the projected average leaf width once per image, because
# it's not stable for small zenith angles θ and "the calculation of 
# [clumping index] Ω with a gap-size distribution is not very sensitive to the
# choice of the elementh width within a reasonable range" (Chen,Chilar 1995, 
# applied optics, section 4.B).
# Wp estimation according to Leblanc et al, 2005, formula (7).
function Wp_estimation(polim::PolarImage, thresh; 
    θ1=WP_EST_θ1, θ2=WP_EST_θ2, Wp_est_n=WP_EST_N)
    
    checkθ1θ2(θ1,θ2)
    gaps = gaplengths(polim, thresh, θ1, θ2)
    Wp =  Wp_estimation(probeprob(gaps); n=Wp_est_n)
end

function Wp_estimation(P; n=WP_EST_N)
    regrcoef = linreg([1.:n], log(P[1:n]))
    Wp = - regrcoef[1] / abs(regrcoef[2])
end

function chencihlar(polim::PolarImage, thresh, θ1::Real, θ2::Real; kwargs...)
    checkθ1θ2(θ1,θ2)
    # Wp estimation not stable for small zenith angles, so estimate Wp from full image
    if θ1 > WP_EST_θ1
        #TODO add Wp_est_n keyword argument
        Wp = Wp_estimation(probeprob(gaplengths(polim, thresh, θ1, θ2)))
    else
        Wp = Wp_estimation(polim, thresh)
    end
    chencihlarF(Wp, polim, thresh, θ1, θ2; kwargs...)    
end

# P method Chen Cihlar, requires Wp input
function chencihlarF(Wp::Float64, polim::PolarImage, thresh, θ1::Real, θ2::Real;
    cutoff=F_CUTOFF, Wp_est_n=WP_EST_N, itermax=CC_ITER_MAX)

    checkθ1θ2(θ1,θ2)
    gaps = gaplengths(polim, thresh, θ1, θ2)
    length(gaps) == 0 && return(Ω = 1)
    pixs = pixels(polim, θ1, θ2)


    # Lp initiation
    Fm = Fmeasured(gaps, pixs)
    Lp = -log(Fm[1])    

    # loop Lp
    Lp_prev = Lp
    counter = 0
    Fmr = Float64[] # initiation required before loop        
    while counter < itermax
        counter += 1
        gaps = reducegaps(gaps, Lp, Wp)        
        Fmr = Fmeasured(gaps, pixs)
        Lp = -log(Fmr[1])
        abs(1-Lp/Lp_prev) < 0.01 && break
        Lp_prev = Lp
    end
    counter == itermax && @warn("chencihlar max iteration reached")
    Ω = log(Fm[1]) / log(Fmr[1]) * (1-Fmr[1])/(1-Fm[1])    
end

# P-method from Chen & Cihlar, IEEE trans. geo rem. sens. 1995
# Gives low clumping factors.
function chencihlarP(polim::PolarImage, thresh, θ1::Real, θ2::Real;n=WP_EST_N)
    checkθ1θ2(θ1,θ2)
    gaps = gaplengths(polim, thresh, θ1, θ2)
    Pl = probeprob(gaps)
    length(Pl) < 2 && return(1)
    maxP = length(Pl)
    n = min(n, maxP)
    Leθ = -linreg([1.:n], log(Pl[1:n]))[1]

    Pstart = ifloor(maxP/4)
    Pstop = ifloor(maxP/3)
    
    Lcθ = -linreg(float([Pstart:Pstop]), log(Pl[Pstart:Pstop]))[1]
    α = 0
    LEp = log((1+α)*Lcθ*exp(-Lcθ)/(sqrt(2(1+α)exp(-(Leθ+Lcθ))-(1+2α)exp(-2Lcθ))-exp(-Lcθ)))
    α = Lcθ*exp(-LEp)/3 #higher order correction
    LEp = log((1+α)*Lcθ*exp(-Lcθ)/(sqrt(2(1+α)exp(-(Leθ+Lcθ))-(1+2α)exp(-2Lcθ))-exp(-Lcθ)))
    Ω = Leθ / (LEp*Lcθ)
end