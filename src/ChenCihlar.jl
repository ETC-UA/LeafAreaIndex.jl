
function gaplengths(polim::PolarImage, thresh, θ1::Real, θ2::Real)
    @checkθ1θ2
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

function probeprob(gaps::Vector{Int})
    sortedgaps = sort(gaps)
    lmax = sortedgaps[end]
    P = zeros(Int, lmax)
    for l = 1:lmax
        firstind = searchsortedfirst(sortedgaps, l)
        s = 0
        for λ in sortedgaps[firstind:end]
            s += λ - l
        end 
        P[l] = s
    end
    return(P/sum(P)) #normalize probability
end

function reducegaps(gaps, Lp, Wp; cutoff=1e-3)
    # theoretical accumulated gap fraction for random foliage (no clumping)
    F = Float64[(1 + Lp*λ/Wp) * exp(-Lp*(1+λ/Wp)) for λ = 0:maximum(gaps)]
    trunc_gaplength = searchsortedfirst(F, cutoff, rev=true)
    sortedgaps = sort(gaps)    
    gapsr = sortedgaps[1:searchsortedfirst(sortedgaps, trunc_gaplength)]
end

function Fmeasured(gaps, pixs)
    gaprange, gapcount = hist(gaps,0:maximum(gaps))
    blacks = length(pixs) - sum(gaps)
    Fm = 1 - cumsum([blacks, gapcount .* gaprange[2:end]])/length(pixs)
end

function Wp_estimation(P; n=3)
    regrcoef = linreg([1.:n], log(P[1:n]))
    Wp = - regrcoef[1] / abs(regrcoef[2])
end

function chencihlar(polim::PolarImage, thresh, θ1::Real, θ2::Real;
    cutoff=1e-3, Wp_est_n=3, itermax=10)

    @checkθ1θ2
    gaps = gaplengths(polim, thresh, θ1, θ2)
    pixs = pixels(polim, θ1, θ2)

    # initiation
    Fm = Fmeasured(gaps, pixs)
    Lp = -log(Fm[1])
    Wp =  Wp_estimation(probeprob(gaps); n=Wp_est_n)
    
    # first iteration with Wp
    gapsr = reducegaps(gapsr, Lp, Wp; cutoff=1e-3)
    Wp = Wp_estimation(probeprob(gapsr); n=Wp_est_n)
    Fmr = Fmeasured(gapsr, pixs)
    Lp = -log(Fmr[1])

    # loop Lp
    Lp_prev = Lp
    counter = 0
    while counter < itermax
        counter += 1
        gapsr = reducegaps(gapsr, Lp, Wp)
        Fmr = Fmeasured(gapsr, pixs)
        Lp = -log(Fmr[1])

        abs(1-Lp/Lp_prev) < 0.01 && break
        Lp_prev = Lp
    end
    counter == itermax && warn("chencihlar max iteration reached")
    Ω = log(Fm[1]) / log(Fmr[1]) * (1-Fmr[1])/(1-Fm[1])
end