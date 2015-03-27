
const MILLER_GROUPS = 10
const LANG_START = 25/180*pi
const LANG_END = 65/180*pi
const LANG_STEPS = 50

# Constant view angle

function zenith57(polim::PolarImage, thresh::Real; ringwidth=RING_WIDTH)
    G = 0.5 
    ring = pixels(polim, 1 - ringwidth, 1 + ringwidth)    
    gf57 = gapfraction(ring, convert(eltype(ring), thresh))
    return( -log(gf57) / G * cos(1))
end
zenith57(polim::PolarImage) = zenith57(polim, threshold(polim))

# Miller's formule with constant leaf angel

miller(polim::PolarImage, thresh::Real) = millergroup(polim, thresh)
miller(polim::PolarImage) = millergroup(polim, threshold(polim))

function millersimple(polim::PolarImage, thresh::Real)
    prevθ = 0.
    s = 0.
    #for faster gapfraction, convert thresh type before loop
    threshT = convert(eltype(polim), thresh) 
    #define inverse function for ρ²
    fρ²θ(ρ²) = polim.cl.fρθ(sqrt(ρ²)) 
    for (ρ², ϕ, px) in rings(polim)
        θ = fρ²θ(ρ²) 
        dθ = θ - prevθ
        P = gapfraction(px, threshT)
        s -= ifelse(P==0., 0., log(P) * cos(θ) * sin(θ) * dθ)
        prevθ = θ
    end
    return(2s)
end
millersimple(polim::PolarImage) = millersimple(polim, threshold(polim))

function millergroup(polim::PolarImage, group::Integer, thresh::Real)
    s = 0.    
    prevθ = 0.
    count = 0
    pixs = eltype(polim)[]
    avgθ = StreamStats.Mean()
    fρ²θ(ρ²) = polim.cl.fρθ(sqrt(ρ²)) 
    #for faster gapfraction, convert thresh type before loop
    threshT = convert(eltype(polim), thresh) 
    for (ρ², ϕ, px) in rings(polim)       
        count += 1        
        StreamStats.update!(avgθ, fρ²θ(ρ²))
        append!(pixs, px)
        
        if count == group
            θ = StreamStats.state(avgθ)
            dθ = θ - prevθ
            P = gapfraction(pixs, threshT)            
            s -= ifelse(P==0., 0., log(P) * cos(θ) * sin(θ) * dθ)        
            
            prevθ = θ
            count = 0
            empty!(pixs)
            empty!(avgθ)
        end            
    end    
    return(2s)
end
millergroup(polim::PolarImage, thresh::Real) = millergroup(polim, MILLER_GROUPS, thresh)
millergroup(polim::PolarImage) = millergroup(polim, MILLER_GROUPS, threshold(polim))

function millerrings(polim::PolarImage, N::Integer, thresh)
    rings = linspace(0, π/2, N+1)
    dθ = π/2/N
    s = 0.
    #for faster gapfraction, convert thresh type before loop
    threshT = convert(eltype(polim), thresh) 
    for i = 1:N
        px = pixels(polim, rings[i], rings[i+1])
        P = gapfraction(px, threshT)
        θ = i * dθ - dθ/2
        s -= ifelse(P==0., 0., log(P) * cos(θ) * sin(θ) * dθ)  
    end
    return(2s)
end

function millerrings(polim::PolarImage, thresh::Real) 
    Nrings = iceil(polim.cl.fθρ(pi/2) / MILLER_GROUPS)
    millerrings(polim, Nrings, thresh)
end
millerrings(polim::PolarImage) = millerrings(polim, threshold(polim))

function lang(polim::PolarImage, thresh::Real) 
    θvals = linspace(LANG_START, LANG_END, LANG_STEPS)
    Kapprox = Float64[contactfreq(polim, θ, thresh) for θ = θvals]
    2*sum(linreg(θvals, Kapprox))
end
