const LANGXIANG = 8

# specialize on slope type
function langxiang(polim::PolarImage, thresh, θ1::Real, θ2::Real, nϕ::Integer)
    langxiang(polim, polim.slope, thresh, θ1, θ2, nϕ)
end

# Without Slope
function langxiang(polim::PolarImage, sl::NoSlope, thresh, θ1::Real, θ2::Real, 
                   nϕ::Integer)
    checkθ1θ2(θ1,θ2)
    segm = segments(polim, θ1, θ2, nϕ)
    
    segm_gapfr = Float64[gapfraction(segm[i], thresh) for i in 1:length(segm)]
    clump_LX = log.(mean(segm_gapfr)) / mean(log.(segm_gapfr))
end

# With Slope
function langxiang(polim::PolarImage, sl::Slope, thresh, θ1::Real, θ2::Real, 
                   nϕ::Integer)
    checkθ1θ2(θ1,θ2)
    
    ρ²indstart = searchsortedfirst(polim.cl.ρ²unique, polim.cl.fθρ(θ1)^2) 
      ρ²indend =  searchsortedlast(polim.cl.ρ²unique, polim.cl.fθρ(θ2)^2) 
    indstart = polim.cl.ρ²unique_ind[ρ²indstart]
      indend = polim.cl.ρ²unique_ind[ρ²indend]

    segmvec = [eltype(polim)[] for i = 1:nϕ]
    τvec = [Float64[] for i = 1:nϕ]
    
    imgsort = polim.imgsort
    ϕsort = polim.cl.ϕsort
    τsort = polim.τsort
    adj = nϕ/2π #adjustment to segment ϕsort
    @inbounds for ind in indstart:indend
        indn = ceil(Int, (ϕsort[ind] + pi) * adj)
        push!(segmvec[indn], imgsort[ind])
        push!(τvec[indn], τsort[ind])
    end
    
    # TODO use weighted average θ instead of simple average
    θ = (θ1 + θ2) / 2
    K = [contactfreqs_iterate(segmvec[i], τvec[i],thresh, θ) for i = 1:nϕ]
    T = exp( -K / cos(θ))
    for i = 1:nϕ
        if T[i] == 0.0            
            T[i] = 1 / length(segmvec[i])
        end
    end
    clump_LX = log.(mean(T)) / mean(log.(T))
end



function langxiang45(polim::PolarImage, thresh, θ1::Real, θ2::Real)
	langxiang(polim, thresh, θ1, θ2, LANGXIANG)
end

# for Chen Cihlar see 'ChenCihlar.jl'