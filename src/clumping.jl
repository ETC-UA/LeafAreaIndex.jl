

function langxiang(thresh, polim::PolarImage, θ1::Real, θ2::Real, n::Integer)
    @checkθ1θ2
    segm = segments(polim, θ1, θ2, n)
    
    segm_gapfr = Float64[gapfraction(segm[i], thresh) for i in 1:length(segm)]
    clump_LX = log(mean(segm_gapfr)) / mean(log(segm_gapfr))
end

function langxiang45(thresh, polim::PolarImage, θ1::Real, θ2::Real)
	langxiang(thresh, polim, θ1, θ2, 8)
end

# function langxiang45(polim::PolarImage, θ1::Real, θ2::Real)
#     @checkθ1θ2
#     thresh = threshold(polim)
#     langxiang45(thresh, polim, θ1, θ2)
# end


function gaplengths(polim::PolarImage, thresh, θ1::Real, θ2::Real)
	@checkθ1θ2
	θ1ind = searchsortedfirst(polim.cl.ρ²unique, polim.cl.fθρ(θ1)^2) 
    θ2ind = searchsortedlast(polim.cl.ρ²unique, polim.cl.fθρ(θ2)^2) 
	
    threshT = convert(eltype(PolarImage), thresh)
    λ = 0
    λvec = Int[]
    for pixel in polim.imgspiral[θ1ind:θ2ind]
        if pixel < threshT
            λ != 0 && push!(λvec, λ)
            λ = 0
        else
            λ += 1
        end
    end
    return(λvec)
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