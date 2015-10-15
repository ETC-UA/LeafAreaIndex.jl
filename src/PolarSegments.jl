# Returns segments of pixelring [θ1, θ2] in n azimuth groups between 0 and 2π
function segments(polim::PolarImage, θ1::Real, θ2::Real, n::Int)
    checkθ1θ2(θ1,θ2)

    ρ²indstart = searchsortedfirst(polim.cl.ρ²unique, polim.cl.fθρ(θ1)^2) 
      ρ²indend =  searchsortedlast(polim.cl.ρ²unique, polim.cl.fθρ(θ2)^2) 
    indstart = polim.cl.ρ²unique_ind[ρ²indstart]
      indend = polim.cl.ρ²unique_ind[ρ²indend]

    segmvec = [eltype(polim)[] for i = 1:n]
    
    imgsort = polim.imgsort
    ϕsort = polim.cl.ϕsort
    adj = n/2π    
    @inbounds for ind in indstart:indend
        indn = iceil((ϕsort[ind]+pi)*adj)
        push!(segmvec[indn], imgsort[ind])
    end
    segmvec
end
