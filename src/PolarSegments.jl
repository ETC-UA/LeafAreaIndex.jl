
function segments(polim::PolarImage, θ1::Real, θ2::Real, n::Int)
    @checkθ1θ2

    θ1ind = searchsortedfirst(polim.cl.ρ²unique, polim.cl.fθρ(θ1)^2) 
    θ2ind = searchsortedlast(polim.cl.ρ²unique, polim.cl.fθρ(θ2)^2) 
    
    segmvec = [eltype(polim)[] for i = 1:n]
    
    imgsort = polim.imgsort
    ϕsort = polim.cl.ϕsort
    adj = n/2π
    @inbounds for ind in polim.cl.ρ²Ncs[θ1ind]:polim.cl.ρ²Ncs[θ2ind]
        indn = iceil((ϕsort[ind]+pi)*adj)
        push!(segmvec[indn], imgsort[ind])
    end
    segmvec
end
