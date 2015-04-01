
function segments(polim::PolarImage, θ1::Real, θ2::Real, n::Int)
    @checkθ1θ2

    θ1ind = searchsortedfirst(polim.cl.ρ²unique, polim.cl.fθρ(θ1)^2) 
    θ2ind = searchsortedlast(polim.cl.ρ²unique, polim.cl.fθρ(θ2)^2) 
    
    segmvec = [eltype(polim)[] for i = 1:n]
    
    imgsort = polim.imgsort
    ϕsort = polim.cl.ϕsort
    adj = n/2pi
    @inbounds for ind in polim.cl.ρ²Ncs[θ1ind]:polim.cl.ρ²Ncs[θ2ind]
        indn = iceil((ϕsort[ind]+pi)*adj)
        push!(segmvec[indn], imgsort[ind])
    end
    segmvec
end

function spiral(cl::CameraLens)
    ρmax = cl.fθρ(π/2)
    #ρ²Ncs = cl.ρ²Ncs
    ρ²sort = cl.ρ²sort
    ϕsort = cl.ϕsort
    ind_prev = 1
    spiralind = Int[]
    ρ²spiralind = Int[]
    ind = 0
    for ρ² in [1:ρmax].^2
        #TODO improve by searching ρ²unique and taking index from ρ²Ncs
        ind = searchsortedfirst(ρ²sort, ρ²)
        ρ²spiralind = sortperm(ArrayViews.view(ϕsort, ind_prev:ind-1))
        append!(spiralind, ind_prev - 1 + ρ²spiralind)
        ind_prev = ind
    end
    # last circle of spiral
    ρ²spiralind = sortperm(ϕsort[ind_prev:end])
    append!(spiralind, ind_prev - 1 + ρ²spiralind)
    
    length(unique(spiralind)) == length(ϕsort) || error("length check")
    return(spiralind)
end