## PolarRings ##

# We create an extra type PolarRings to iterate over pixels with increasing ρ.
# Each iteration outputs a tuple (ρ²,ϕ,pixels::Vector{eltype(polim.data)})
type PolarRings{T,A}
    polim::PolarImage{T,A}
    ind_first::Int
    ind_last::Int
    ind::Vector{Int}
end 

function PolarRings(polim::PolarImage, θ1::Real, θ2::Real)
    θ1 < 0 && throw(DomainError())        
    θ2 > π/2 && throw(DomainError())
    θ2 < θ1 && error("θ2 < θ1")
    ind_first = searchsortedfirst(polim.cl.ρ²unique, polim.cl.fθρ(θ1)^2) 
    ind_last = searchsortedlast(polim.cl.ρ²unique, polim.cl.fθρ(θ2)^2) 
    ind = [0, polim.cl.ρ²Ncs] #add 0 for start of first range
    PolarRings{eltype(polim.img),typeof(polim.img)}(polim, ind_first, ind_last, ind)
end

# Defining start, done and next methods will make PolarRings an iterator
Base.start(pr::PolarRings) = pr.ind_first

Base.done(pr::PolarRings, i) = i == pr.ind_last+1

function Base.next(pr::PolarRings, i) 
    vecstart = pr.ind[i]+1
    vecend = pr.ind[i+1]
    ρ² = pr.polim.cl.ρ²unique[i]
    ϕvec = sub(pr.polim.cl.ϕsort,vecstart:vecend)
    pixs = sub(pr.polim.imgsort, vecstart:vecend)
    ((ρ²,ϕvec, pixs), i+1)
end

rings(polim::PolarImage, θ1::Real, θ2::Real) = PolarRings(polim, θ1, θ2)
rings(polim::PolarImage, θ::Real) = rings(polim, zero(θ), θ)
rings(polim::PolarImage) = rings(polim, pi/2)
