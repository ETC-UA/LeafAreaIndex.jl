## PolarPixels ##

#pixels(polim::PolarImage, θ::Real) = pixels(polim, zero(θ), θ)
pixels(polim::PolarImage) = pixels(polim, 0.0, π/2)

# for just polar distance, returning a view is sufficient and much faster
function pixels(polim::PolarImage, θ1::Real, θ2::Real) 
    checkθ1θ2(θ1,θ2)
    ind_first = searchsortedfirst(polim.cl.ρ²sort, polim.cl.fθρ(θ1)^2)
    ind_last = searchsortedlast(polim.cl.ρ²sort, polim.cl.fθρ(θ2)^2) 
    ArrayViews.view(polim.imgsort, ind_first:ind_last)
end

# # Similar to PolarRings, but this will output only pixels. Use pixels(polim,θ1,θ2)
# type PolarPixels{T,A}
#     ind_first::Int
#     ind_last::Int
#     polim::PolarImage{T,A}
# end

# # PolarPixels constructor
# function PolarPixels(polim::PolarImage, θ1::Real, θ2::Real)
#     checkθ1θ2(θ1,θ2)
#     ind_first = searchsortedfirst(polim.cl.ρ²sort, polim.cl.fθρ(θ1)^2)
#     ind_last = searchsortedlast(polim.cl.ρ²sort, polim.cl.fθρ(θ2)^2) 
#     PolarPixels{eltype(polim.img), typeof(polim.img)}(ind_first, ind_last, polim)
# end

# # PolarPixels is an iterator as soon as start, done and next are defined
# Base.start(pc::PolarPixels) = pc.ind_first
# Base.done(pc::PolarPixels, i) = i == pc.ind_last+1
# Base.next{T,A}(pc::PolarPixels{T,A}, i) = (pc.polim.imgsort[i]::T, i+1)

# Base.length(pc::PolarPixels) = pc.ind_last - pc.ind_first + 1
# Base.endof(pc::PolarPixels) = pc.polim.imgsort[pc.ind_last]

# #pixels{T,A}(polim::PolarImage{T,A}, θ1::Real, θ2::Real) = collect(T, PolarPixels(polim,θ1,θ2))

