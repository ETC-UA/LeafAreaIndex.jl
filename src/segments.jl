
# TODO only used in clumping.jl? Move it there or rewrite it there for view indices only.
"Return segments of pixelring [θ1, θ2] in n azimuth groups between 0 and 2π"
function segments(polim::PolarImage, θ1::Real, θ2::Real, n::Int)
    checkθ1θ2(θ1,θ2)
    ind_first, ind_last = firstlastind(polim, θ1, θ2)
    segmvec = [eltype(polim)[] for i = 1:n]
    imgsort = polim.imgsort
    ϕsort = getϕsort(polim)

    if isa(polim.mask, NoMask)
        adj = n/2π    
        for ind in ind_first:ind_last
            ϕ = ϕsort[ind]
            indn = ceil(Int, (ϕ + pi) * adj)
            push!(segmvec[indn], imgsort[ind])
        end
    else
        _, ϕmin, ϕmax = params(polim.mask)
        maskrange = ϕmax - ϕmin
        adj = n / (2π - maskrange)
        for ind in ind_first:ind_last
            ϕ = ϕsort[ind]
            ϕ += ifelse(ϕ < ϕmax, maskrange, 0)
            indn = ceil(Int, (ϕ + pi - maskrange) * adj)
            push!(segmvec[indn], imgsort[ind])
        end
    end
    segmvec
end

function firstlastind(polim::PolarImage{T, NoSlope, Mask}, θ1::Real, θ2::Real) where T
    checkθ1θ2(θ1,θ2)
    ind_first = searchsortedfirst(polim.mask.mask_ρ²sort, polim.cl.fθρ(θ1)^2)
    ind_last  = searchsortedlast( polim.mask.mask_ρ²sort, polim.cl.fθρ(θ2)^2)
    ind_first, ind_last
end

getϕsort(polim::PolarImage{T,<:SlopeInfo, NoMask}) where T = polim.cl.ϕsort
getϕsort(polim::PolarImage{T,<:SlopeInfo, Mask})   where T = polim.mask.mask_ϕsort

"Return the first and last index for the sorted image to abtain pixels in a given zenith range."
function firstlastind(polim::PolarImage{T, NoSlope, NoMask}, θ1::Real, θ2::Real) where T
    firstlastind(polim.cl, θ1, θ2)
end
function firstlastind(polim::PolarImage{T, Slope, NoMask}, θ1::Real, θ2::Real) where T
    α, ε = params(polim.slope)
    θmax = min(θ2, pi/2 - α)
    firstlastind(polim.cl, θ1, θmax)
end
function firstlastind(polim::PolarImage{T, NoSlope, Mask}, θ1::Real, θ2::Real) where T
    checkθ1θ2(θ1,θ2)
    ind_first = searchsortedfirst(polim.mask.mask_ρ²sort, polim.cl.fθρ(θ1)^2)
    ind_last  = searchsortedlast( polim.mask.mask_ρ²sort, polim.cl.fθρ(θ2)^2)
    ind_first, ind_last
end

"Return all image pixels between two given zenith angles."
function pixels(polim::PolarImage, θ1::Real, θ2::Real)
    ind_first, ind_last = firstlastind(polim, θ1, θ2)
    view(polim.imgsort, ind_first:ind_last)
end
pixels(polim::PolarImage) = pixels(polim, 0, π/2)

