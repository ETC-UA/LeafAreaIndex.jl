
# used in clumping (langxiang) and (currently unused) slope method of Espana
"Return segments of pixelring [θ1, θ2] in n azimuth groups between 0 and 2π"
function segments(polim::PolarImage, θ1::Real, θ2::Real, n::Int)
    checkθ1θ2(θ1,θ2)
    ind_first, ind_last = firstlastind(polim, θ1, θ2)
    segmvec = [eltype(polim)[] for i = 1:n]
    imgsort = polim.imgsort
    ϕsort = getϕsort(polim)

    if !hasmask(polim)
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

"Return the first and last index for the sorted image to abtain pixels in a given zenith range."
function firstlastind(polim::PolarImage, θ1::Real, θ2::Real)
    checkθ1θ2(θ1,θ2)
    
    ρ²sort = hasmask(polim) ? polim.mask.mask_ρ²sort : polim.cl.ρ²sort
    if hasslope(polim)
        α, ε = params(polim.slope)
        θ2 = min(θ2, pi/2 - α)
    end
    ind_first = searchsortedfirst(ρ²sort, polim.cl.fθρ(θ1)^2)
    ind_last  = searchsortedlast( ρ²sort, polim.cl.fθρ(θ2)^2)
    ind_first, ind_last
end

"Return all image pixels between two given zenith angles."
function pixels(polim::PolarImage, θ1::Real, θ2::Real)
    ind_first, ind_last = firstlastind(polim, θ1, θ2)
    view(polim.imgsort, ind_first:ind_last)
end
pixels(polim::PolarImage) = pixels(polim, 0, π/2)

