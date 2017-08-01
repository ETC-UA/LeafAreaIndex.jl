abstract type MaskInfo end
struct NoMask <: MaskInfo; end

struct MaskParams
    θ :: Float64
    ϕmin :: Float64
    ϕmax :: Float64
end
function MaskParams(θ::Real, ϕmin::Real, ϕmax::Real)
    # masking a view angel smaller than pi/4 seems excessive
    @assert pi/4 <= θ < pi/2 
    @assert -π < ϕmin < ϕmax < π
    Mask(θ, ϕmin, ϕmax)
end

struct Mask <: MaskInfo
    params :: MaskParams
    # The same as those in CameraLens, but without the masked pixels
    mask_sort_ind :: Vector{Int}   # indices to sort according to ρ², then ϕ
    mask_ρ²sort :: Vector{UInt32}  # ρ² sorted by sort_ind
    mask_ϕsort  :: Vector{Float64} # ϕ  sorted by sort_ind
    mask_spiral_ind :: Vector{Int} # indices to sort according to spiral: first per pixel ρ², then per ϕ
end
params(m::MaskParams) = m.θ, m.ϕmin, m.ϕmax
params(m::Mask) = params(m.params)

function Mask(cl::CameraLens, mp::MaskParams)
    θmask, ϕmin, ϕmax = params(mp)

    # Brute force approach, only calculated once per CameraLens + Mask 
    # combination. For coordinates (θ,ϕ) with θ>θmask, check every ϕ
    # and add to the two indices vectors and the two sorted vectors.
    ind_θmask = searchsortedfirst(cl.ρ²sort, cl.fθρ(θmask)^2)
    mask_sort_ind   = cl.sort_ind[1:ind_θmask-1]
    mask_ρ²sort     = cl.ρ²sort[1:ind_θmask-1]
    mask_ϕsort      = cl.ϕsort[1:ind_θmask-1]
    mask_spiral_ind = cl.spiral_ind[1:ind_θmask-1]
    for i = ind_θmask:length(cl.sort_ind)
        ϕ = cl.ϕsort[i]
        if !(ϕmin < ϕ < ϕmax)
            push!(mask_sort_ind, cl.sort_ind[i])
            push!(mask_ρ²sort, cl.ρ²sort[i])
            push!(mask_ϕsort, cl.ϕsort[i])
            push!(mask_spiral_ind, cl.spiral_ind[i])
        end
    end
    Mask(mp, mask_sort_ind, mask_ρ²sort, mask_ϕsort, mask_spiral_ind)
end

function Mask(cl::CameraLens, θ::Real, ϕmin::Real, ϕmax::Real)
    params = MaskParams(θ, ϕmin, ϕmax)
    Mask(cl, params)
end

function Base.show(io::IO, mask::Mask) 
    θmask, ϕmin, ϕmax = params(mask)
    println("Mask for θ >", θmask, " and ϕ in [", ϕmin, ",", ϕmax,"].")
end