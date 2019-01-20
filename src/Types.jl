const APPROX_N = 100 

#####
##### CameraLensParams
#####

"""
CameraLensParams stores the parameters of a combination of Lens + Camera.
Specifically, the image size, center of lens, maximum visible radius `ρmax` (in pixels) 
and the parameter of the lens projection function in a vector `A` such that
    `ρ = ρmax * sum(A[i] * θ^i for i = 1:length(A))`
    with `ρ` the distance to center of image in pixels and `θ` the viewing angle
"""
struct CameraLensParams
    imgsize    :: Tuple{Int,Int}
    lenscenter :: Tuple{Int,Int}
    ρmax       :: Int
    projcoefs  :: Vector{Float64}
end

Base.size(clp::CameraLensParams) = clp.imgsize

"Check consistency of CameraLensParams"
function check(clparams::CameraLensParams)
    Parameters.@unpack imgsize, lenscenter, ρmax, projcoefs = clparams
    
    height, width = imgsize
    ci, cj = lenscenter
    
    @assert height <= width "image detected in portrait mode"
    @assert width > 0 && height > 0 && ci > 0 && cj > 0

    abs(ci/height - 0.5) > 0.1 && @warn("cj ($ci) more than 10% away from center ($(height/2)).")
    abs(cj/width  - 0.5) > 0.1 && @warn("ci ($cj) more than 10% away from center ($(width/2).")

    @assert abs(ci/height - 0.5) < 0.2 "ci ($ci) more than 20% away from center ($(height/2)."
    @assert abs(cj/width - 0.5)  < 0.2 "cj ($cj) more than 20% away from center ($(width/2))."
    
    fθρ = make_proj(clparams)
    fρθ = make_proj_inv(fθρ)
    check_proj(fθρ, fρθ)
end

"""
Create the polynomial projection functions from given coefficients.
Returns function fθρ: θ [0, π/2] -> ρ [0, ρmax].
"""
function make_proj(clparams::CameraLensParams)
    Parameters.@unpack imgsize, lenscenter, ρmax, projcoefs = clparams

    A = projcoefs
    (length(A) == 1 || length(A) > 4) && return θ->ρmax * sum(A[i] * θ^i for i = 1:length(A))
    # @evalpoly is much more efficient if number of coefficients is known
    if length(A) == 2
        @inbounds fθρ(θ) = ρmax * @evalpoly(θ, 0, A[1], A[2])
        # fθρ(θ) = ρmax * (A[1]*θ + A[2]*θ^2)
    elseif length(A) == 3
        @inbounds fθρ(θ) = ρmax * @evalpoly(θ, 0, A[1], A[2], A[3])
    elseif length(A) == 4
        @inbounds fθρ(θ) = ρmax * @evalpoly(θ, 0, A[1], A[2], A[3], A[4])
    else
        fθρ(θ) = ρmax * sum(A[i] * θ^i for i = 1:length(A))
    end
    return fθρ
end

"""
Approximate the inverse projection functions given a projection function fθρ,
with fθρ: θ [0, π/2] -> ρ [0, ρmax].
Returns function fρθ through 3rd degree interpolation.
"""
function make_proj_inv(fθρ)
    θs = range(0, stop=π/2, length=APPROX_N)
    ρterms = [fθρ.(θ)^p for θ in θs, p in [1, 2, 3]]
    a, b, c = ρterms \ θs
    return ρ -> @evalpoly(ρ, 0, a, b, c)
end

function check_proj(fθρ, fρθ)
    @assert fθρ(0.0) >= 0 
    @assert fθρ(pi/2) >= 2 
    all(diff(map(fθρ, range(0, stop=pi/2, length=100))) .> 0) || error("Incorrectly defined fθρ; fθρ not monotonic")

    fρθ(0) < 0.0 && error("Incorrectly defined fρθ; fρθ(0) < 0.")
    fρθ(2) > pi/2 && error("Incorrectly defined fρθ; fρθ(2)> π/2")
    all(diff(map(fρθ, range(0, stop=pi/2, length=100))) .> 0) || error("Incorrectly defined fρθ; fρθ not monotonic")

    for θ in range(0, stop=pi/2, length=100)
         @assert isapprox(θ, fρθ(fθρ(θ)), atol=1e-3)
    end
end

function Base.show(io::IO, clp::CameraLensParams)
    println(io, "Parameters for CameraLens object with:")
    println(io, "\t image size:", clp.imgsize)
    println(io, "\t lens center:", clp.lenscenter)
    println(io, "\t max radius:", clp.ρmax)
    println(io, "\t projection function coefficients:", clp.projcoefs)
end

#####
##### CameraLensParams
#####

# Note UInt32 is sufficient to store ρ². FIXME premature optimization?
# mutable because we will only calculate `spiral_ind` when necessary (lazily)
# we need 
# * (inverse) projection function to go from view angle to polar radius in pixels
# * sort_ind to sort the image according to radius
# * ρ²sorted and ϕsorted to retrieve the indices to create circular segments
# * spiral indices for sorting the image according to a spiral for Chen Cihlar clumping
mutable struct CameraLens
    params     :: CameraLensParams
    fθρ        :: Function      # projection function θ → ρ (θ in [0,π/2], ρ in pixels)
    fρθ        :: Function      # inverse projection function ρ → θ
    sort_ind   :: Vector{Int}   # indices to sort according to ρ², then ϕ
    ρ²sort   :: Vector{UInt32}  # ρ² sorted by sort_ind
    ϕsort    :: Vector{Float64} # ϕ  sorted by sort_ind
    spiral_ind :: Union{Missing, Vector{Int}} # indices to sort according to spiral: first per pixel ρ², then per ϕ
end

"""
Construct the CameraLens object given CameraLensParams that contains the necessary
precomputed transformations to work with images taken from a particular Camera+Lens
combination.
"""
function CameraLens(clparams::CameraLensParams; storage_path::Union{Missing,String}=missing)
    check(clparams)
    ind, ρ²sort, ϕsort = get_sortedρ²ϕ(clparams; storage_path=storage_path)
    fθρ = make_proj(clparams)
    fρθ = make_proj_inv(fθρ)
    CameraLens(clparams, fθρ, fρθ, ind, ρ²sort, ϕsort, missing)
end

"Convenience constructor for CameraLens given parameters."
function CameraLens(imgsize::Tuple{Int,Int}, 
                    lenscenter::Tuple{Int,Int}, 
                    ρmax::Int,
                    projcoefs::Vector{Float64}; 
                    storage_path::Union{Missing,String}=missing)
    clparams = CameraLensParams(imgsize, lenscenter, ρmax, projcoefs)
    return CameraLens(clparams; storage_path=storage_path)
end


# TODO tests show that storage size is more than 100Mb, therefore not much
# quicker to use deepcopy (or serialize to workers)
"Retrieve the precomputed transformations from disk if `storage_path` given, else compute."
function get_sortedρ²ϕ(clparams::CameraLensParams; 
                       storage_path::Union{Missing,String}=missing)
    ismissing(storage_path) && (return make_sortedρ²ϕ(clparams))
    
    #isfile(storage_path) || throw(SystemError("storage_path not found: $storage_path."))
    isfile(storage_path) || close(JLD2.jldopen(storage_path, "w"))
    past_hashes = JLD2.jldopen(storage_path, "r") do file
        keys(file)
    end
    clhash = string(hash(clparams))
    if clhash in past_hashes
        return FileIO.load(storage_path, clhash)
    else
        res = make_sortedρ²ϕ(clparams)
        JLD2.jldopen(storage_path, "r+") do file
            file[clhash] = res
        end
        return res
    end
end

"Compute the image transformations for a given Camera+Lens."
function make_sortedρ²ϕ(clparams::CameraLensParams)
    Parameters.@unpack imgsize, lenscenter, ρmax, projcoefs = clparams

    ci, cj = lenscenter
    ρ² = zeros(UInt32, imgsize)
    ϕ = zeros(Float64, imgsize)
    
    for j = 1:imgsize[2]
        @fastmath @inbounds @simd for i = 1:imgsize[1]
            x = j - cj
            y = ci - i
            ρ²[i, j] = x^2 + y^2
            ϕ[i, j]  = atan(y, x)
        end
    end
    ϕ[ci, cj] = 0.0 #fix ϕ at center

    # sort by increasing ρ², then by increasing ϕ (ϕ ∈ [-π/2, π/2])
    ind = sortperm(reshape(ρ² + ϕ / 10, length(ρ²))) #this takes most time
    ρ²sort = ρ²[ind]
    ρ²indmax = searchsortedlast(ρ²sort, ρmax^2)
    ρ²sort = ρ²sort[1:ρ²indmax]
    ind = ind[1:ρ²indmax]
    ϕsort = ϕ[ind]
    return ind, ρ²sort, ϕsort
end

Base.size(cl::CameraLens) = size(cl.params)
hasspiral(cl::CameraLens) = !ismissing(cl.spiral_ind)

function Base.show(io::IO, cl::CameraLens)
    println(io, "CameraLens based on:")
    show(io, cl.params)
end

# "Generic constructor for testing"
function CameraLens(M::AbstractMatrix)
    height, width = size(M)
    ci, cj = ceil(Int, height/2), ceil(Int, width/2) #lenscenter
    ρmax = min(ci, cj, width-ci, height-cj)
    CameraLens(size(M), (ci, cj), ρmax, [1/(pi/2)])
end 

#####
##### Slope
#####

"""
Parameters (α, ϵ) that define a slope.
The slope inclination α is in radians [0, π/3].
The aspect ϵ is the downhill direction relative to the geographic north,
in clockwise direction on a map (i.e. in counter-clockwise direction on 
the picture). In radians [0, 2π]. Eg pi/2 means downhill to the East.
"""
struct SlopeParams
    α :: Float64
    ϵ :: Float64 
    function SlopeParams(α::Real, ϵ::Real)
        @assert 0 <= α < pi/3 
        @assert 0 <= ϵ <= 2π 
        new(float(α), float(ϵ))
    end
end

struct Slope
    params :: SlopeParams
    cosτsort  :: Vector{Float64} # incidence angle τ for sloped images, sorted as imgsort
end

function Slope(cl::CameraLens, sp::SlopeParams)
    α, ϵ = params(sp)
    θs = cl.fρθ.(sqrt.(cl.ρ²sort))
    cosτ = cos.(θs).*cos(α) + cos.(ϵ.-cl.ϕsort).*sin.(θs).*sin(α)
    Slope(sp, cosτ)
end
Slope(cl::CameraLens, α::Real, ϵ::Real) = Slope(cl, SlopeParams(α, ϵ))

params(sp::SlopeParams) = sp.α, sp.ϵ
params(sl::Slope) = params(sl.params)

# # auxiliary methods for slope correction (España et al 2007.)
# # used in `millergroup` and `millersimple`
# function slope_adj(sl::Slope, θ::Real, ϕ::Real)
#     # the cos(sl.α) term adjust for cartographic surface projection
#      1 / (cos(sl.α) * (1 + tan(θ)*cos(sl.ϵ-ϕ)*tan(sl.α)))
# end
# slope_adj(sl::Slope, θ::Real, ϕv::AbstractArray) = [slope_adj(sl, θ, ϕ) for ϕ in ϕv]

function Base.show(io::IO, sl::Slope) 
    α, ϵ = params(sl)
    println("Slope with ground inclination ", α, " for aspect direction ", ϵ)
end

#####
##### Mask
#####

"Parameters for masking the area larger than θ and between azimuth angles ϕmin and ϕmax (in range -π..π)"
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

# The same as those in CameraLens, but without the masked pixels
struct Mask
    params :: MaskParams
    mask_sort_ind   :: Vector{Int}   
    mask_ρ²sort     :: Vector{UInt32}  
    mask_ϕsort      :: Vector{Float64}
    mask_spiral_ind :: Union{Missing, Vector{Int}}
end

function Mask(cl::CameraLens, mp::MaskParams)
    θmask, ϕmin, ϕmax = params(mp)
    calc_spiral = hasspiral(cl)

    # Brute force approach, only calculated once per CameraLens + Mask 
    # combination. For coordinates (θ,ϕ) with θ>θmask, check every ϕ
    # and add to the two indices vectors and the two sorted vectors.
    ind_θmask = searchsortedfirst(cl.ρ²sort, cl.fθρ(θmask)^2)
    mask_sort_ind = cl.sort_ind[1:ind_θmask-1]
    mask_ρ²sort   = cl.ρ²sort[1:ind_θmask-1]
    mask_ϕsort    = cl.ϕsort[1:ind_θmask-1]
    calc_spiral && (mask_spiral_ind = cl.spiral_ind[1:ind_θmask-1])
    for i = ind_θmask:length(cl.sort_ind)
        ϕ = cl.ϕsort[i]
        if !(ϕmin < ϕ < ϕmax)
            push!(mask_sort_ind, cl.sort_ind[i])
            push!(mask_ρ²sort, cl.ρ²sort[i])
            push!(mask_ϕsort, cl.ϕsort[i])
            calc_spiral && push!(mask_spiral_ind, cl.spiral_ind[i])
        end
    end
    calc_spiral || (mask_spiral_ind = missing)
    Mask(mp, mask_sort_ind, mask_ρ²sort, mask_ϕsort, mask_spiral_ind)
end

"Convenience constructor for Mask"
function Mask(cl::CameraLens, θ::Real, ϕmin::Real, ϕmax::Real)
    Mask(cl, MaskParams(θ, ϕmin, ϕmax))
end

params(m::MaskParams) = m.θ, m.ϕmin, m.ϕmax
params(m::Mask) = params(m.params)

function Base.show(io::IO, mask::Mask) 
    θmask, ϕmin, ϕmax = params(mask)
    println("Mask for θ >", θmask, " and ϕ in [", ϕmin, ",", ϕmax,"].")
end

#####
##### PolarImage
#####

# The original image is required for EdgeDetection
# imgsort: image sorted by ρ², then ϕ
# imgspiral: image sorted by spiral_ind (lazily)
"Type to contain an image and its polar transformations"
mutable struct PolarImage{T}
    cl        :: CameraLens
    img       :: AbstractArray{T, 2}
    imgsort   :: Vector{T}       
    imgspiral :: Union{Missing, Vector{T}}
    slope     :: Union{Missing, Slope}
    mask      :: Union{Missing, Mask}
end

# PolarImage constructors. Note that masking and slope are mutually exclusive (for now no crops on slopes considered).
function PolarImage(img::AbstractMatrix, cl::CameraLens, 
        slope::Union{Missing, Slope, SlopeParams}=missing, 
        mask ::Union{Missing, Mask,  MaskParams }=missing)
    hasslope = !ismissing(slope)
    hasmask  = !ismissing(mask)
    if hasslope & hasmask
        error("Masking (for crops) together slope is not yet supported. ")
    end
    @assert size(img) == size(cl)

    imgsort   = img[cl.sort_ind]
    imgspiral = hasspiral(cl) ? img[cl.spiral_ind] : missing
    isa(slope, SlopeParams) && (slope = Slope(cl, slope))
    isa(mask,  MaskParams ) && (mask  = Mask(cl, mask))
    PolarImage{eltype(img)}(cl, img, imgsort, imgspiral, slope, mask)
end

hasslope(polim::PolarImage) = !ismissing(polim.slope)
hasmask(polim::PolarImage)  = !ismissing(polim.mask)
getϕsort(polim::PolarImage) = hasmask(polim) ? polim.mask.mask_ϕsort : polim.cl.ϕsort

Base.eltype(polim::PolarImage{T}) where T = T
Base.length(pm::PolarImage) = length(getϕsort(pm))
Base.size(pm::PolarImage) = size(pm.cl)

function Base.show(io::IO, polim::PolarImage)
    slope = hasslope(polim) ? "with" : "without"
    mask  = hasmask(polim) ? "with" : "without"
    print(io::IO, "PolarImage $slope slope and $mask mask")
end

# TODO add get methods for (lazy) spiral accessors