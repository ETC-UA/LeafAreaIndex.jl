## PolarImage ##
const APPROX_N=100

"Type to contain an image and its polar transform"
struct PolarImage{T, S <:SlopeInfo, M <: MaskInfo}
    cl :: CameraLens
    slope :: S
    mask :: M
    img                      # original image, required for EdgeDetection
    imgsort :: Vector{T}     # image sorted by ρ², then ϕ
    imgspiral :: Vector{T}   # image spiral sorted
end

# PolarImage constructors. Note that masking and slope are mutually exclusive (for now no crops on slopes considered).
function PolarImage(img::AbstractMatrix, cl::CameraLens, sl::SlopeInfo)
    @assert size(img) == size(cl)
    imgsort   = img[cl.sort_ind]
    imgspiral = img[cl.spiral_ind]
    PolarImage{eltype(img), typeof(sl), NoMask}(cl, sl, NoMask(), img, imgsort, imgspiral)
end

PolarImage(img::AbstractMatrix, cl::CameraLens) = PolarImage(img, cl, NoSlope())

function PolarImage(img::AbstractMatrix, cl::CameraLens, mask::Mask)
    @assert size(img) == size(cl)
    imgsort   = img[mask.mask_sort_ind]
    imgspiral = img[mask.mask_spiral_ind]
    PolarImage{eltype(img), NoSlope, Mask}(cl, NoSlope(), mask, img, imgsort, imgspiral)
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

Base.eltype(polim::PolarImage{T}) where T = T
Base.length(pm::PolarImage) = length(getϕsort(pm))
Base.size(pm::PolarImage) = size(pm.cl)

function Base.show(io::IO, polim::PolarImage)
    slope = ifelse(isa(polim.slope, NoSlope), "without", "with")
    mask  = ifelse(isa(polim.mask,  NoMask ), "without", "with")
    print(io::IO, "PolarImage $slope slope and $mask mask")
end

# generic constructor for testing
# genPolarImage(M) = PolarImage(M, gencalibrate(M))

# slope(polim::PolarImage) = slope(polim.slope)
# aspect(polim::PolarImage) = aspect(polim.slope)
# has_slope(polim::PolarImage) = isa(polim.slope, Slope)
