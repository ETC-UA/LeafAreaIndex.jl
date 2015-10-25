## PolarImage ##
const APPROX_N=100

# Type to contain an image and its polar transform
type PolarImage{T, A <: AbstractMatrix}
    cl::CameraLens
    slope::SlopeInfo
    img::A              #original image
    imgsort::Vector{T}  #image sorted by ρ², then ϕ
    imgspiral::Vector{T}#image spiral sorted
    τsort::Vector{Float64}# incidence angle τ for sloped images, sorted as imgsort
end 

# PolarImage constructor
function PolarImage(img::AbstractMatrix, cl::CameraLens)
    size(img) == size(cl) || error("Size image ($(size(img))) does not correspond with size CameraLens ($(size(cl))).")
    imgsort = img[cl.sort_ind]
    imgspiral = img[cl.spiral_ind]
    τsort = create_τsort(NoSlope(), cl)
    PolarImage{eltype(img), typeof(img)}(cl, NoSlope(), img, imgsort, imgspiral, τsort)
end
function PolarImage(img::AbstractMatrix, cl::CameraLens, slope::SlopeInfo)
    size(img) == size(cl) || error("Size image ($(size(img))) does not correspond with size CameraLens ($(size(cl))).")
    imgsort = img[cl.sort_ind]
    imgspiral = img[cl.spiral_ind]
    τsort = create_τsort(slope, cl)
    PolarImage{eltype(img), typeof(img)}(cl, slope, img, imgsort, imgspiral,τsort)
end

# removed Memoize.@memoize b/c memory will not be freed. 
function create_τsort(sl::NoSlope, cl::CameraLens)
    # create θ in case of no slope
    τsort = similar(cl.ϕsort)
    ρ²sort = cl.ρ²sort 
    fρ²θ = approx_fρ²θ(cl)
    @inbounds for i = 1:length(τsort)
        τsort[i] = fρ²θ(τsort[i])
    end
    τsort
end
# removed Memoize.@memoize b/c memory will not be freed.
function create_τsort(sl::Slope, cl::CameraLens)
    τsort = deepcopy(cl.ϕsort)
    ρ²sort = cl.ρ²sort    
    #fρ²θ(ρ²) = cl.fρθ(sqrt(ρ²)) 
    fρ²θ = approx_fρ²θ(cl)
    α = sl.α
    ε = sl.ε
    @inbounds for i = 1:length(τsort)
        θ = fρ²θ(ρ²sort[i])
        ϕ = τsort[i]
        τsort[i] = acos(cos(θ)*cos(α) + sin(θ)*cos(ε-ϕ)*sin(α))
    end
    τsort
end

function approx_fρ²θ(cl::CameraLens; N=APPROX_N)    
    # fit θ(x) = a*sqrt(x) + bx + c*x^3/2 with x=ρ² for fast inverse projection function

    ρ² = linspace(0, cl.ρ²sort[end], N)
    A = [ρ²ᵢ^p for ρ²ᵢ in ρ², p in [0.5, 1, 1.5]]
    a,b,c = A \ map(cl.fρθ, sqrt(ρ²))
    
    f = FastAnonymous.@anon x -> (a*sqrt(x) + b*x + c*x^1.5)
    return f
end


slope(polim::PolarImage) = slope(polim.slope)
aspect(polim::PolarImage) = aspect(polim.slope)
has_slope(polim::PolarImage) = isa(polim.slope, Slope)

# generic constructor for testing
genPolarImage(M) = PolarImage(M, gencalibrate(M))

Base.eltype{T}(polim::PolarImage{T}) = T
Base.length(pm::PolarImage) = length(pm.cl.ρ²sort)
Base.size(pm::PolarImage) = size(pm.img)

function Base.show(io::IO, polim::PolarImage)
    with_out = ifelse(isa(polim.slope, NoSlope), "without", "with")
    print(io::IO, "PolarImage $with_out slope")
end
Base.writemime(io::IO, ::MIME"text/plain", cl::CameraLens) = show(io, cl)

