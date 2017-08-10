## CameraLens ##

"""
Type including distance to center ρ and polar angle ϕ for each pixel.
Useful for quickly accessing polar rings and segments, because of memory layout.
This is Camera + Lens dependent, but it's the same for each of its pictures,
so ρ² and ϕ only need to be (pre)calculated per camera+lens once.
We assume the lens has at least 180ᵒ field of view.
"""
struct CameraLens
    height :: Int # number of rows of image = size(image, 1)
    width :: Int # number of columns of image = size(image, 2)
    ci :: Int    # lens center is located at image[ci,cj]
    cj :: Int
    fθρ :: Function           # projection function θ → ρ (θ in [0,π/2], ρ in pixels)
    fρθ :: Function           # inverse projection function ρ → θ
    sort_ind :: Vector{Int}   # indices to sort according to ρ², then ϕ
    # Note UInt32 is sufficie nt to store ρ².
    ρ²sort :: Vector{UInt32}  # ρ² sorted by sort_ind
    ϕsort  :: Vector{Float64} # ϕ  sorted by sort_ind
    spiral_ind :: Vector{Int} # indices to sort according to spiral: first per pixel ρ², then per ϕ
end
 
"CameraLens constructor function."
function CameraLens(height::Int, width::Int, ci::Int, cj::Int, fθρ::Function, fρθ::Function)

    check_calibration_inputs(height, width, ci, cj, fθρ, fρθ)

    ρ² = zeros(UInt32, height, width)
    ϕ = zeros(Float64, height, width)
    for j = 1:width
        @fastmath @inbounds @simd for i = 1:height
            x = j - cj
            y = ci - i
            ρ²[i,j] = x^2 + y^2
            ϕ[i,j] = atan2(y, x)
        end
    end
    ϕ[ci, cj] = 0.0 #fix ϕ at center

    # sort by increasing ρ², then by increasing ϕ (ϕ ∈ [-π/2, π/2])
    ind = sortperm(reshape(ρ² + ϕ / 10, length(ρ²))) #this takes most time
    ρ²sort = ρ²[ind]
    ρmax = round(Int, fθρ(π/2))
    ρ²indmax = searchsortedlast(ρ²sort, ρmax^2)
    ρ²sort = ρ²sort[1:ρ²indmax]
    ind = ind[1:ρ²indmax]
    ϕsort = ϕ[ind]

    spiral_ind = spiralindex(ind, ρmax, ρ²sort, ϕsort)

    CameraLens(height,width,ci,cj,fθρ,fρθ,ind,ρ²sort,ϕsort,spiral_ind)
end

function firstlastind(cl::CameraLens, θ1::Real, θ2::Real)
    checkθ1θ2(θ1,θ2)
    ind_first = searchsortedfirst(cl.ρ²sort, cl.fθρ(θ1)^2)
    ind_last  = searchsortedlast( cl.ρ²sort, cl.fθρ(θ2)^2)
    ind_first, ind_last
end

function check_calibration_inputs(height, width, ci::Int, cj::Int, fθρ::Function, fρθ::Function)
    @assert height <= width
    @assert width > 0 && height > 0 && ci > 0 && cj > 0

    abs(ci/height - 0.5) > 0.1 && warn("cj ($ci) more than 10% away from center ($(height/2)).")
    abs(cj/width  - 0.5) > 0.1 && warn("ci ($cj) more than 10% away from center ($(width/2).")

    abs(ci/height - 0.5) > 0.2 && error("ci ($ci) more than 20% away from center ($(height/2).")
    abs(cj/width - 0.5) > 0.2 && error("cj ($cj) more than 20% away from center ($(width/2)).")

    @assert fθρ(0.0) >= 0
    @assert fθρ(pi/2) >= 2 
    all(diff(map(fθρ, linspace(0, pi/2, 100))) .> 0) || error("Incorrectly defined fθρ; fθρ not monotonic")

    fρθ(0) < 0. && error("Incorrectly defined fρθ; fρθ(0) < 0.")
    fρθ(2) > pi/2 && error("Incorrectly defined fρθ; fρθ(2)> π/2")
    all(diff(map(fρθ, linspace(0, pi/2, 100))) .> 0) || error("Incorrectly defined fρθ; fρθ not monotonic")

    for θ in linspace(0, pi/2, 100)
        @assert θ ≈ fρθ(fθρ(θ))
    end
end

"""
Spiral indices: combine ρ² for a ring of single pixel width and then sort on ϕ.
Assumes implicitely an offset of -π, where spiral jumps to next ring.
Used for gap lengths in Chen Chilar clumping.
"""
function spiralindex(ind, ρmax, ρ²sort, ϕsort)
    spiralind = similar(ind)
    ρ²spiralind = Int[] #temp array for ring of single pixel ρ² indices
    
    # loop per pixel to get spirals of 1 pixel width
    ind_prev = 1   
    for ρ in 1:ρmax
        ρ²startind = searchsortedfirst(ρ²sort, ρ^2)
        
        ρ²spiralind = sortperm(view(ϕsort, ind_prev:ρ²startind-1))
        spiralind[ind_prev:ρ²startind-1] = ind_prev - 1 + ρ²spiralind
        ind_prev = ρ²startind
    end
    # last circle of spiral
    ρ²spiralind = sortperm(ϕsort[ind_prev:end])
    spiralind[ind_prev:end] = ind_prev - 1 + ρ²spiralind
    # finally the spiral_ind is to be used on original image, not on ϕsort (or ρsort)
    spiral_ind = ind[spiralind]
end

function Base.show(io::IO, cl::CameraLens)
    println(io::IO, "CameraLens object with:")
    println(io::IO, "\t size:", size(cl))
    println(io::IO, "\t center i,j: ",cl.ci, ", ", cl.cj)
end

Base.size(cl::CameraLens) = (cl.height, cl.width)
Base.length(cl::CameraLens) = prod(size(cl))

"Generic constructor for testing"
function CameraLens(M::AbstractMatrix)
    height, width = size(M)
    ci, cj = ceil(Int, height/2), ceil(Int, width/2) #center point
    ρmax = min(ci, cj, width-ci, height-cj)
    fθρ(θ) = θ / (π/2) * ρmax
    fρθ(ρ) = ρ / ρmax * π/2
    CameraLens(height, width, ci, cj, fθρ, fρθ)
end
