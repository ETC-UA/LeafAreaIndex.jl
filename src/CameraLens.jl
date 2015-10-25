## CameraLens ##

"""
Type including distance to center ρ and polar angle ϕ for each pixel.
Useful for quickly accessing polar rings and segments, because of memory layout.
This is Camera + Lens dependent, but it's the same for each of its pictures,
so ρ² and ϕ only need to be (pre)calculated per camera+lens once.
We assume the lens has at least 180ᵒ field of view.
"""
immutable CameraLens
    size1::Int #number of rows of image = size(image, 1)
    size2::Int #number of columns of image = size(image, 2)
    ci::Int    #lens center is located at image[ci,cj]
    cj::Int
    fθρ::Function           #projection function θ → ρ (θ in [0,π/2], ρ in pixels)
    fρθ::Function           #inverse projection function ρ → θ
    #ρ²::Array{UInt32,2}     #squared distance to center for each point
    #ϕ::Array{Float64,2}     #azimuth angle [-π,π]
    sort_ind::Vector{Int}   #indices to sort according to ρ², then ϕ
    # Note UInt32 is sufficient to store ρ².
    spiral_ind::Vector{Int} #indices to sort according to spiral: ρ² per pixel width, then ϕ
    ρ²sort::Vector{UInt32}  #ρ² sorted by sort_ind
    ϕsort::Vector{Float64}  #ϕ sorted by sort_ind
    ρ²unique::Vector{UInt32}#unique elements of ρ²
    ρ²unique_ind::Vector{UInt32}   #(start) indices of each unique ρ² for use in ρ²sort
end

"Auxiliary function to count the number of unique ρ² and their index for CameraLens."
function count_unique(ρ²sort)
    issorted(ρ²sort) || error("ρ²sort not sorted")
    ρ²unique = unique(ρ²sort)
    ρ²uniquecounts = zeros(UInt32, length(ρ²unique)-1)
    ρ²uniquecounts[1] = 1
    prevρ² = ρ²sort[1]
    cnt = 1
    ind = 1
    for i = 2:length(ρ²sort)
        ρ² = ρ²sort[i]
        if ρ² == prevρ²
            cnt += 1
        else
            ρ²uniquecounts[ind] = cnt
            ind += 1
            cnt = 1
        end
        prevρ² = ρ²
    end
    ρ²unique, unshift!(cumsum(ρ²uniquecounts) + 1, 1)
end

"""
Spiral indices: combine ρ² for a ring of single pixel width and then sort on ϕ.
Assumes implicitely an offset of -π, where spiral jumps to next ring.
Used for gap lengths in Chen Chilar clumping.
"""
function spiralindex(ind, ρmax, ρ²unique, ρ²unique_ind, ϕsort)
    spiralind = similar(ind)
    ρ²spiralind = Int[] #temp array for ring of single pixel ρ² indices
    ind_prev = 1    
    for ρ in 1:ρmax
        singleρ² = ρ^2
        firstρ²ind = searchsortedfirst(ρ²unique, singleρ²)       
        ρ²startind = ρ²unique_ind[firstρ²ind]
        ρ²spiralind = sortperm(ArrayViews.view(ϕsort, ind_prev:ρ²startind-1))
        # TODO use @devec?
        spiralind[ind_prev:ρ²startind-1] = ind_prev - 1 + ρ²spiralind
        ind_prev = ρ²startind
    end
    # last circle of spiral
    ρ²spiralind = sortperm(ϕsort[ind_prev:end])
    spiralind[ind_prev:end] = ind_prev - 1 + ρ²spiralind
    # finally the spiral_ind is to be used on original image, not on ϕsort (or ρsort)
    spiral_ind = ind[spiralind]
end

using Base.Test
function check_calibration_inputs(size1, size2, ci::Int, cj::Int, fθρ::Function,
                                  fρθ::Function)

    size1 < 0 && error(BoundsError())
    size2 < 0 && error(BoundsError())
    ci < 0 && error(BoundsError())
    cj < 0 && error(BoundsError())

    abs(ci/size1 - 0.5) > 0.1 && warn("ci ($ci) more than 10% away from center ($(size1/2).")
    abs(cj/size2 - 0.5) > 0.1 && warn("cj ($cj) more than 10% away from center ($(size1/2)).")

    abs(ci/size1 - 0.5) > 0.2 && error("ci ($ci) more than 20% away from center ($(size1/2).")
    abs(cj/size2 - 0.5) > 0.2 && error("cj ($cj) more than 20% away from center ($(size1/2)).")

    fθρ(0.) < 0 && error("Incorrectly defined fθρ; fθρ(0) < 0.")
    fθρ(pi/2) < 2 && error("Incorrectly defined fθρ; fθρ(π/2) < 2")
    all(diff(map(fθρ, linspace(0, pi/2, 100))) .> 0) || error("Incorrectly defined fθρ; fθρ not monotonic")

    fρθ(0) < 0. && error("Incorrectly defined fρθ; fρθ(0) < 0.")
    fρθ(2) > pi/2 && error("Incorrectly defined fρθ; fρθ(2)> π/2")
    all(diff(map(fρθ, linspace(0, pi/2, 100))) .> 0) || error("Incorrectly defined fρθ; fρθ not monotonic")

    for θ in linspace(0, pi/2, 100)
        @test_approx_eq θ fρθ(fθρ(θ))
    end

    return true
end

"CameraLens constructor function."
function CameraLens(size1::Int, size2::Int, ci::Int, cj::Int, fθρ::Function, fρθ::Function)

    check_calibration_inputs(size1, size2, ci, cj, fθρ, fρθ)

    ρ² = zeros(UInt32, size1, size2)
    ϕ = zeros(Float64, size1, size2)
    for j = 1:size2
        @fastmath @inbounds @simd for i = 1:size1
            x = j - cj
            y = ci - i
            ρ²[i,j] = x^2 + y^2
            ϕ[i,j] = atan2(y, x)
        end
    end
    ϕ[ci, cj] = 0.0 #fix ϕ at center

    # sort by increasing ρ², then by increasing ϕ
    ind = sortperm(reshape(ρ² + ϕ / 10, length(ρ²))) #this takes most time
    ρ²sort = ρ²[ind]
    ρmax = round(Int, fθρ(π/2))
    ρ²indmax = searchsortedlast(ρ²sort, ρmax^2)
    ρ²sort = ρ²sort[1:ρ²indmax]
    ind = ind[1:ρ²indmax]
    ϕsort = ϕ[ind]

    # For fast indexing in ϕsort and ρ²sort by increasing ρ² we calculate the
    # indices where each unique ρ² starts. This is uesed by PolarRings to
    # iterate over rings with increasing ρ² values.
    ρ²unique, ρ²unique_ind =  count_unique(ρ²sort)

    # Calculate the spiral index
    spiral_ind = spiralindex(ind, ρmax, ρ²unique, ρ²unique_ind, ϕsort)

    CameraLens(size1,size2,ci,cj,fθρ,fρθ,ind,spiral_ind,ρ²sort,ϕsort,ρ²unique,ρ²unique_ind)
end

function Base.show(io::IO, cl::CameraLens)
    print(io::IO, "CameraLens object with:\n")
    print(io::IO, "\t size: (", cl.size1,", ",cl.size2,")\n")
    print(io::IO, "\t center i,j: ",cl.ci, ", ",cl.cj,"\n")
end
Base.writemime(io::IO, ::MIME"text/plain", cl::CameraLens) = show(io, cl)

Base.size(cl::CameraLens) = (cl.size1, cl.size2)
Base.length(cl::CameraLens) = prod(size(cl))

"Generic constructor for testing"
function CameraLens(M::AbstractMatrix)
    size1, size2 = size(M)
    ci, cj = ceil(Int, size1/2), ceil(Int, size2/2) #center point
    ρmax = min(ci, cj, size2-ci, size2-cj)
    fθρ(θ) = θ / (π/2) * ρmax
    fρθ(ρ) = ρ / ρmax * π/2
    CameraLens(size1, size2, ci, cj, fθρ, fρθ)
end
