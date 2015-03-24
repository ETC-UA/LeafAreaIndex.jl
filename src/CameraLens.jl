## CameraLens ##

# Type including distance to center ρ and polar angle ϕ for each pixel.
# This is Camera + Lens dependent, but it's the same for each of its pictures, 
# so ρ² and ϕ only need to be (pre)calculated per camera+lens once.
# (Note Uint32 sufficient to store ρ²).
type CameraLens
    size1::Uint32 #number of rows in picture (heigth)
    size2::Uint32 #number of columns in picture (width)
    ci::Uint32    #lens center is located at image[ci,cj]
    cj::Uint32
    fθR::Function           #projection function θ->pixels (θ in [0,π/2])
    fRθ::Function           #inverse projection function pixels->θ
    #ρ²::Array{Uint32,2}     #squared distance to center for each point
    #ϕ::Array{Float64,2}     #azimuth angle [-π,π]
    sort_ind::Vector{Int}   #indices to sort according to ρ², then ϕ
    ρ²sort::Vector{Uint32}  #sorted ρ²
    ϕsort::Vector{Float64}  #sorted ϕ
    ρ²unique::Vector{Uint32}#unique elements of ρ²
    ρ²Ncs::Vector{Uint32}   #cumsum of occurrences of each ρ²unique in ρ²
end

# CameraLens constructor function
function calibrate(size1, size2, ci, cj, fθR, fRθ)
    size1 < 0 && error(BoundsError())
    size2 < 0 && error(BoundsError())
    cj < 0 && error(BoundsError())
    ci < 0 && error(BoundsError())

    ρ² = zeros(Uint32, size1, size2)
    ϕ = zeros(Float64, size1, size2)
    for j = 1:size2, i = 1:size1
        x = j - cj
        y = ci - i
        ρ²[i,j] = x^2 + y^2
        ϕ[i,j] = atan2(y,x)
    end
    ϕ[ci,cj] = 0. #fix ϕ at center

    ind = sortperm(reshape(ρ²+ϕ/10, length(ρ²)))
    ρ²sort = ρ²[ind]
    Rmax = fθR(π/2)
    indmax = searchsortedlast(ρ²sort, Rmax^2)
    ρ²sort = ρ²sort[1:indmax]
    ind = ind[1:indmax]
    ϕsort = ϕ[ind]
    
    ρ²unique = unique(ρ²sort)
    _,ρ²N = hist(ρ²sort, unique(ρ²sort))
    ρ²Ncs = cumsum(prepend!(ρ²N, [1])) #add center b/c falls outside first bin 
    CameraLens(uint32(size1),uint32(size2), uint32(ci),uint32(cj),fθR,fRθ,ind,ρ²sort,ϕsort,ρ²unique,ρ²Ncs)
end

# generic constructor for testing
function gencalibrate(M::AbstractMatrix)
    size1, size2 = size(M)
    ci,cj = iceil(size1/2), iceil(size2/2) #center point
    Rmax = min(ci, cj, size2-ci, size2-cj) 
    fθR(θ) = θ / (π/2) * Rmax
    fRθ(R) = R / Rmax * π/2
    CameraLens(size1, size2, ci, cj, fθR, fRθ)
end