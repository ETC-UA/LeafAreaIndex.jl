## CameraLens ##

# Type including distance to center ρ and polar angle ϕ for each pixel.
# This is Camera + Lens dependent, but it's the same for each of its pictures, 
# so ρ² and ϕ only need to be (pre)calculated per camera+lens once.
# (Note Uint32 sufficient to store ρ²).
type CameraLens
    size1::Int #number of rows in picture (heigth)
    size2::Int #number of columns in picture (width)
    ci::Int    #lens center is located at image[ci,cj]
    cj::Int
    fθρ::Function           #projection function θ → ρ (θ in [0,π/2], ρ in pixels)
    fρθ::Function           #inverse projection function ρ → θ
    #ρ²::Array{Uint32,2}     #squared distance to center for each point
    #ϕ::Array{Float64,2}     #azimuth angle [-π,π]
    sort_ind::Vector{Int}   #indices to sort according to ρ², then ϕ
    spiral_ind::Vector{Int} #indices to sort according to spiral: ρ² per pixel width, then ϕ
    ρ²sort::Vector{Uint32}  #ρ²
    ϕsort::Vector{Float64}  #ϕ sorted by sort_ind
    ρ²unique::Vector{Uint32}#unique elements of ρ²
    ρ²Ncs::Vector{Uint32}   #cumsum of occurrences of each ρ²unique in ρ²
end

# CameraLens constructor function
function calibrate(size1, size2, ci, cj, fθρ, fρθ)
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

    # sort by increasing ρ², then by increasing ϕ
    ind = sortperm(reshape(ρ²+ϕ/10, length(ρ²)))
    ρ²sort = ρ²[ind]
    ρmax = fθρ(π/2)
    ρ²indmax = searchsortedlast(ρ²sort, ρmax^2)
    ρ²sort = ρ²sort[1:ρ²indmax]
    ind = ind[1:ρ²indmax]
    ϕsort = ϕ[ind]

    ρ²unique = unique(ρ²sort)
    # TODO check if faster to use `StatsBase.mapcount` instead of `hist`
    _,ρ²N = hist(ρ²sort, unique(ρ²sort))
    ρ²Ncs = cumsum(prepend!(ρ²N, [1])) #add center b/c falls outside first bin 

    # spiral indices: combine ρ² for single pixel width and then sort on ϕ.
    # Assumes implicitely an offset of -π, where spiral jumps to next ring.
    spiralind = similar(ind)
    ρ²spiralind = Int[]
    ind_prev = 1
    #ρ²startind = 0
    for singleρ² in [1:ρmax].^2
        #TODO improve by searching ρ²unique and taking index from ρ²Ncs
        ρ²startind = searchsortedfirst(ρ²sort, singleρ²)
        ρ²spiralind = sortperm(ArrayViews.view(ϕsort, ind_prev:ρ²startind-1))
        # TODO use @devec
        spiralind[ind_prev:ρ²startind-1] = ind_prev - 1 + ρ²spiralind        
        ind_prev = ρ²startind
    end
    # last circle of spiral
    ρ²spiralind = sortperm(ϕsort[ind_prev:end])
    spiralind[ind_prev:end] = ind_prev - 1 + ρ²spiralind  
    spiral_ind = ind[spiralind]      

    CameraLens(size1,size2,ci,cj,fθρ,fρθ,ind,spiral_ind,ρ²sort,ϕsort,ρ²unique,ρ²Ncs)
end

# generic constructor for testing
function gencalibrate(M::AbstractMatrix)
    size1, size2 = size(M)
    ci,cj = iceil(size1/2), iceil(size2/2) #center point
    ρmax = min(ci, cj, size2-ci, size2-cj) 
    fθρ(θ) = θ / (π/2) * ρmax
    fρθ(ρ) = ρ / ρmax * π/2
    CameraLens(size1, size2, ci, cj, fθρ, fρθ)
end

