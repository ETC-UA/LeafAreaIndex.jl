## PolarImage ##

# Type to contain an image and its polar transform
type PolarImage{T, A <: AbstractMatrix}
    cl::CameraLens
    img::A              #original image
    imgsort::Vector{T}  #image sorted by ρ², then ϕ
end 

# PolarImage constructor
function PolarImage(img::AbstractMatrix, cl::CameraLens)
    imgsort = img[cl.sort_ind]
    PolarImage{eltype(img), typeof(img)}(cl, img, imgsort)
end

# generic constructor for testing
genPolarImage(M) = PolarImage(M, gencalibrate(M))

# convenience method for checking memory size of PolarImage object for testing
Base.sizeof(polim::PolarImage) = sum([sizeof(getfield(polim,n)) for n in names(polim)])
