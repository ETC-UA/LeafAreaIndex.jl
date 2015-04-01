## PolarImage ##

# Type to contain an image and its polar transform
type PolarImage{T, A <: AbstractMatrix}
    cl::CameraLens
    img::A              #original image
    imgsort::Vector{T}  #image sorted by ρ², then ϕ
    imgspiral::Vector{T}#image spiral sorted
end 

# PolarImage constructor
function PolarImage(img::AbstractMatrix, cl::CameraLens)
    imgsort = img[cl.sort_ind]
    imgspiral = img[cl.spiral_ind]
    PolarImage{eltype(img), typeof(img)}(cl, img, imgsort, imgspiral)
end

# generic constructor for testing
genPolarImage(M) = PolarImage(M, gencalibrate(M))

Base.eltype{T}(polim::PolarImage{T}) = T
Base.length(pm::PolarImage) = pm.cl.ρ²Ncs[end]

