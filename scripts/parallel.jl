# include this file @everywhere for local pic processing
using Images #first include Images to avoid race conditions
using LeafAreaIndex

function listrawfiles(dir::AbstractString; raw_ext=LeafAreaIndex.RAW_EXT)
    allfiles = readdir(dir)
    filenames = ASCIIString[]
    filepaths = ASCIIString[]

    for file in allfiles
        if lowercase(splitext(file)[2]) ∈ raw_ext
            push!(filenames, file)
            push!(filepaths, joinpath(dir, file))
        end
    end
    filenames, filepaths
end

function LeafAreaIndex.PolarImage(file::AbstractString, cl::CameraLens)
    img = rawblueread(file)
    size(img) == size(cl) || error("wrong image dimensions")
    LeafAreaIndex.isoverexposed(img) && warn("Overexposed image $file")
    
    if Images.height(img) > Images.width(img) # Rotate portrait mode.
        img["spatialorder"] = reverse(img["spatialorder"])
    end
    
    PolarImage(img, cl)
end

# Always use EllipsOpt and langxiang45 methods
type LAIresult{T<:LeafAreaIndex.ThresholdMethod}
	file::ASCIIString
    LAI::Float64
    thresh::Float64
    thresh_method::T
    clump::Float64
    overexp::Bool
end

function getLAI(file::AbstractString, cl::CameraLens)
    polim = PolarImage(file, cl)
    binars = [RidlerCalvard(), EdgeDetection(), MinimumThreshold()]
    
    oe = LeafAreaIndex.isoverexposed(polim)

    lais = LAIresult[]
    for binar in binars
    	thresh = threshold(polim, binar)
    	Ω = langxiang45(polim, thresh, 0, pi/2)
    	LAIe = inverse(polim, thresh)
    	filename = last(splitdir(file))
    	push!(lais, LAIresult(filename, LAIe / Ω, Float64(thresh), binar, Ω, oe))
	end
	lais
end

