module LeafAreaIndex

import FixedPointNumbers, StreamStats, ArrayViews, Optim, Memoize, LsqFit

VERSION < v"0.4-" && using Docile
@docstrings

export calibrate, PolarImage, pixels, gapfraction, contactfreq, threshold,
        zenith57, miller, lang, langxiang45

macro checkθ1θ2()
    quote θ1 < 0 && throw(DomainError())
        θ2 > π/2 && throw(DomainError())
        θ2 < θ1 && error("θ2 < θ1")
    end
end

include("CameraLens.jl")
include("Slope.jl")
include("PolarImage.jl")
include("PolarRings.jl")
include("PolarPixels.jl")
include("PolarSegments.jl")
include("thresholding.jl")
include("gapfraction.jl")
include("inversion.jl")
include("clumping.jl")
include("ChenCihlar.jl")

end # module
