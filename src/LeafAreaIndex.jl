module LeafAreaIndex

import FixedPointNumbers, ArrayViews, Optim, Memoize, LsqFit

VERSION < v"0.4-" && using Docile
@docstrings

export calibrate, PolarImage, pixels, gapfraction, contactfreq, threshold,
       edge_threshold, minimum_threshold, RidlerCalvard, zenith57, miller, lang, 
       langxiang45

# reoccuring argument checks
macro checkθ1θ2()
    quote θ1 < 0 && throw(DomainError())
        θ2 > π/2 && throw(DomainError())
        θ2 < θ1 && error("θ2 < θ1")
    end
end

# Specialed type for immutable streaming sum, based on PR #18 from StreamStats.jl
immutable StreamMean
    streamsum::Float64
    len::Int
end
StreamMean() = StreamMean(0.,0)
update(sm::StreamMean, term) = StreamMean(sm.streamsum + term, sm.len+1)
Base.mean(sm::StreamMean) = sm.streamsum / sm.len
Base.empty!(sm::StreamMean) = StreamMean()


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
