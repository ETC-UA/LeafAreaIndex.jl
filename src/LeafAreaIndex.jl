module LeafAreaIndex

import FixedPointNumbers, ArrayViews, Optim, Memoize, LsqFit, FastAnonymous
using Compat

export calibrate, PolarImage, Slope,
    pixels, gapfraction, contactfreq,
    threshold, edge_threshold, minimum_threshold, RidlerCalvard,
    zenith57, miller, lang, ellips_LUT, ellips_opt,
    langxiang45, chencihlar,
    calibrate_center, calibrate_projfun

# Keep compatibility with v0.3 where possible
VERSION < v"0.4-" && using Docile

# reoccuring argument checks
function checkθ1θ2(θ1, θ2)
    θ1 < 0   && throw(DomainError())
    θ2 > π/2 && throw(DomainError())
    θ2 < θ1  && error("θ2 < θ1")
end

"Specialed type for immutable streaming sum, based on PR #18 from StreamStats.jl"
immutable StreamMean
    streamsum::Float64
    len::Int
end
StreamMean() = StreamMean(0.0, 0)
update(sm::StreamMean, term) = StreamMean(sm.streamsum + term, sm.len + 1)
Base.mean(sm::StreamMean) = sm.streamsum / sm.len
Base.empty!(sm::StreamMean) = StreamMean()

# A fast histogram method derived from julialang PR #8952. Used in thresholding
# and slope adjustment.
function fasthist(img::AbstractVector, edg::Range)
    n = length(edg) - 1
    histcount = zeros(Int, n)
    step(edg) <= 0 && error("step(edg) must be positive")
    for pixel in img
        f = (pixel - first(edg)) / step(edg)
        if 0 < f <= n
            histcount[iceil(f)] += 1
        end
    end
    histcount
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
include("calibration.jl")

end # module
