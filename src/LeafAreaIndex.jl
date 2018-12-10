module LeafAreaIndex

using Optim: optimize, minimizer
import FileIO, Netpbm, Optim, ColorTypes

# in calibration.jl
using Graphics: Point, norm
using DataFrames

export rawblueread, rawcolourread,
    CameraLens, PolarImage, Slope, Mask,
    pixels, gapfraction, contactfreqs, RedMax,
    threshold, EdgeDetection, MinimumThreshold, RidlerCalvard,
    inverse, Zenith57, Miller, Lang, EllipsLUT, EllipsOpt,
    langxiang45, chencihlar,
    calibrate_center, calibrate_projfun

# reoccuring argument checks
function checkθ1θ2(θ1, θ2)
    θ1 < 0   && throw(DomainError())
    θ2 > π/2 && throw(DomainError())
    @assert θ1 < θ2
end

"Specialed type for immutable streaming sum, based on PR #18 from StreamStats.jl"
struct StreamMean
    streamsum::Float64
    len::Int
end
StreamMean() = StreamMean(0.0, 0)
update(sm::StreamMean, term) = StreamMean(sm.streamsum + term, sm.len + 1)
Base.mean(sm::StreamMean) = sm.streamsum / sm.len
Base.empty!(sm::StreamMean) = StreamMean()
Base.length(sm::StreamMean) = sm.len

# A fast histogram method derived from julialang PR #8952. Used in thresholding
# and slope adjustment.
function fasthist(img::AbstractVector, edg::Range)
    n = length(edg) - 1
    histcount = zeros(Int, n)
    step(edg) <= 0 && error("step(edg) must be positive")
    for pixel in img
        f = (pixel - first(edg)) / step(edg)
        if 0 < f <= n
            histcount[ceil(Int, f)] += 1
        end
    end
    histcount
end

include("CameraLens.jl")
include("Slope.jl")
include("Mask.jl")
include("PolarImage.jl")
include("thresholding.jl")
include("gapfraction.jl")
include("crops.jl")
include("inversion.jl")
include("clumping.jl")
include("ChenCihlar.jl")
include("calibration.jl")
include("io.jl")

end # module
