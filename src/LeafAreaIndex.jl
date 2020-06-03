module LeafAreaIndex

using Optim: optimize, minimizer
using StatsBase: mean, midpoints
import FileIO, JLD2, Netpbm, Parameters, StatsBase, Optim

# import StatsBase #for `midpoints`
import ColorTypes #in crops.jl for AbstractRGB

# in calibration.jl
using Graphics: Point, norm
using DataFrames

export rawblueread, rawcolourread,
    CameraLens, CameraLensParams, PolarImage,
    Slope, SlopeParams, Mask, MaskParams,
    RedMax,#pixels, gapfraction, contactfreqs,
    threshold, #EdgeDetection, MinimumThreshold, RidlerCalvard,
    inverse, #Zenith57, Miller, Lang, EllipsLUT, EllipsOpt,
    langxiang #chencihlar,
    #calibrate_center, calibrate_projfun

include("Types.jl")
include("thresholding.jl")
include("segments.jl")
include("gapfraction.jl")
include("inversion.jl")
include("clumping.jl")
include("crops.jl")
include("ChenCihlar.jl")
include("calibration.jl")
include("io.jl")

# reoccuring argument checks
function checkθ1θ2(θ1, θ2)
    θ1 < 0   && throw(DomainError())
    θ2 > π/2 && throw(DomainError())
    @assert θ1 < θ2
end

#TODO replace with OnlineStats.jl
"Specialed type for immutable streaming sum, based on PR #18 from StreamStats.jl"
struct StreamMean
    streamsum::Float64
    len::Int
end
StreamMean() = StreamMean(0.0, 0)
update(sm::StreamMean, term) = StreamMean(sm.streamsum + term, sm.len + 1)
StatsBase.mean(sm::StreamMean) = sm.streamsum / sm.len
Base.empty!(sm::StreamMean) = StreamMean()
Base.length(sm::StreamMean) = sm.len

# A fast histogram method derived from julialang PR #8952. Used in thresholding
# and slope adjustment.
function fasthist(img::AbstractVector, edg::AbstractRange)
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

end # module
