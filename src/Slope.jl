abstract type SlopeInfo end
struct NoSlope <: SlopeInfo; end

struct SlopeParams
    # Slope inclination α in radians [0, π/3]
    α :: Float64
    # Aspect ϵ is the downhill direction relative to the geographic north,
    # in clockwise direction on a map (i.e. in counter-clockwise direction on 
    # the picture). In radians [0, 2π]. Eg pi/2 means downhill to the East.
    ϵ :: Float64 
    function SlopeParams(α::Real, ϵ::Real)
        @assert 0 <= α < pi/3 
        @assert 0 <= ϵ <= 2π 
        new(float(α), float(ϵ))
    end
end

struct Slope <: SlopeInfo
    params :: SlopeParams
    τsort :: Vector{Float64} # incidence angle τ for sloped images, sorted as imgsort
end
params(sp::SlopeParams) = sp.α, sp.ϵ
params(sl::Slope) = params(sl.params)
params(sl::NoSlope) = (0.0, 0.0)

function Slope(cl::CameraLens, sp::SlopeParams)
    α, ϵ = params(sp)

    τsort = deepcopy(cl.ϕsort)
    ρ²sort = cl.ρ²sort    
    fρ²θ(ρ²) = cl.fρθ(sqrt(ρ²)) 
    #fρ²θ = approx_fρ²θ(cl) # approximate for speedup
    
    @inbounds @simd for i in eachindex(τsort)
        θ = fρ²θ(ρ²sort[i])
        ϕ = τsort[i]
        τsort[i] = acos(cos(θ)*cos(α) + sin(θ)*cos(ϵ-ϕ)*sin(α))
    end
    Slope(sp, τsort)
end
Slope(cl::CameraLens, α::Real, ϵ::Real) = Slope(cl, SlopeParams(α, ϵ))

function Base.show(io::IO, sl::Slope) 
    α, ϵ = params(sl)
    println("Slope with ground inclination", α, " for aspect direction ", ϵ)
end

# "Fit θ(x) = a*sqrt(x) + bx + c*x^3/2 with x=ρ² for fast inverse projection function."
# function approx_fρ²θ(cl::CameraLens; N=APPROX_N)    
#     ρ² = range(0, stop=cl.ρ²sort[end], length=N)
#     A = [ρ²ᵢ^p for ρ²ᵢ in ρ², p in [0.5, 1, 1.5]]
#     a,b,c = A \ cl.fρθ.(sqrt.(ρ²))
#     #return x -> a*sqrt(x) + b*x + c*x^1.5)
#     d = c / a
#     return x -> b*x + a*sqrt(x)*(1 + d*x)
# end

# slope(sl::Slope) = sl.α
# slope(sl::NoSlope) = 0.0

# aspect(sl::Slope) = sl.ϵ
# aspect(sl::NoSlope) = 0.0

# # auxiliary methods for slope correction (España et al 2007.)
# # used in `millergroup` and `millersimple`
# function slope_adj(sl::Slope, θ::Real, ϕ::Real)
#     # the cos(sl.α) term adjust for cartographic surface projection
#      1 / (cos(sl.α) * (1 + tan(θ)*cos(sl.ϵ-ϕ)*tan(sl.α)))
# end
# slope_adj(sl::Slope, θ::Real, ϕv::AbstractArray) = [slope_adj(sl, θ, ϕ) for ϕ in ϕv]
# slope_adj(sl::NoSlope, θ::Real, ϕ::Real) = 1
# slope_adj(sl::NoSlope, θ::Real, ϕv::AbstractArray) = ones(length(ϕv))


# function create_τsort(sl::NoSlope, cl::CameraLens)
#     # create θ in case of no slope
#     τsort  = similar(cl.ϕsort)
#     ρ²sort = cl.ρ²sort 
#     fρ²θ   = approx_fρ²θ(cl)
#     @inbounds for i = 1:length(τsort)
#         τsort[i] = fρ²θ(τsort[i])
#     end
#     τsort
# end