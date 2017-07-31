abstract type SlopeInfo end
struct NoSlope <: SlopeInfo; end

struct SlopeParams
     # Slope inclination α in radians [0, π/3]
    α :: Float64
    # Aspect ε is the downhill direction relative to the geographic north,
    # in clockwise direction on a map (i.e. in counter-clockwise direction on 
    # the picture). In radians [0, 2π]. Eg pi/2 means downhill to the East.
    ε :: Float64 
end
function SlopeParams(α::Real, ε::Real)
    @assert 0 <= α < pi/3 
    @assert 0 <= ε <= 2π 
    SlopeParams(α, ε)
end

struct Slope <: SlopeInfo
    params :: SlopeParams
    τsort :: Vector{Float64} # incidence angle τ for sloped images, sorted as imgsort
end


params(sp::SlopeParams) = sp.α, sp.ε
params(sl::Slope) = params(sl.params)

function Slope(cl::CameraLens, sp::SlopeParams)
    α, ε = params(sp)

    τsort = deepcopy(cl.ϕsort)
    ρ²sort = cl.ρ²sort    
    fρ²θ(ρ²) = cl.fρθ(sqrt(ρ²)) 
    #fρ²θ = approx_fρ²θ(cl) # approximate for speedup
    
    @inbounds @simd for i in eachindex(τsort)
        θ = fρ²θ(ρ²sort[i])
        ϕ = τsort[i]
        τsort[i] = acos(cos(θ)*cos(α) + sin(θ)*cos(ε-ϕ)*sin(α))
    end
    τsort
end
Slope(cl::CameraLens, α::Real, ε::Real) = Slope(cl, SlopeParams(α, ϵ))

# "Fit θ(x) = a*sqrt(x) + bx + c*x^3/2 with x=ρ² for fast inverse projection function."
# function approx_fρ²θ(cl::CameraLens; N=APPROX_N)    
#     ρ² = linspace(0, cl.ρ²sort[end], N)
#     A = [ρ²ᵢ^p for ρ²ᵢ in ρ², p in [0.5, 1, 1.5]]
#     a,b,c = A \ cl.fρθ.(sqrt.(ρ²))
#     #return x -> a*sqrt(x) + b*x + c*x^1.5)
#     d = c / a
#     return x -> b*x + a*sqrt(x)*(1 + d*x)
# end

# slope(sl::Slope) = sl.α
# slope(sl::NoSlope) = 0.0

# aspect(sl::Slope) = sl.ε
# aspect(sl::NoSlope) = 0.0

# # auxiliary methods for slope correction (España et al 2007.)
# # used in `millergroup` and `millersimple`
# function slope_adj(sl::Slope, θ::Real, ϕ::Real)
#     # the cos(sl.α) term adjust for cartographic surface projection
#      1 / (cos(sl.α) * (1 + tan(θ)*cos(sl.ε-ϕ)*tan(sl.α)))
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