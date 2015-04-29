
abstract SlopeInfo

immutable NoSlope <: SlopeInfo; end

immutable Slope <: SlopeInfo
    # Slope inclination α in radians [0, π/3]
    α::Float64 
    # Aspect ε is the downhill direction relative to the geographic north,
    # in clockwise direction on a map (i.e. in counter-clockwise direction on 
    # the picture). In radians [0, 2π]. Eg pi/2 means downhill to the East.
    ε::Float64 
end

#Constructor
function Slope(α::Real, ε::Real)
    α < 0 || α > pi/3 && throw(DomainError())
    ε < 0 || ε > 2π && throw(DomainError())
    Slope(α,ε)
end

slope(sl::Slope) = sl.α
slope(sl::NoSlope) = 0.

aspect(sl::Slope) = sl.ε
aspect(sl::NoSlope) = 0.

# auxiliary methods for slope correction (España et al 2007.)
# used in `millergroup` and `millersimple`
slope_adj(sl::NoSlope, θ::Real, ϕ::Float64) = 1.
function slope_adj(sl::Slope, θ::Real, ϕ::Float64) 
    # the cos(sl.α) term adjust for cartographic surface projection
     1 / (cos(sl.α) * (1 + tan(θ)*cos(sl.ε-ϕ)*tan(sl.α)))
end
slope_adj(sl::NoSlope, θ::Real, ϕv::AbstractArray) = ones(length(ϕv))
function slope_adj(sl::Slope, θ::Real, ϕv::AbstractArray)
    adj = ones(length(ϕv))
    for i = 1:length(adj)
        adj[i] = slope_adj(sl, θ, ϕv[i])
    end
    adj
    #Float64[slope_adj(sl, θ, ϕ) for ϕ in ϕv] #allocates 
end