
abstract SlopeInfo

immutable NoSlope <: SlopeInfo; end

immutable Slope <: SlopeInfo
    # Slope inclination α in radians [0, π/3]
    α::Float64 
    # Aspect ε is the downhill direction relative to the geographic north,
    # in clockwise direction on a map, i.e. in counter-clockwise direction on 
    # the picture. In radians [0, 2π]
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