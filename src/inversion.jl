# This file contains the inversion procedures for LAI (and ALIA) 
# calculation. The default is the ellipsoidal optimization method.


# Constants

"Ring width for zenith57."
const RING_WIDTH = 5 / 180 * π
"""Number of grouped consecutive ρ² rings for Miller's approach
 (each typically with 4 or 8 pixels)."""
const MILLER_GROUPS = 10
"Minimum number of rings for PolarImages with slope."
const NRINGS_MIN_SLOPE = 25#10
"""Start view angle for use in Lang's regression method from Weiss 
et al 2004 paragraph 2.2.2.2."""
const LANG_START = 25/180*pi
"""Stop view angle for use in Lang's regression method from Weiss 
et al 2004 paragraph 2.2.2.2."""
const LANG_END = 65/180*pi
"Maximum viewing angle used."
const θMAX = π/2
"Number of points parameters for Lookup Table method as in Weiss et al 2004."
const LUT_POINTS = 10_000 # number of points in LUT
"Sample number for LAI result median for Lookup Table method as in Weiss et al 2004."
const LUT_NMEDIAN = 25 # number of points to sample (median) LAI from Weiss 2004
"""Default LAI starting point for `ellips_opt` method. 
Better to use a simple method such as `zenith57` or `lang` for a realistic
LAI starting point for keyword argument `LAI_init`."""
const DEFAULT_LAI_INIT = 3.0
# min LAI value for warning and box constrained optimization. Was 0.2 for upward LAI, but reduced to 0.05 for crops.
const LAI_MIN = 0.05
const LAI_MAX = 10.0

abstract type InversionMethod end
struct EllipsOpt   <: InversionMethod end
struct Zenith57    <: InversionMethod end
struct Lang        <: InversionMethod end
struct Miller      <: InversionMethod end
struct MillerRings <: InversionMethod end
struct EllipsLUT   <: InversionMethod end

## Generic function
## ----------------

"Generic inversion function. Default method is EllipsOpt."
function inverse(polim::PolarImage, thresh = threshold(polim), im::InversionMethod=EllipsOpt(); kwargs...)
    inverse(polim, thresh, im; kwargs...)
end

# "Default inversion function for set of gapfractions using EllipsOpt."
# function inverse(θedges::AbstractArray, θmid::Vector{Float64}, K::Vector{Float64};
#                   kwargs...)
#     inverse(θedges, θmid, K, EllipsOpt(); kwargs...)
# end

"Calculate the default number of rings for integration of a Polar Image. 
Because the logarithm of the gap fraction is calculated, there's a delicate trade off
between more rings for more algorithmic accuracy and more pixels per ring 
(i.e. less rings) for more accurate log gap fraction per ring."
function Nrings_def(polim::PolarImage)::Int
    N =  polim.cl.fθρ(pi/2) / MILLER_GROUPS
    # hasslope(polim) && (N = max(NRINGS_MIN_SLOPE, N / INCIDENCE_GROUPS))
    return ceil(Int, N)
end
function maxviewangle(polim::PolarImage, θmax=θMAX)
    if hasslope(polim)
        α, ε = params(polim.slope)
        return min(θmax, pi/2 - α)
    else
        return θmax
    end
end

## Constant angle
## --------------

""" At constant zenith view angle of 1 rad the gapfraction is almost independent
of the leaf angle distribution."""
function inverse(polim::PolarImage, thresh, ::Zenith57; ringwidth = RING_WIDTH, args...)
    θedges, θmid, K = contactfreqs(polim, 1 - ringwidth, 1 + ringwidth, 1, thresh)
    G = 0.5
    K[1] / G
end

function inverse(θedges::AbstractArray, θmid::Vector{Float64}, K::Vector{Float64}, ::Zenith57)
    2 * K[1] # i.e. K[1] / 0.5
end


## Lang's method
## -------------

"Approximate the contact frequency with a first order regression."
function inverse(polim::PolarImage, thresh, ::Lang; θ1::Real = LANG_START, θ2::Real = LANG_END,
              Nrings = round(Int, Nrings_def(polim)*(θ2-θ1)/(pi/2)), kwargs...)
    checkθ1θ2(θ1, θ2)
    θedges, θmid, K = contactfreqs(polim, θ1, θ2, Nrings, thresh)
    inverse(θedges, θmid, K, Lang())
end

function inverse(θedges::AbstractArray, θmid::Vector{Float64}, K::Vector{Float64}, ::Lang)
    2 * sum(linreg(θmid, K))
end

## Miller's formula
## ----------------

"""
Miller's formula assumes a constant leaf angle. We divide the zenith range 
in N rings each with an equal number of pixels.
"""
function inverse(polim::PolarImage, thresh, ::Miller; kwargs...) 
    inverse(polim, thresh, MillerRings(); kwargs...)
end

"Create a N rings and integrate. This is the most robust way for Miller's method."
function inverse(polim::PolarImage, thresh, ::MillerRings;
                     N::Integer = Nrings_def(polim), θmax = θMAX, kwargs...)
    θedges, θmid, K = contactfreqs(polim, 0, θmax, N, thresh)
    inverse(θedges, θmid, K, MillerRings())
end

function inverse(θedges::AbstractArray, θmid::Vector{Float64}, K::Vector{Float64}, ::MillerRings)
    dθ = diff(θedges)
    2 * sum(K .* sin(θmid) .* dθ)
end

## Ellipsoidal
## -----------

"Ellipsoidal projection function, formula (A.6) Thimonier et al. 2010"
G_ellips(θᵥ, χ) = cos(θᵥ) * sqrt(χ^2 + tan(θᵥ)^2) / (χ + 1.702 * (χ + 1.12)^-0.708)

"function to convert ALIA to ellipsoidal parameter χ, formula (30) in Wang et al 2007"
ALIA_to_x(ALIA) = (ALIA / 9.65).^-0.6061 - 3

"The model that links ALIA and LAI with contact frequency K"
function model_ellips(θmid::Vector, params::Vector)
    alia, LAI = params
    K = [LAI * G_ellips(θᵥ, ALIA_to_x(alia)) for θᵥ in θmid]
end

"""
Assuming an ellipsoidal leaf angle distribution (with one parameter), we use the inverse model 
to estimate the average leaf inclination angle (ALIA) and Leaf Area Index (LAI)
from the observed gap fraction per view zenith angle.

This method uses a curve fitting technique to find the optimal values for the leaf angle
distribution parameter ALIA and the LAI.
"""
function inverse(θedges::AbstractArray, θmid::Vector{Float64}, K::Vector{Float64}, ::EllipsOpt;
                 LAI_init::Float64 = DEFAULT_LAI_INIT)

    LAI_init == DEFAULT_LAI_INIT && @warn("Default value detected for `LAI_init`, it's better to use an estimate from `zenith57`.")
    if !(LAI_MIN < LAI_init < LAI_MAX)
        @warn("LAI starting point ($LAI_init) for optimization outside limits [$LAI_MIN, $LAI_MAX] and will be reset.")
        LAI_init = DEFAULT_LAI_INIT # max(min(LAI_init, 9), 0.2)
    end

    # Find an initial value for ALIA
    fitfunalia(alia) = sum((model_ellips(θmid, [alia, LAI_init]) .- K).^2)
    aliares = optimize(fitfunalia, 0.1, pi/2 - 0.1)
    ALIA_init = minimizer(aliares)

    # Optimize both ALIA and LAI at same time
    fitfun(x) = sum((K .- model_ellips(θmid, x)).^2)
    initial = [ALIA_init, LAI_init]
    fitdf = Optim.OnceDifferentiable(fitfun, initial; autodiff = :forward)

    try
        res = optimize(fitdf, initial, Optim.LBFGS())
        ALIA, LAI = minimizer(res)
        return LAI
    catch y
        if isa(y, DomainError)
            # in case LsqFit.curve_fit does not converge and wanders out of the
            # parameter space, use box constrained optimization.
            lower = [0.05, LAI_MIN]
            upper = [pi/2 - 0.05, LAI_MAX]            
            res = optimize(fitfun, lower, upper, initial, Optim.Fminbox(Optim.LBFGS()))
            ALIA, LAI = minimizer(res)
            return LAI
        else
            throw(y)
        end
    end
    error("inverse EllipsOpt did not terminate normally")
end

function inverse(polim::PolarImage, thresh, ::EllipsOpt;                     
                    Nrings = Nrings_def(polim), 
                    θmax = maxviewangle(polim),
                    Nτ = INCIDENCE_GROUPS)
    
    θedges, θmid, K = contactfreqs(polim, 0.0, θmax, Nrings, thresh;Nτ=Nτ)
    # Find inital value for LAI for optimization
    LAI_init = inverse(polim, thresh, Zenith57())
    inverse(θedges, θmid, K, EllipsOpt(); LAI_init=LAI_init)
end


## Ellipsoidal LUT
## ---------------

"Holds each element of the ellipsoidal Lookup Table method."
struct LUTel
    alia::Float64
    LAI::Float64
    modelled::Array{Float64}
end

"""Populates a Lookup Table (Weiss 2004) with random contact frequencies for a 
range of ALIA and LAI values"""
function populateLUT(θmid::Vector{Float64}; Nlut = LUT_POINTS)
    LAI_max = LAI_MAX
    LUT = Array{LUTel}(Nlut)
    alia_max = pi/2 - .001 #against possible instability at π/2
    # As in paper [Weiss2004], populate LUT randomly.
    # TODO consider using Sobol pseudorandom numbers for consistency.
    for i = 1:Nlut
        alia = rand() * alia_max
        LAI = rand() * LAI_max
        modelled = model_ellips(θmid, [alia, LAI])
        LUT[i] = LUTel(alia, LAI, modelled)
    end
    return LUT
end

"""
Assuming an ellipsoidal leaf angle distribution (with one parameter), we use the 
inverse model to estimate the average leaf inclination angle (ALIA) and Leaf Area Index (LAI)
from the observed gap fraction per view zenith angle.

This method uses a Lookup Table (LUT) for the leaf angle distribution parameter 
ALIA and the LAI, following [Weiss2004](http://www.researchgate.net/profile/Inge_Jonckheere/publication/222931516_Review_of_methods_for_in_situ_leaf_area_index_(LAI)_determination_Part_II._Estimation_of_LAI_errors_and_sampling/links/09e4150cefe5a4fea5000000.pdf).
"""
function inverse(θedges::AbstractArray, θmid::Vector{Float64}, K::Vector{Float64}, ::EllipsLUT;
                    Nlut::Int = LUT_POINTS, Nmed::Integer=LUT_NMEDIAN, kwargs...)

    LUT = populateLUT(θmid; Nlut = Nlut)

    # create fitness function against observed contact frequencies
    LUTdiff = zeros(Float64, Nlut)
    for i in eachindex(LUT)
        # errfun as keyword argument does not accept FastAnonymous, so use SSE.
        #LUTdiff[i] = sum(map(errfun, LUT[i].modelled .- K))
        LUTdiff[i] = sum((LUT[i].modelled .- K).^2)
    end
    # get median of closest parameters
    LUTsortind = sortperm(LUTdiff)
    #LUT_alia = median([el.alia for el in LUT[LUTsortind[1:Nmed]]])
    LUT_LAI = median([el.LAI for el in LUT[LUTsortind[1:Nmed]]])
end

function inverse(polim::PolarImage, thresh, ::EllipsLUT;
                    Nrings = Nrings_def(polim), θmax = θMAX, Nlut::Int = LUT_POINTS, 
                    Nmed::Integer = LUT_NMEDIAN, kwargs...)

    θedges, θmid, K = contactfreqs(polim, 0.0, θmax, Nrings, thresh)
    inverse(θedges, θmid, K, EllipsLUT(); Nlut = Nlut, Nmed = Nmed)
end
