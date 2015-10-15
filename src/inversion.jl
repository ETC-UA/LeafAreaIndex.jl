# Constants

"Ring width for zenith57."
const RING_WIDTH = 5 / 180 * π
"""Number of grouped consecutive ρ² rings for Miller's approach
 (each typically with 4 or 8 pixels)."""
const MILLER_GROUPS = 10
"""Start and stop view angle for use in Lang's regression method from Weiss 
et al 2004 paragraph 2.2.2.2."""
const LANG_START = 25/180*pi
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

abstract InversionMethod
type Zenith57    <: InversionMethod end
type Lang        <: InversionMethod end
type Miller      <: InversionMethod end
type MillerGroup <: InversionMethod end
type MillerNaive <: InversionMethod end
type MillerRings <: InversionMethod end
type EllipsOpt   <: InversionMethod end
type EllipsLUT   <: InversionMethod end


## Generic function
## ----------------

"Generic inversion function. Default method is EllipsOpt."
function inverse(polim::PolarImage, thresh = threshold(polim); kwargs...)
    inverse(polim, thresh, EllipsOpt(); kwargs...)
end

function inverse(polim::PolarImage, im::InversionMethod; kwargs...)
    inverse(polim, threshold(polim), im; kwargs...)
end

"Default inversion function for set of gapfractions using EllipsOpt."
function inverse(θedges::AbstractArray, θmid::Vector{Float64}, K::Vector{Float64};
                  kwargs...)
    inverse(θedges, θmid, K, EllipsOpt(); kwargs...)
end

"Default number of rings for integration of a Polar Image."
function Nrings_def(polim::PolarImage)
    if isa(polim.slope, NoSlope)
        return iceil(polim.cl.fθρ(pi/2) / MILLER_GROUPS)
    end
    return iceil(polim.cl.fθρ(pi/2) / AZIMUTH_GROUPS / MILLER_GROUPS)    
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
              Nrings = Nrings_def(polim)*(θ2-θ1)/(pi/2), kwargs...)
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
Miller's formula assumes a constant leaf angle. We implement 3 methods:
* the naive one taking the integration for each unique ρ (`millersimple`)
* grouping consecutive ρ distances (`millergroup`)
* dividing the zenith range in N rings each with an equal number of pixels (default).
"""
function inverse(polim::PolarImage, thresh, ::Miller; kwargs...) 
    inverse(polim, thresh, MillerRings(); kwargs...)
end

"""Because the naive method takes too little pixels per iteration, it distorts
the gap fraction. It is only given for reference."""
function inverse(polim::PolarImage, thresh, ::MillerNaive; θmax = θMAX, kwargs...)
    prevθ = 0.0
    s = 0.0
    #define inverse function for ρ²
    fρ²θ(ρ²) = polim.cl.fρθ(sqrt(ρ²))

    for (ρ², ϕ, px) in rings(polim, 0, θmax)
        θ = fρ²θ(ρ²)
        dθ = θ - prevθ
        # slope adjustment should be done per ϕ group, but because we only have
        # enough pixels per iteration for a single gap fraction, we take the
        # mean.
        adj = mean(slope_adj(polim.slope, θ, ϕ))
        logP = loggapfraction(px, thresh)
        s -= logP  * cos(θ) * adj * sin(θ) * dθ
        prevθ = θ
    end
    2s
end

"""Group a number of consecutive ρ²-rings together and then integrate.
 dθ is incorrect for first and last, but the cos or sin will reduce these terms."""
function inverse(polim::PolarImage, thresh, ::MillerGroup;
                 group::Integer = MILLER_GROUPS, θmax::Real = θMAX, kwargs...)
    s = 0.
    prevθ = 0.
    count = 0
    pixs = eltype(polim)[]
    ϕs = Float64[] # keep ϕ for slope adjustment
    avgθ = StreamMean()
    # lens projection function from ρ² to θ
    fρ²θ(ρ²) = polim.cl.fρθ(sqrt(ρ²))
    for (ρ², ϕ, px) in rings(polim, 0., θmax)
        count += 1
        avgθ = update(avgθ, fρ²θ(ρ²))
        append!(ϕs, ϕ)
        append!(pixs, px)

        if count == group
            θ = mean(avgθ)
            dθ = θ - prevθ
            logP = loggapfraction(pixs, thresh)
            adj = mean(slope_adj(polim.slope, θ, ϕ))
            # TODO rewrite with contactfreqs?
            s -= logP * cos(θ) * adj * sin(θ) * dθ

            prevθ = θ
            count = 0
            empty!(pixs)
            empty!(ϕs)
            empty!(avgθ)
        end
    end
    return(2s)
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
function model_ellips(θmid::Vector{Float64}, params::Vector{Float64})
    alia, L = params
    K = Float64[L * G_ellips(θᵥ, ALIA_to_x(alia)) for θᵥ in θmid]
end

"""
Assuming an ellipsoidal leaf angle distribution (with one parameter), we use the inverse model 
to estimate the average leaf inclination angle (ALIA) and Leaf Area Index (LAI)
from the observed gap fraction per view zenith angle.

This method uses a curve fitting technique to find the optimal values for the leaf angle
distribution parameter ALIA and the LAI.
"""
function inverse(θedges::AbstractArray, θmid::Vector{Float64}, K::Vector{Float64}, ::EllipsOpt;
                 LAI_init::Float64 = DEFAULT_LAI_INIT, kwargs...)

    LAI_init == DEFAULT_LAI_INIT && warn("Default value detected for `LAI_init`, it's better to use an estimate from `zenith57`.")

    # Find an initial value for ALIA
    fitfunalia(alia) = sum((model_ellips(θmid, [alia, LAI_init]) .- K).^2)
    aliares = Optim.optimize(fitfunalia, 0.1, pi/2 - 0.1)
    ALIA_init = aliares.minimum

    try
        res = LsqFit.curve_fit(model_ellips, θmid, K, [ALIA_init, LAI_init])
        ALIA, LAI = res.param
        return LAI
    catch y
        if isa(y, DomainError)
            # in case LsqFit.curve_fit does not converge and wanders out of the
            # parameter space, use Optim.fminbox.
            # TODO register MinFinder and use minfinder
            lower = [0.1, 0.2]
            upper = [pi/2 - 0.1, 9]
            fitfun(x) = sum((K .- model_ellips(θmid,x)).^2)
            fitdf = fitdf = Optim.DifferentiableFunction(fitfun)
            res = Optim.fminbox(fitdf, [ALIA_init, LAI_init], lower, upper)
            ALIA, LAI = res.minimum
            return LAI
        else
            throw(y)
        end
    end
    error("inverse EllipsOpt did not terminate normally")
end

function inverse(polim::PolarImage, thresh, ::EllipsOpt;                     
                    Nrings = Nrings_def(polim), θmax = θMAX, kwargs...)
    
    θedges, θmid, K = contactfreqs(polim, 0.0, θmax, Nrings, thresh)
    # Find inital value for LAI for optimization
    LAI_init = zenith57(polim, thresh)
    inverse(θedges, θmid, K, EllipsOpt(); LAI_init=LAI_init)
end


## Ellipsoidal LUT
## ---------------

"Holds each element of the ellipsoidal Lookup Table method."
immutable LUTel
    alia::Float64
    LAI::Float64
    modelled::Array{Float64}
end

"""Populates a Lookup Table (Weiss 2004) with random contact frequencies for a 
range of ALIA and LAI values"""
function populateLUT(θmid::Vector{Float64}; Nlut = LUT_POINTS)
    LAI_max = 9.0
    LUT = Array(LUTel, Nlut)
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
