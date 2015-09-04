# default ring width for zenith57
const RING_WIDTH = 5/180*π
# number of grouped consecutive ρ² rings for Miller's approach
# (each typically with 4 or 8 pixels)
const MILLER_GROUPS = 10
# start and stop view angle for use in Lang's regression method from Weiss et al
# 2004 paragraph 2.2.2.2
const LANG_START = 25/180*pi
const LANG_END = 65/180*pi
# maximum viewing angle
const θMAX = π/2
# Lookup Table parameters as in paper
const LUT_POINTS = 10_000 # number of points in LUT
const LUT_NMEDIAN = 25 # number of points to sample (median) LAI from
# Fast anonymous functions do not work as constant keyword argument
#const LUT_ERROR_FUN = FastAnonymous.@anon (x->x^2) # for sum of squares error function
#const LUT_ERROR_FUN = x->x^2 # for sum of squares error function



##
## Constant angle
##

# At constant zenith view angle of 1 rad the gapfraction is almost independent
# of the leaf angle distribution.
function zenith57(polim::PolarImage, thresh; ringwidth=RING_WIDTH)
    G = 0.5
    θedges, θmid, K = contactfreqs(polim, 1-ringwidth, 1+ringwidth, 1, thresh)
    K57 = K[1]
    K57 / G 
end
zenith57(polim::PolarImage) = zenith57(polim, threshold(polim))

##
## Lang's method
##

# Approximate the contact frequency with a first order regression
lang(polim::PolarImage, thresh::Real) = lang(polim, thresh, LANG_START, LANG_END)

function lang(polim::PolarImage, thresh::Real, θ1::Real, θ2::Real)
    @checkθ1θ2
    Nrings = iceil(polim.cl.fθρ(pi/2) / MILLER_GROUPS * (θ2-θ1)/(pi/2))
    lang(polim, thresh, θ1, θ2, Nrings)
end

function lang(polim::PolarImage, thresh::Real, θ1::Real, θ2::Real, N::Integer)
    @checkθ1θ2
    θedges, θmid, K = contactfreqs(polim, θ1, θ2, N, thresh)
    2*sum(linreg(θmid, K))
end

##
## Miller's formula
##

# Miller's formula assuming constant leaf angle. We implement 3 methods:
# * the naive one taking the integration for each unique ρ
# * grouping consecutive ρ distances
# * dividing the zenith range in N rings with equal number of pixels (default)
miller(polim::PolarImage, thresh::Real) = millerrings(polim, thresh)
miller(polim::PolarImage) = miller(polim, threshold(polim))

# Because the naive method takes too little pixels per iteration, it distorts
# the gap fraction. It is only given for reference.
function millersimple(polim::PolarImage, thresh::Real; θmax=θMAX)
    prevθ = 0.
    s = 0.
    #define inverse function for ρ²
    fρ²θ(ρ²) = polim.cl.fρθ(sqrt(ρ²))

    for (ρ², ϕ, px) in rings(polim, 0, θmax)
        θ = fρ²θ(ρ²)
        dθ = θ - prevθ
        # slope adjustment should be done per ϕ group, but  because we only have
        # enough pixels per iteration for a single gap fraction, we take the
        # mean.
        adj = mean(slope_adj(polim.slope, θ, ϕ))
        logP = loggapfraction(px, thresh)
        s -= logP  * cos(θ) * adj * sin(θ) * dθ
        prevθ = θ
    end
    return(2s)
end

# Group a number of consecutive ρ²-rings together and then integrate.
# dθ is incorrect for first and last, but the cos or sin will reduce these terms.
function millergroup(polim::PolarImage, group::Integer, thresh::Real; θmax=θMAX)
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
millergroup(polim::PolarImage, thresh::Real) = millergroup(polim, MILLER_GROUPS, thresh)

# Finally, create a N rings and integrate. This is the most robust way.
function millerrings(polim::PolarImage, N::Integer, thresh; θmax=θMAX)
    θedges, θmid, K = contactfreqs(polim, 0., θmax, N, thresh)
    dθ = diff(θedges)
    2 * sum(K .* sin(θmid) .* dθ)
end

function millerrings(polim::PolarImage, thresh::Real; kwargs...)
    if isa(polim.slope, NoSlope)
        Nrings = iceil(polim.cl.fθρ(pi/2) / MILLER_GROUPS)
    else
        Nrings = iceil(polim.cl.fθρ(pi/2) / AZIMUTH_GROUPS / MILLER_GROUPS)
    end
    millerrings(polim, Nrings, thresh; kwargs...)
end


##
## Ellipsoidal
##

# Assuming an ellipsoidal leaf angle distribution, we use the inverse model to
# estimate the average leaf inclination angle (ALIA) and Leaf Area Index (LAI)
# from the observed gap fraction per view zenith angle.

# First method using optimization
function ellips_opt(polim::PolarImage, Nrings::Integer, thresh; θmax=θMAX)

    # Ellipsoidal projection function, formula (A.6) Thimonier et al. 2010
    G(θᵥ, χ) = cos(θᵥ) * sqrt(χ^2 + tan(θᵥ)^2) / (χ+1.702*(χ+1.12)^-0.708)

    # function to convert ALIA to ellipsoidal parameter χ,
    # formula (30) in Wang et al 2007
    ALIA_to_x(ALIA) = (ALIA/9.65).^-0.6061 - 3

    θedges, θmid, K = contactfreqs(polim, 0., θmax, Nrings, thresh)


    function model(θmid, params)
        alia, L = params
        Float64[L * G(θᵥ, ALIA_to_x(alia)) for θᵥ in θmid]
    end

    # starting point for optimization
    L_init = zenith57(polim, thresh)
    fitfunalia(alia) = sum((model(θmid, [alia, L_init]) .- K).^2)
    aliares = Optim.optimize(fitfunalia, 0.1, pi/2-.1)
    ALIA_init = aliares.minimum

    try
        res = LsqFit.curve_fit(model, θmid, K, [ALIA_init, L_init])
        ALIA, LAI = res.param
        return LAI
    catch y
        if isa(y, DomainError)
            # in case LsqFit.curve_fit does not converge and wanders out of the
            # parameter space, use Optim.fminbox.
            # TODO register MinFinder and use minfinder
            lower = [0.1, 0.2]
            upper = [pi/2-0.1, 9]
            fitfun(x) = sum((K .- model(θmid,x)).^2)
            fitdf = fitdf = Optim.DifferentiableFunction(fitfun)
            res = Optim.fminbox(fitdf, [ALIA_init, L_init], lower, upper)
            ALIA, LAI = res.minimum
            return LAI
        else
            throw(y)
        end
    end
    error("ellips_opt did not terminate normally")
end

function ellips_opt(polim::PolarImage, thresh::Real; θmax=θMAX)
    Nrings = iceil(polim.cl.fθρ(pi/2) / MILLER_GROUPS)
    ellips_opt(polim, Nrings, thresh; θmax=θmax)
end

# Second method using Lookup Table (LUT)
# type to hold each element of the LUT
immutable LUTel
    alia::Float64
    LAI::Float64
    modelled::Array{Float64}
end

function ellips_LUT(polim::PolarImage, Nrings::Integer, thresh; θmax=θMAX,
                    Nlut::Integer=LUT_POINTS, Nmed::Integer=LUT_NMEDIAN)

    # TODO precalculate LUT
    # TODO refactor because same begin as `ellips_opt`?

    # Ellipsoidal projection function, formula (A.6) Thimonier et al. 2010
    G(θᵥ, χ) = cos(θᵥ) * sqrt(χ^2 + tan(θᵥ)^2) / (χ+1.702*(χ+1.12)^-0.708)

    # function to convert ALIA to ellipsoidal parameter χ,
    # formula (30) in Wang et al 2007
    ALIA_to_x(ALIA) = (ALIA/9.65).^-0.6061 - 3

    θedges, θmid, K = contactfreqs(polim, 0., pi/2, Nrings, thresh)

    function model(θmid, params)
        alia, L = params
        Float64[L * G(θᵥ, ALIA_to_x(alia)) for θᵥ in θmid]
    end

    # As in paper, populate LUT randomly.
    # TODO consider using Sobol pseudorandom numbers for consistency.
    function populateLUT(Nlut, model; LAI_max = 9.)
        LUT = Array(LUTel, Nlut)
        alia_max = pi/2 - .001 #against possible instability at π/2
        for i = 1:Nlut
            alia = rand() * alia_max
            LAI = rand() * LAI_max
            modelled = model(θmid, [alia, LAI])
            LUT[i] = LUTel(alia, LAI, modelled)
        end
        LUT
    end

    LUT = populateLUT(Nlut, model)

    # create fitness function against observed contact frequencies
    LUTdiff = zeros(Float64, Nlut)
    for i = 1:length(LUT)
        # errfun as keyword argument does not accept FastAnonymous, so use SSE.
        #LUTdiff[i] = sum(map(errfun, LUT[i].modelled .- K))
        LUTdiff[i] = sum((LUT[i].modelled .- K).^2)
    end
    # get median of closest parameters
    LUTsortind = sortperm(LUTdiff)
    #LUT_alia = median([el.alia for el in LUT[LUTsortind[1:Nmed]]])
    LUT_LAI = median([el.LAI for el in LUT[LUTsortind[1:Nmed]]])
end

function ellips_LUT(polim::PolarImage, thresh::Real; kwargs...)
    Nrings = iceil(polim.cl.fθρ(pi/2) / MILLER_GROUPS)
    ellips_LUT(polim, Nrings, thresh; kwargs...)
end
