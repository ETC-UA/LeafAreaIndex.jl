const LANGXIANG = 8

# langxiang(polim::PolarImage, thresh,  θ1::Real, θ2::Real) = langxiang(polim, thresh, θ1, θ2, LANGXIANG)

# specialize on slope type
function langxiang(polim::PolarImage, thresh, θ1::Real, θ2::Real, nϕ::Integer=LANGXIANG)
    checkθ1θ2(θ1,θ2)

    ## WITHOUT SLOPE ##
    if !hasslope(polim) 
        segm = segments(polim, θ1, θ2, nϕ)
        return log(mean(gapfraction.(segm, thresh))) / mean(loggapfraction.(segm, thresh))
    end

    ## WITH  SLOPE ##
    # like `segments` but also with τsort
    ind_first, ind_last = firstlastind(polim, θ1, θ2)
    segmvec = [eltype(polim)[] for i = 1:nϕ]
    τvec = [Float64[] for i = 1:nϕ]
    imgsort = polim.imgsort
    ϕsort = getϕsort(polim)
    τsort = polim.slope.τsort
    adj = nϕ/2π #adjustment to segment ϕsort
    @inbounds for ind in ind_first:ind_last
        indn = ceil(Int, (ϕsort[ind] + pi) * adj)
        push!(segmvec[indn], imgsort[ind])
        push!(τvec[indn], τsort[ind])
    end
    
    # TODO use weighted average θ instead of simple average
    θ = (θ1 + θ2) / 2
    K = [contactfreqs_iterate(segmvec[i], τvec[i], thresh, θ) for i = 1:nϕ]
    T = exp.( -K / cos(θ))
    for i = 1:nϕ
        if T[i] == 0.0            
            T[i] = 1 / length(segmvec[i])
        end
    end
    clump_LX = log.(mean(T)) / mean(log.(T))
end

# for Chen Cihlar see 'ChenCihlar.jl'