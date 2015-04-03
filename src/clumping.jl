

function langxiang(polim::PolarImage, thresh, θ1::Real, θ2::Real, n::Integer)
    @checkθ1θ2
    segm = segments(polim, θ1, θ2, n)
    
    segm_gapfr = Float64[gapfraction(segm[i], thresh) for i in 1:length(segm)]
    clump_LX = log(mean(segm_gapfr)) / mean(log(segm_gapfr))
end

function langxiang45(polim::PolarImage, thresh, θ1::Real, θ2::Real)
	langxiang(polim, thresh, θ1, θ2, 8)
end
