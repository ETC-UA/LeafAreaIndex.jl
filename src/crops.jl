abstract type CropMethod end
struct RedMax <: CropMethod end

function isredmax(p::T) where T <: ColorTypes.AbstractRGB
    pmax = max(p.r, p.g, p.b)
    # disallow masked and overexposed pixels, assume plants reflect better than dirt
    pmax == 1 && return false
    return p.r == pmax
end

gapfraction(pixs::AbstractArray, thresh::RedMax) = sum(isredmax(p) for p in pixs) / length(pixs)