# LeafAreaIndex

[![Build Status](https://travis-ci.org/ETC-UA/LeafAreaIndex.jl.svg?branch=master)](https://travis-ci.org/ETC-UA/LeafAreaIndex.jl)

Tools to work with [hemispherical pictures](http://en.wikipedia.org/wiki/Hemispherical_photography) for the determination of [Leaf Area Index (LAI)](http://en.wikipedia.org/wiki/Leaf_area_index).

This package introduces the PolarImage type with a few convenient methods to access certain parts of the image in polar coordinates.

## Quick introduction

You construct a PolarImage from a Calibration type and an Image (or in general, a Matrix). For the calibration you need the image size, the coordinates of the lens center and the (inverse) projection function. 
(The projection function maps polar distance ρ on the image to the zenith angle θ of the scene and is usually not linear.)

    mycameralens = calibrate(width, height, ci, cj, fθρ, fρθ)
    polarimg = PolarImage(img, mycameralens)

You can now construct an *iterator* to access a specific zenith range, eg between 30ᵒ and 60ᵒ. It will return the pixels on each ring in the range by increasing ρ² in a tuple with a vector of polar angles ϕ and a vector of corresponding pixels.

    
    for (ρ², ϕ, px) in rings(polarimg, pi/6, pi/3)
        # do something with each ρ², ϕ, px variable
    end

If you just want the pixels in a zenith range, `pixels(polarimg, pi/6, pi/3)` will return a vector with pixels very fast. A shortcut `pixels(polarimg)` is translated to `pixels(polarimg, 0, pi/2)`.

The `segments` function can further splits these ring pixels in n segments (eg. for clumping calculation). It returns a vector with n elements, each a vector with segment pixels.

Some methods for automatic thresholding and LAI determination:

    thresh = threshold(pixels(polarimg))
    gf30 = gapfraction(pixels(polarimg, pi/6 - 5/180*pi, pi/6 + 5/180*pi))
    LAI57 = zenith57(polarimg, thresh)
    LAImiller = miller(polarimg, thresh)
    LAIlang = lang(polarimg, thresh)
    seg = segments(polim, 55/180*pi, 60/180*pi, 8)
    seg_gapfr = Float64[LAI.gapfraction(seg[i], thresh) for i = 1:length(seg)]
    clump_LX = log(mean(seg_gapfr)) / mean(log(seg_gapfr))

This package is in full development. Please feel free to submit any issue through the issue tracker. Code suggestions are very much appreciated through a pull request.



