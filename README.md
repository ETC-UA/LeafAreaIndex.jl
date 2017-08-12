# LeafAreaIndex

[![Build Status](https://travis-ci.org/ETC-UA/LeafAreaIndex.jl.svg?branch=master)](https://travis-ci.org/ETC-UA/LeafAreaIndex.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/ETC-UA/LeafAreaIndex.jl?branch=master&svg=true)](https://ci.appveyor.com/project/ETCUA/LeafAreaIndex-jl/branch/master)

Tools to work with [hemispherical pictures](http://en.wikipedia.org/wiki/Hemispherical_photography) for the determination of [Leaf Area Index (LAI)](http://en.wikipedia.org/wiki/Leaf_area_index).

View the full documentation on (https://etc-ua.github.io/LeafAreaIndex.jl).

# Quick introduction

Install the package through

    Pkg.clone("https://github.com/ETC-UA/LeafAreaIndex.jl")

The basic type used by this package is a PolarImage. You construct a PolarImage from a CameraLens type and an Image (or in general, an AbstractMatrix). Note that for LAI calculations typically only the blue channel of the image is used.

You can load the image eg. with the Images package:

    using Images
    img = imread("image.jpg")
    imgblue = blue(img) #take the blue channel

or in case you have the raw image from the camera, we provide a more accurate, dedicated function to extract the pixels from the blue channel (using `dcraw` under the hood):

    using LeafAreaIndex
    imgblue = rawblueread("image.NEF")

A CameraLens type is constructed given an image size, the coordinates of the lens center and the (inverse) projection function. 
(The projection function maps polar distance ρ [in pixels] on the image to the zenith angle θ [in radians] of the scene and is usually not linear.)

The basic PolarImage type is then constructed:

    using LeafAreaIndex
    mycameralens = CameraLens(height, width, centeri, centerj, funcθρ, funcρtoθ)
    polarimg = PolarImage(imgblue, mycameralens)

The first processing step is automatical thresholding (default method Ridler Calvard):

    thresh = threshold(polarimg)

In the second step LAI is estimated through the inversion model. The default method assumes an ellipsoidal leave angle distribution and uses an optimization method.

    LAI = inverse(polarimg, thresh)

Finally, the clumping factor can be estimated with the method of Lang Xiang with 45ᵒ segments between view angles θ1 and θ2:

    clump = langxiang45(polarimg, thresh, 0, pi/2)

## Further methods

For images taken (always vertically upwards) on a domain with a *slope* of eg 10ᵒ and sloping downward to the East, you must include this information in your PolarImage with the `Slope(inclination, direction)` function:

    myslope = Slope(10/180*pi, pi/2)
    polarimg = PolarImage(imgblue, mycameralens, myslope)
    
For downward taken (crop) images, create a `mask` to cut out the photographer's shoes and use the `RedMax()` method instead of thresholding to separate soil from (green) plant material

    mymask = Mask(pi/3, -2*pi/3, -pi/3)
    polarimg = PolarImage(imgblue, mycameralens, myslope)
    LAI = inverse(polarimg, RedMax())

Besides the default Ridler Calvard method, two more automatic *thresholding*methods Edge Detection and Minimum algorithm can be used:
    
    thresh  = threshold(polarimg, RidlerCalvard())
    thresh2 = threshold(polarimg, EdgeDetection())
    thresh3 = threshold(polarimg, MinimumThreshold())

Further *LAI* estimation methods for the inversion model are available: 
    * The `EllipsLUT` also assumes an ellipsoidal leaf angle distribution, but uses a Lookup Table approach instead of optimization approach.
    * The `Zenith57` method uses a ring around the view angle of 57ᵒ (1 rad) where the ALIA influence is minimal;
    * The `Miller` method integrates several zenith rings assuming a constant leaf angle; and
    * The `Lang` method uses a first order regression on the Miller method.

    LAI  = inverse(polarimg, thresh, EllipsOpt())
    LAI2 = inverse(polarimg, thresh, EllipsLUT())
    LAI3 = inverse(polarimg, thresh, Zenith57())
    LAI4 = inverse(polarimg, thresh, Miller())
    LAI5 = inverse(polarimg, thresh, Lang())

For the *clumping* factor, besides the method from Lang & Xiang, also the method from Chen & Chilar is available:

    clump2 = chencihlar(polarimg, thresh, 0, pi/2)


## Lower level methods

Under the hood several lower level methods are used to access pixels and calculated gapfractions. We suggest to look at the code for their definition and usage.

To access the pixels in a particular zenith range, `pixels(polarimg, pi/6, pi/3)` will return a vector with pixels quickly, sorted by increasing ρ (and then by polar angles ϕ for identical ρ). A shortcut `pixels(polarimg)` is translated to `pixels(polarimg, 0, pi/2)`.

The `segments` function can further split these ring pixels in n segments (eg. for clumping calculation). It returns a vector with n elements, each again a vector with the segment pixels.

For the *gapfraction*, we suggest (see online documentation) to use the contact frequencies $K(\theta_V) = -\ln[T(\theta_v)] \cos\theta_V$ for LAI inversion calculations, with $T$ the gapfraction and $\theta_V$ the view angle. The input N determines the number of rings between view angles θ1 and θ2 for a polimg with a certain treshold. The function returns a vector with angle edges of the rings, the weighted average midpoint angle for each ring and the contact frequency for each ring.

    θedges, θmid, K = contactfreqs(polimg, θ1, θ2, N, thresh)

In case of problems or suggestion, don't hesitate to submit an issue through the issue tracker or code suggestions through a pull request.

## Documentation
View the full documentation on (https://etc-ua.github.io/LeafAreaIndex.jl).
