# Calibration

The necessary calibration of the camera & lens setup requires two parts: the lens projection function and the lens center. Both are required to map the pixels of the picture to corresponding polar coordinates of the scene.

The lens projection function is assumed to have the same shape for every lens type, but different pixel radius for the 180$^o$ viewing circle of each specific lens with camera setup. This means you can simply take the parameters of the generic projection function for your lens type and multiply by the visible radius. 

The lens center is assumed different for each different lens with camera, even between same lens type and camera type. However, determining the lens center is straightforward.

To determine pixel coordinates of calibration images we recommend to use the [ImageView julia package](https://github.com/timholy/ImageView.jl) or the free software [paint.net](http://www.getpaint.net/).

## Lens Center

The easiest way to determine the lens center is the prick a few holes in the lens cap (or in an attached alumium foile on top of the lens) and use  pictures of the holes turned several times. These holes should each lie on a circle with the lens center in the middle. 

We provide an [IJulia](https://github.com/JuliaLang/IJulia.jl) notebook with the necessary straightforward calculations [here](https://github.com/ETC-UA/LeafAreaIndex.jl/blob/master/calibration/CenterCalibration.ipynb)

## Lens Projection Function

The lens projection function maps the view angle θ of the scene to the polar distance ρ to the center. The view angle θ lies in the interval [0, π/2]. The polar distance can be assumed to lie in an interval between [0, 1] for a generic lens projection function, and can then be multiplied by the maximum radius pixels of the specific camera setup. The `calibration` expects the function for the actual pixels, so multiplied by the maximum radius. The maximum radius can be easily determined by sampling a few points on the horizon of a level test picture and measure the distance to the lens center.

Also the inverse projection function is required by `calibration` method and simply maps back from radius to zenith angle.

Note that for both projection functions the calibration method you can provide any function between two scalars in the appropriate domain, not just a polynomial function.

We follow the setup from the [CAN_EYE user guide](http://www6.paca.inra.fr/can-eye/Documentation-Publications/Documentation) and again provide an [IJulia notebook](https://github.com/ETC-UA/LeafAreaIndex.jl/blob/master/calibration/ProjFunctionCalibration2.ipynb) with the implementation.