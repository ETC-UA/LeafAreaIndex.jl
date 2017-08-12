# Clumping

The inversion models assumes a random distribution of the foliage that is often not observed in forest canopies. Usually the outcome of the inversion procedure is considered the *effective* LAI ``L_e`` and a clumping correction ``\Omega`` is applied to obtain the true LAI ``L = L_e / \Omega``. Often ``\Omega`` is considered dependent on the view zenith angle ``\Omega(\theta)``.

We provide two standard methods to calculate the clumping for a specific zenith range: the [Lang Xiang method](http://www.sciencedirect.com/science/article/pii/016819238690033X) and the [Chen Cihlar](http://faculty.geog.utoronto.ca/Chen/Chen's%20homepage/assets/article/Quantifying%20the%20effect%20of%20canopy-IEEE.pdf).


## Lang Xiang

The Lang Xiang method divides the zenith ring in azimuthal slices and takes the logarithmic average to calculate the correction (overline denotes the mean):
```math
\Omega_{LX} = \frac{\ln \overline T}{\overline{\ln T}}
```

Use `langxiang(polimg, thresh, θ1, θ2, nϕ)` to obtain the clumping correction for a zenith range between θ1 and θ2 with nϕ the number of azimuthal slices. 

## Chen Cihlar

The Chen Cihlar method was originally proposed not for hemisperical pictures but for direct (one-dimensional) gap fraction measurements over a line under the canopy. We adopt the method to hemisperical pictures by taking a spiral of a single pixel width, similar to a [Fermat's spiral](https://en.wikipedia.org/wiki/Fermat%27s_spiral), as a basis to calculate the gap fraction (inspiration from private communication with P. Schleppi).

The principle for the Chen Cihlar method is to remove large gaps from the measured (cumulative) gap size distribution ``F_m`` to get a reduced distribution ``F_{mr}`` until it corresponds to the theoretical (cumulative) gap size distribution ``F`` from a random canopy.

The theoretical cumulative gap size distribution ``F`` depends on the projected LAI ``L_p`` and the projected characteristic leaf width ``W_p``:

```math
F(\lambda, \theta) = \big[ 1+L_p(\theta)\lambda/W_p(\theta) \big] \exp\Big(-L_p(\theta) -L_p(\theta)\lambda/W_p(\theta) \Big)
```

where ``\lambda`` represents the gap size.

``W_p`` can be estimated from the probe probabilty ``P``, i.e. "the probability that a probe of length ``l`` falls completely within a sunfleck on the forect floor" ([Leblanc et al. 2005](http://digitalcommons.unl.edu/cgi/viewcontent.cgi?article=1005&context=nasapub)):

```math
W_p = - \frac{\ln P(0,\theta)}{\big|\partial \ln P(l, \theta) / \partial l\big|_{l=0}}
```

and 
```math
P(l) = 1/L_t \sum_{i} H(\lambda_i - l)(\lambda_i - l)
```
with ``\lambda_i`` all the measured gaps, ``H`` the [Heaviside step function](https://en.wikipedia.org/wiki/Heaviside_step_function) and ``L_t`` for normalization.

``L_p(\theta)`` is initially estimated from the gap fraction ``T(\theta) = -\ln T(\theta) = -\ln F(0,\theta)``. Then by reducing ``F_m`` to ``F_{mr}`` by excluding those gaps too improbably according the theoretical gap size distributin ``F``, a new value ``L_p`` is set by ``L_p = -\ln F_{mr}`` and this procedure is repeated a few steps until the reduced distribution matches the theoretical one.

The clumping correction is then defined as
```math
\Omega_{CC} = \frac{\ln F_m}{\ln F_{mr}}\frac{1-F_{mr}}{1-F_m}
```
