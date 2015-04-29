# LAI inversion

## Model

We estimate LAI from the gap fraction through inversion of the Poisson model:

$$ T(\theta_V) = \exp\big(-k(\theta_V) L_e\big) = \exp\Big(-\frac{G(\theta_V)L_e}{\cos\theta_V}\Big)$$
with $T(\theta_V)$: the mean __t__ransmission or gap fraction,

$\theta_V$ the __v__iew zenit angle,

$L_e$: the __e__ffective __L__eaf Area Index with $L_e = L*\Omega$ and $\Omega$ the cluming factor/correction,

$k(\theta_V)$: the extinction coefficient that depends on

$G(\theta_V)$: the projection function, which in its turn depends on a leaf angle distribution.

