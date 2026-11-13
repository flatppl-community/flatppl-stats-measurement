# Measurement reference — 0.1

Load `main.flatppl` to access the functions below. The library uses FlatPPL
language compatibility `0.1`; each source module also loads independently.

| Function | Result |
|---|---|
| [`independent_errors(means, sigma)`](#independent_errors) | Normalized independent Gaussian observation law |
| [`correlated_errors(means, sigma, correlation)`](#correlated_errors) | Normalized multivariate Gaussian observation law |
| [`calibration(standards, offset, gain, sigma)`](#calibration) | Normalized affine-response observation law |
| [`gaussian_below_probability(limit, mean, sigma)`](#gaussian_below_probability) | Scalar lower-event probability |
| [`gaussian_above_probability(limit, mean, sigma)`](#gaussian_above_probability) | Scalar upper-event probability |
| [`gaussian_interval_probability(lower, upper, mean, sigma)`](#gaussian_interval_probability) | Scalar interval-event probability |

## Conventions

All means, observations, limits, and calibration parameters are finite real
numbers. Every `sigma` is a strictly positive finite **standard deviation** in
the same units as its corresponding observation. Arrays have fixed lengths.
The caller supplies aligned dimensions and parameters in the stated domains;
the library performs no clipping or regularization.

The three observation laws retain all Gaussian normalization factors, including
scale and determinant terms when uncertainty is fitted. Their densities use
product Lebesgue measure. They add no prior, profile, marginalization, or
calibration fit. All independence statements are conditional on their supplied
parameters; shared fitted parameters may induce marginal dependence.

Event probabilities are dimensionless values, not likelihood objects or
point-observation densities. Use `log(p)` for a log likelihood factor, or
`flag ~ Bernoulli(p)` with an observed flag of one. A zero flag has likelihood
`1-p`. Multiple reported events can be multiplied only when their underlying
readings are conditionally independent. A rounded value and detection flag
derived from the same reading require their joint event probability instead.

<a id="independent_errors"></a>
## `independent_errors(means, sigma)`

Arguments: `means=mu[N]`, with `N>=1`, and either scalar `sigma=s>0` or an aligned
vector `sigma=s[N]`. Returns a normalized law on real vectors `[N]`.

```math
Y_i\mid\mu,s\ \overset{\mathrm{ind}}\sim N(\mu_i,s_i^2),\qquad
\log p(y)=-\frac12\sum_i\left[\left(\frac{y_i-\mu_i}{s_i}\right)^2
+\log(2\pi)+2\log s_i\right].
```

For scalar `sigma`, the same standard deviation applies to every component.
This is a response-error model conditional on the true means. It does not
integrate uncertainty on the means or model noisy explanatory variables.

<a id="correlated_errors"></a>
## `correlated_errors(means, sigma, correlation)`

Arguments: aligned `means=mu[N]`, `sigma=s[N]`, and correlation matrix
`correlation=R[N,N]`. `R` must be finite, symmetric, positive definite, and have
unit diagonal; `s_i>0`. Returns a normalized law on real vectors `[N]`.

```math
D=\mathrm{diag}(s),\qquad \Sigma=DRD,\qquad
Y\mid\mu,s,R\sim N_N(\mu,\Sigma).
```

`sigma` is always a vector here. In particular, off-diagonal covariance is
`s_i*s_j*R_ij`, not `R_ij` or its square. Zero off-diagonal correlations recover
the independent law; negative correlations retain their sign. Perfectly
correlated or singular matrices are outside the domain. To marginalize to a
subset of components, select both scales and the matching correlation submatrix.

<a id="calibration"></a>
## `calibration(standards, offset, gain, sigma)`

Arguments: finite exact standard values `standards=x[N]`, scalar `offset=a`,
scalar `gain=b`, and positive scalar or aligned vector `sigma=s`. Returns a
normalized law on response vectors `[N]`.

```math
Y_i\mid a,b,x,s\ \overset{\mathrm{ind}}\sim N(a+b x_i,s_i^2).
```

The offset has response units and the gain has response units per standard unit.
The standards are treated as exact. For correlated responses, use
`correlated_errors(offset .+ gain .* standards, sigma, correlation)`. Uncertain
standards require explicit latent true standards and their own observation law.
The gain is not restricted to be positive.

The [posterior example](../examples/playground.flatppl) uses independent priors
`offset ~ Normal(0,1)` and `gain ~ Normal(2,0.5)`. With design matrix
`X=[1,x]`, weights `W=diag(s^-2)`, prior mean `m0=(0,2)` and prior precision
`P0=diag(1,4)`, its analytic posterior is

```math
V=(P_0+X^T W X)^{-1},\qquad m=V(P_0m_0+X^TWy).
```

The supplied data give `m=(0.6921544209215443,1.9650077760497666)` and
`V=diag(0.011207970112079703,0.006998444790046656)`.

<a id="gaussian_below_probability"></a>
## `gaussian_below_probability(limit, mean, sigma)`

Arguments: scalar limit `c`, mean `mu`, and standard deviation `s>0`. Returns a
scalar probability for a latent measurement `Y ~ Normal(mu,s)`.

```math
P(Y\le c\mid\mu,s)=\Phi\!\left(\frac{c-\mu}{s}\right).
```

Use this for a reading reported only as below a known limit (left censoring).
It is not a distribution conditional on selection below that limit.

<a id="gaussian_above_probability"></a>
## `gaussian_above_probability(limit, mean, sigma)`

Arguments: scalar limit `c`, mean `mu`, and standard deviation `s>0`. Returns

```math
P(Y\ge c\mid\mu,s)=\Phi\!\left(\frac{\mu-c}{s}\right).
```

Use this for a reading reported only as above a limit (right censoring), or as
the probability in a detection-flag Bernoulli law. The implementation evaluates
the reflected CDF rather than subtracting a CDF close to one. The distinction
between strict and non-strict endpoints has zero probability for a Gaussian.

<a id="gaussian_interval_probability"></a>
## `gaussian_interval_probability(lower, upper, mean, sigma)`

Arguments: scalar finite limits `l<=u`, mean `mu`, and standard deviation `s>0`.
Returns

```math
P(l\le Y\le u\mid\mu,s)
=\Phi\!\left(\frac{u-\mu}{s}\right)
-\Phi\!\left(\frac{l-\mu}{s}\right).
```

For intervals entirely above the mean, both CDFs are reflected into the lower
tail before subtraction. This avoids subtracting two numbers close to one.
Both left and right intervals from eight to nine standard deviations are
checked against independent numerical integration of the Gaussian PDF.

A value rounded to nearest cell width `h>0`, reported as `r`, has event
likelihood `gaussian_interval_probability(r-h/2,r+h/2,mu,s)` under that stated
rounding rule. A value rounded down uses `[r,r+h]`. These are probability masses
on the rounded grid; do not divide by `h` to form their likelihood. A zero-width
interval has zero probability, and is not a point-measurement density.

## Censoring, truncation, and finite precision

For an observed numerical value known to have entered the sample only when
`Y>=c`, the appropriate selected-data density is `p(y)/P(Y>=c)` on that selected
support. The event helper supplies the denominator; it does not apply it
implicitly. Counts of missed observations or detection flags can change the
sampling design and must be modeled explicitly.

Helpers use the engine's `invprobit` (standard-normal CDF). Tail accuracy is
limited by that implementation and the backend's floating-point type. Very
narrow intervals can lose relative precision when two CDF values are nearly
equal, even after reflection; very small probabilities can underflow to zero.
The library supplies probabilities rather than claiming a stable log-CDF or
log-interval implementation. Taking `log` of an underflowed zero produces
negative infinity. No probability floor or epsilon is added.

## Language basis and validation

The constructors follow the FlatPPL specification's
[Normal](https://github.com/flatppl/flatppl-design/blob/main/docs/08-distributions.md#normal)
and [MvNormal](https://github.com/flatppl/flatppl-design/blob/main/docs/08-distributions.md#mvnormal)
definitions; CDFs use
[`invprobit`](https://github.com/flatppl/flatppl-design/blob/main/docs/07-functions.md#elementary-functions).
The calibration example follows the explicit prior, forward kernel,
`likelihoodof`, and `bayesupdate` measure pipeline.

`node tests/check.cjs` compares public law densities with scalar Gaussian and
conditional-normal oracles, event masses with direct PDF quadrature, and the
example posterior with Gaussian conjugacy. It also checks the covariance of
16,384 joint draws with seed 47710 and analytical Monte Carlo standard errors.
Compiled values and derivatives can be exercised with `tests/runtime.flatppl`.
