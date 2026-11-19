# flatppl-stats-measurement

[Function reference](docs/reference.md) · [Live example](https://live.flatppl.org/#model=https%3A%2F%2Fraw.githubusercontent.com%2Fflatppl-community%2Fflatppl-stats-measurement%2Fmain%2Fexamples%2Fplayground.flatppl&target=Pi_post)

Library **0.1** provides independent and correlated Gaussian observations and affine calibration:

$$Y_i\sim N(\mu_i,\sigma_i^2),\qquad Y\sim N_n(\mu,DRD),\qquad \mu_i=a+g x_i,\quad D=\operatorname{diag}(\sigma).$$

Means and responses have shape `[n]`; correlations have shape `[n,n]`. Standard deviations are positive and finite. Correlated errors require `sigma[n]` and a positive-definite correlation matrix $R$. The laws have support $\mathbb R^n$ and retain their scale and determinant normalization.

```flatppl
flatppl_compat = "0.1"
measurement = load_module("main.flatppl")
readings ~ measurement.calibration([-1, 0, 2], 0.7, 2, [0.2, 0.3, 0.4])
```

For rounded or censored data, event helpers return $\Phi((u-\mu)/\sigma)-\Phi((l-\mu)/\sigma)$ and one-sided masses. These are probabilities, not point densities. Priors and uncertainty on calibration standards remain explicit.

The [calibration example](examples/playground.flatppl) has two parameters and an analytic Gaussian posterior. The [coarsened-observation example](examples/coarsened.flatppl) covers rounding and detection.

See the reference for definitions, domains, numerical limits, and validation commands.
