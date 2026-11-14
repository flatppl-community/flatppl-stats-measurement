# flatppl-stats-measurement

[Function reference](docs/reference.md) · [Live example](https://live.flatppl.org/#model=https%3A%2F%2Fraw.githubusercontent.com%2Fflatppl-community%2Fflatppl-stats-measurement%2Fmain%2Fexamples%2Fplayground.flatppl&target=Pi_post)

Gaussian measurement laws, correlated uncertainties, affine calibration, and rounded or censored observations for FlatPPL. Library version **0.1**.

```flatppl
flatppl_compat = "0.1"
measurement = load_module("main.flatppl")
readings ~ measurement.calibration([-1, 0, 2], 0.7, 2, [0.2, 0.3, 0.4])
```

Uncertainty arguments are standard deviations. Correlations are separate from marginal scales. Observation laws are normalized. Event helpers return probabilities. Priors, calibration fits, and uncertainty on standards remain explicit.

The [calibration example](examples/playground.flatppl) has two parameters and an analytic Gaussian posterior. The [coarsened-observation example](examples/coarsened.flatppl) covers rounding and detection.

See the reference for definitions, domains, numerical limits, and validation commands.
