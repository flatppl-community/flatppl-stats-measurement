# flatppl-stats-measurement

Gaussian measurement errors, correlated uncertainties, affine calibration, and
event probabilities for rounded or censored readings in FlatPPL 0.1.

```flatppl
flatppl_compat = "0.1"
measurement = load_module("main.flatppl")
readings ~ measurement.calibration([-1, 0, 2], 0.7, 2, [0.2, 0.3, 0.4])
```

All uncertainty arguments are standard deviations. Correlations are supplied
separately from marginal scales. Constructors return normalized observation
laws; event helpers return scalar probabilities to use as likelihood factors.
No prior, calibration fit, or uncertainty on the standards is implicit.

The [calibration example](https://live.flatppl.org/#model=https%3A%2F%2Fraw.githubusercontent.com%2FBJMCox%2Fflatppl-stats-measurement%2Fmain%2Fexamples%2Fplayground.flatppl&target=Pi_post)
infers two parameters and has an analytic Gaussian posterior. Its library import
is pinned to an immutable commit. The [coarsened example](examples/coarsened.flatppl)
shows a rounding cell and a detection flag. See the
[reference](docs/reference.md) for domains, formulas, and finite-precision limits.

With a JavaScript engine checkout and Node 24 or newer:

```sh
FLATPPL_JS_DIR=/path/to/flatppl-js node tests/check.cjs
```

Tests use conditional-normal density factorization, direct PDF quadrature of
both tails, a conjugate posterior, and seeded covariance checks. The command
installs nothing. [tests/runtime.flatppl](tests/runtime.flatppl) provides explicit
inputs and deterministic outputs for compiled-backend checks.
