// Independent density, quadrature, and conjugate-posterior oracles.
const assert = require('node:assert/strict');
const path = require('node:path');
const { spawnSync } = require('node:child_process');
const root = path.resolve(__dirname, '..');
const evaluate = request => {
  const result = spawnSync(process.execPath, [path.join(__dirname, 'evaluate.cjs')], {
    input: JSON.stringify(request), encoding: 'utf8', maxBuffer: 8 * 1024 * 1024,
  });
  assert.equal(result.status, 0, result.stderr);
  return JSON.parse(result.stdout);
};
const close = (actual, expected, rtol = 2e-11, atol = 2e-13) => {
  assert.ok(Math.abs(actual - expected) <= atol + rtol * Math.abs(expected),
    `${actual} differs from independent oracle ${expected}`);
};
const normalLog = (x, mu, sigma) => -0.5 * ((x - mu) / sigma) ** 2
  - Math.log(sigma) - 0.5 * Math.log(2 * Math.PI);
const expressions = [], expected = [];
const add = (expression, value) => { expressions.push(expression); expected.push(value); };

// Diagonal laws retain scale-dependent normalizers and heteroscedastic errors.
add('logdensityof(a.independent_errors([1,2],0.3),[1.1,1.8])',
  normalLog(1.1, 1, 0.3) + normalLog(1.8, 2, 0.3));
add('logdensityof(a.independent_errors([1,2],[0.3,0.5]),[1.1,1.8])',
  normalLog(1.1, 1, 0.3) + normalLog(1.8, 2, 0.5));
add('logdensityof(a.calibration([-1,0,2],0.7,2,[0.2,0.3,0.4]),[-1.2,0.5,4.8])',
  normalLog(-1.2, -1.3, 0.2) + normalLog(0.5, 0.7, 0.3) + normalLog(4.8, 4.7, 0.4));

// A 2D conditional-normal factorization is independent of matrix inversion.
for (const rho of [-0.6, 0, 0.4]) {
  const z1 = (1.1 - 1) / 0.3;
  const score = normalLog(1.1, 1, 0.3)
    + normalLog(1.8, 2 + rho * 0.5 * z1, 0.5 * Math.sqrt(1 - rho ** 2));
  add(`logdensityof(a.correlated_errors([1,2],[0.3,0.5],rowstack([[1,${rho}],[${rho},1]])),[1.1,1.8])`, score);
}

// Composite Simpson quadrature of the PDF does not call any CDF implementation.
const integral = (lo, hi) => {
  const n = 8192, h = (hi - lo) / n;
  const pdf = x => Math.exp(-0.5 * x * x) / Math.sqrt(2 * Math.PI);
  let sum = pdf(lo) + pdf(hi);
  for (let i = 1; i < n; ++i) sum += (i % 2 ? 4 : 2) * pdf(lo + i * h);
  return sum * h / 3;
};
for (const [lo, hi] of [[-1.2, 0.8], [-6, -5], [5, 6], [-9, -8], [8, 9], [0.1, 0.1001], [1, 1]]) {
  add(`a.gaussian_interval_probability(${lo},${hi},0,1)`, integral(lo, hi));
}
for (const z of [-8, -1.3, 0, 1.3, 8]) {
  add(`a.gaussian_below_probability(${z},0,1)`, integral(-14, z));
  add(`a.gaussian_above_probability(${z},0,1)`, integral(z, 14));
}
add('a.gaussian_interval_probability(-0.1,0.8,0.3,0.5)', integral(-0.8, 1));
const actual = evaluate(expressions);
actual.forEach((value, i) => close(value, expected[i], 2e-10, 1e-28));

// In the example, X'WX is diagonal by symmetric standards and uncertainties.
const standards = [-2, -1, 0, 1, 2], data = [-3.3, -1.2, 0.7, 2.6, 4.7];
const sigma = [0.3, 0.2, 0.25, 0.2, 0.3];
const p00 = 1 + sigma.reduce((s, x) => s + 1 / x ** 2, 0);
const p11 = 4 + sigma.reduce((s, x, i) => s + standards[i] ** 2 / x ** 2, 0);
const mean = [data.reduce((s, y, i) => s + y / sigma[i] ** 2, 0) / p00,
  (8 + data.reduce((s, y, i) => s + standards[i] * y / sigma[i] ** 2, 0)) / p11];
const points = [[0.6, 2], [0.4, 1.8], [0.8, 2.2]];
const source = 'flatppl_compat="0.1"\ne=load_module("../examples/playground.flatppl")\n'
  + points.map(([offset, gain], i) =>
    `p${i}=logdensityof(e.Pi_post,record(offset=${offset},gain=${gain}))`).join('\n');
const scores = evaluate({ source, bindings: ['p0', 'p1', 'p2'] });
const gaussianScore = ([a, b]) => -0.5 * (p00 * (a - mean[0]) ** 2 + p11 * (b - mean[1]) ** 2);
scores.forEach((score, i) => close(score - scores[0], gaussianScore(points[i]) - gaussianScore(points[0])));

// The coarsened example scores event masses once, without density normalization.
const eventScore = evaluate({
  source: 'flatppl_compat="0.1"\ne=load_module("../examples/coarsened.flatppl")\ns=logdensityof(e.L,record(mean=1.1))',
  bindings: ['s'],
})[0];
close(eventScore, Math.log(integral(-0.75, -0.25)) + Math.log(integral(1 / 3, 14)), 2e-10);

// Joint generation must preserve the declared covariance, including its sign.
const count = 16384, rho = -0.6, seed = 47710;
const [sampled] = evaluate({
  source: 'flatppl_compat="0.1"\na=load_module("../main.flatppl")\nx ~ a.correlated_errors([1,2],[0.3,0.5],rowstack([[1,-0.6],[-0.6,1]]))',
  bindings: ['x'], sampleCount: count,
});
const draws = Object.values(sampled.samples), sampleMeans = [0, 0];
for (let i = 0; i < count; ++i) {
  sampleMeans[0] += draws[2 * i] / count;
  sampleMeans[1] += draws[2 * i + 1] / count;
}
const sd = [0.3, 0.5];
sampleMeans.forEach((m, j) => assert.ok(Math.abs(m - (j + 1)) < 6 * sd[j] / Math.sqrt(count)));
let covariance = 0;
for (let i = 0; i < count; ++i) covariance += (draws[2*i] - 1) * (draws[2*i+1] - 2) / count;
const covarianceSE = sd[0] * sd[1] * Math.sqrt((1 + rho ** 2) / count);
assert.ok(Math.abs(covariance - rho * sd[0] * sd[1]) < 6 * covarianceSE);
console.log(JSON.stringify({
  numericChecks: actual.length, posteriorMean: mean,
  posteriorCovariance: [[1 / p00, 0], [0, 1 / p11]],
  samplingSeed: seed, sampleCount: count, sampleMeans,
  meanStandardErrors: sd.map(x => x / Math.sqrt(count)), covariance, covarianceSE,
}));
