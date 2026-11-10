# Options Pricing & Risk Engine — Three Pricers That Have to Agree

A C++20 European options pricing engine built around cross-verification:
a **Black–Scholes closed form**, a **Monte Carlo pricer** with antithetic
variates, and an **implied-volatility solver** — each one checked against
the others, plus a full set of analytic Greeks (delta, gamma, vega, theta,
rho) verified by finite differences.

Any one of these implementations is a homework exercise. The point of this
repo is that they are forced to agree: Monte Carlo must converge to the
closed form, the solver must invert the closed form back to the input
volatility, and every Greek must match a numerical derivative of the price.

## Quick start

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build                       # 98 scenario checks
python3 scripts/model.py ./build/dump_cases  # differential vs Python model
./build/dump_cases                           # the case table as CSV
```

## What's here

| Component | File | Notes |
|---|---|---|
| Black–Scholes price + Greeks | `include/ope/black_scholes.hpp` | `erfc`-based normal CDF, accurate in the tails; continuous dividend yield |
| Monte Carlo | `include/ope/monte_carlo.hpp` | Risk-neutral GBM, terminal-price sampling (European payoff needs no time stepping), antithetic variates |
| Implied vol | `include/ope/implied_vol.hpp` | Newton steps on vega kept inside a bisection bracket; returns NaN for prices outside the no-arbitrage band instead of a wrong number |
| RNG | `include/ope/rng.hpp` | splitmix64 + Box–Muller: same seed, same stream, bit-for-bit reproducible runs |

## Verification

- **98 scenario checks** (`tests/test_options.cpp`): put–call parity to
  1e-10 across 6 parameter sets; every Greek against central finite
  differences of the pricer; Monte Carlo within 4 standard errors of the
  closed form on 8 configurations; implied-vol round-trips recovering vols
  from 0.05 to 1.50 to within 1e-8; the textbook benchmark case
  (S=42, K=40, T=0.5, r=10%, σ=20% → call 4.759, put 0.809).
- **160 differential checks** (`scripts/model.py`): an independent Python
  reference — prices from a separate implementation of the formulas, Greeks
  by finite differences only, implied vol by plain bisection (a different
  algorithm than the C++ hybrid solver) — agrees on all 16 dumped cases.

## Measured findings

- Monte Carlo (100k antithetic draws per case) landed within **2.46
  standard errors** of the Black–Scholes price on the worst of 16 cases —
  inside the 4-SE acceptance band everywhere, with no parameter set needing
  special handling.
- Antithetic variates cut the Monte Carlo standard error to **52.6%** of
  the plain estimator at the same draw count (S=K=100, T=1, σ=0.25) —
  roughly a 4× path saving for the same accuracy.
- The implied-vol solver recovers input volatilities to **< 1e-8** across
  the full 0.05–1.50 range and a moneyness sweep, including deep
  out-of-the-money strikes where vega is tiny and pure Newton iteration
  stalls — the bisection bracket is what carries those cases.

See `DESIGN.md` for design decisions and numerical notes.
