# Design Notes — Options Pricing & Risk Engine

## Why three pricers

Each method fails differently. The closed form is exact but only exists
because Black–Scholes assumes lognormal dynamics; Monte Carlo makes no
closed-form assumption but is noisy; an implied-vol solver is only as good
as the pricer it inverts. Wiring them against each other turns three
error-prone components into a system that detects its own bugs:

- If the Monte Carlo drifts from the closed form beyond sampling error,
  one of them is wrong — parity and the textbook case say which.
- If the solver cannot round-trip a price the closed form produced, the
  solver (or vega) is wrong.
- If an analytic Greek disagrees with a finite difference of the pricer,
  the Greek is wrong, because the pricer is already pinned by parity.

## Numerical choices

- **Normal CDF via `std::erfc`.** The classic polynomial approximations
  lose accuracy exactly where options get interesting — deep in/out of the
  money, where N(d) is near 0 or 1. `0.5 * erfc(-x / sqrt(2))` keeps full
  double precision in both tails.
- **Terminal-price Monte Carlo.** A European payoff depends only on S_T,
  so each path is a single normal draw: S_T = S·exp((r−q−σ²/2)T + σ√T·Z).
  Time-stepping would add cost and discretization error and buy nothing.
- **Antithetic variates.** Pairing z with −z cancels the first-order
  sampling error of a monotone payoff. Measured effect in this build: the
  standard error falls to 52.6% of the plain estimator at equal draws.
- **Deterministic RNG.** splitmix64 is small, fast, and identical on every
  platform, so a Monte Carlo check that passes here passes everywhere —
  the convergence tests are deterministic assertions, not coin flips.
- **Implied vol: bracketed Newton.** Newton on vega converges in a handful
  of steps near the money, but far from the money vega → 0 and raw Newton
  jumps wildly. The solver keeps a [1e-6, 5.0] bisection bracket, accepts a
  Newton step only inside it, and bisects otherwise. Prices outside the
  no-arbitrage band (below discounted intrinsic value, above the discounted
  spot/strike bound) return NaN: no volatility produces them, and a loud
  failure beats a plausible-looking wrong answer.

## Conventions

- Vega and rho are per 1.00 unit change (100 vol points / 100 rate points),
  the raw partial derivatives — desks scale them per point at display time.
- Theta is per year of time *passing* (i.e. −∂V/∂T), the sign a trader
  means by "theta decay".
- All rates are continuously compounded; the dividend yield is continuous.

## What this is not

European exercise only, constant vol and rates, no transaction costs. An
American-options pricer (binomial tree or Longstaff–Schwartz) and a local /
stochastic vol model are the natural next steps; the verification harness
(parity-style identities + finite differences + an independent reference
model) is designed to carry over to them.
