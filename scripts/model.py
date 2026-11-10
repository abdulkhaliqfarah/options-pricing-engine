#!/usr/bin/env python3
"""Differential test: an independent Python reference model vs the C++ engine.

The Python model is written straight from the Black-Scholes formulas using
only the standard library. Greeks are verified by finite differences of the
Python pricer -- not by re-implementing the analytic Greek formulas -- and
the implied vol is recovered with a plain bisection solver, a different
algorithm from the C++ Newton/bisection hybrid. Monte Carlo is checked by
requiring the C++ estimate to sit within 4 of its own reported standard
errors of the exact Python price, and the Python MC to do the same.

Usage: python3 scripts/model.py <path-to-dump_cases-binary>
Exit code 0 = agreement on every case.
"""

import csv
import io
import math
import random
import subprocess
import sys


def norm_cdf(x):
    return 0.5 * math.erfc(-x / math.sqrt(2.0))


def bs_price(spot, strike, t, rate, div, vol, is_call):
    if t <= 0.0:
        return max(spot - strike, 0.0) if is_call else max(strike - spot, 0.0)
    vol_sqrt_t = vol * math.sqrt(t)
    d1 = (math.log(spot / strike) + (rate - div + 0.5 * vol * vol) * t) / vol_sqrt_t
    d2 = d1 - vol_sqrt_t
    disc = math.exp(-rate * t)
    div_disc = math.exp(-div * t)
    if is_call:
        return spot * div_disc * norm_cdf(d1) - strike * disc * norm_cdf(d2)
    return strike * disc * norm_cdf(-d2) - spot * div_disc * norm_cdf(-d1)


def fd_greek(spot, strike, t, rate, div, vol, is_call, which):
    """Central finite difference of the Python pricer, per input unit."""
    p = dict(spot=spot, strike=strike, t=t, rate=rate, div=div, vol=vol)
    steps = {"spot": 1e-4 * spot, "vol": 1e-4, "t": 1e-4, "rate": 1e-4}
    h = steps[which]
    up = dict(p)
    dn = dict(p)
    up[which] += h
    dn[which] -= h
    f = lambda d: bs_price(d["spot"], d["strike"], d["t"], d["rate"], d["div"], d["vol"], is_call)
    diff = (f(up) - f(dn)) / (2 * h)
    if which == "t":
        diff = -diff  # theta is per year passing
    return diff


def fd_gamma(spot, strike, t, rate, div, vol, is_call):
    h = 1e-4 * spot
    f = lambda s: bs_price(s, strike, t, rate, div, vol, is_call)
    return (f(spot + h) - 2 * f(spot) + f(spot - h)) / (h * h)


def implied_vol_bisect(target, spot, strike, t, rate, div, is_call):
    lo, hi = 1e-6, 5.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if bs_price(spot, strike, t, rate, div, mid, is_call) > target:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def mc_price_py(spot, strike, t, rate, div, vol, is_call, draws=60000, seed=7):
    rng = random.Random(seed)
    drift = (rate - div - 0.5 * vol * vol) * t
    diffusion = vol * math.sqrt(t)
    disc = math.exp(-rate * t)
    total = 0.0
    total_sq = 0.0
    for _ in range(draws):
        z = rng.gauss(0.0, 1.0)
        st = spot * math.exp(drift + diffusion * z)
        payoff = max(st - strike, 0.0) if is_call else max(strike - st, 0.0)
        # antithetic pair: reuse z with -z
        st2 = spot * math.exp(drift - diffusion * z)
        payoff2 = max(st2 - strike, 0.0) if is_call else max(strike - st2, 0.0)
        sample = disc * 0.5 * (payoff + payoff2)
        total += sample
        total_sq += sample * sample
    mean = total / draws
    var = max(total_sq / draws - mean * mean, 0.0)
    return mean, math.sqrt(var / draws)


def main():
    if len(sys.argv) < 2:
        print("usage: model.py <dump_cases-binary>")
        return 2
    out = subprocess.run([sys.argv[1]], capture_output=True, text=True, check=True).stdout
    rows = list(csv.DictReader(io.StringIO(out)))
    failures = 0
    checks = 0

    def check(cond, label):
        nonlocal failures, checks
        checks += 1
        if not cond:
            failures += 1
            print(f"FAIL {label}")

    for row in rows:
        is_call = row["type"] == "call"
        s = float(row["spot"]); k = float(row["strike"]); t = float(row["t"])
        r = float(row["rate"]); q = float(row["div"]); v = float(row["vol"])
        label = f'{row["type"]} S={s} K={k} T={t} vol={v}'

        exact = bs_price(s, k, t, r, q, v, is_call)
        check(abs(float(row["price"]) - exact) <= 1e-9 * max(1.0, abs(exact)),
              f"{label} price")
        check(abs(float(row["delta"]) - fd_greek(s, k, t, r, q, v, is_call, "spot")) <= 1e-5,
              f"{label} delta")
        check(abs(float(row["gamma"]) - fd_gamma(s, k, t, r, q, v, is_call)) <= 1e-5,
              f"{label} gamma")
        check(abs(float(row["vega"]) - fd_greek(s, k, t, r, q, v, is_call, "vol")) <= 1e-4,
              f"{label} vega")
        check(abs(float(row["theta"]) - fd_greek(s, k, t, r, q, v, is_call, "t")) <= 1e-3,
              f"{label} theta")
        check(abs(float(row["rho"]) - fd_greek(s, k, t, r, q, v, is_call, "rate")) <= 1e-4,
              f"{label} rho")
        cpp_se = float(row["mc_se"])
        check(abs(float(row["mc_price"]) - exact) < 4.0 * cpp_se,
              f"{label} cpp MC within 4 SE")
        py_mc, py_se = mc_price_py(s, k, t, r, q, v, is_call)
        check(abs(py_mc - exact) < 4.0 * py_se, f"{label} python MC within 4 SE")
        iv_py = implied_vol_bisect(exact, s, k, t, r, q, is_call)
        check(abs(float(row["iv"]) - iv_py) <= 1e-8, f"{label} implied vol")
        check(abs(float(row["iv"]) - v) <= 1e-8, f"{label} iv recovers vol")

    print(f"{checks} differential checks, {failures} failures "
          f"across {len(rows)} cases")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
