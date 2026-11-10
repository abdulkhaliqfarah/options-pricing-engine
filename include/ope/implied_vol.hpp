#pragma once

#include <cmath>
#include <limits>

#include "ope/black_scholes.hpp"

namespace ope {

// Recovers the volatility that makes bs_price(inputs, vol) equal the target
// market price. Newton steps (using vega) are taken while they stay inside a
// bisection bracket, so the solver keeps Newton's speed near the root and
// bisection's guarantee far from it. Returns NaN when the target price is
// outside the no-arbitrage band for the given inputs -- no volatility can
// produce it, and silently returning a wrong vol would be worse than failing.
inline double implied_vol(double target_price, BsInputs in, double tol = 1e-12,
                          int max_iters = 100) {
    double t = in.time_years;
    double disc = std::exp(-in.rate * t);
    double div_disc = std::exp(-in.div_yield * t);
    double lower, upper;
    if (in.type == OptionType::Call) {
        lower = std::max(in.spot * div_disc - in.strike * disc, 0.0);
        upper = in.spot * div_disc;
    } else {
        lower = std::max(in.strike * disc - in.spot * div_disc, 0.0);
        upper = in.strike * disc;
    }
    if (!(target_price > lower) || !(target_price < upper)) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    double lo = 1e-6, hi = 5.0;
    in.vol = 0.25;  // initial guess
    for (int i = 0; i < max_iters; ++i) {
        double price = bs_price(in);
        double diff = price - target_price;
        if (std::fabs(diff) < tol) return in.vol;
        if (diff > 0.0) {
            hi = in.vol;
        } else {
            lo = in.vol;
        }
        double vega = bs_greeks(in).vega;
        double next = (vega > 1e-12) ? in.vol - diff / vega : 0.5 * (lo + hi);
        if (!(next > lo && next < hi)) next = 0.5 * (lo + hi);
        in.vol = next;
    }
    return in.vol;
}

}  // namespace ope
