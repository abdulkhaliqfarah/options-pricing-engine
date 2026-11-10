#pragma once

#include <cmath>
#include <cstdint>

#include "ope/black_scholes.hpp"
#include "ope/rng.hpp"

namespace ope {

struct McResult {
    double price;
    double std_error;
    std::uint64_t paths;  // terminal samples used (antithetic pairs count twice)
};

// Monte Carlo price of a European option under risk-neutral GBM. Only the
// terminal price matters for a European payoff, so each path is one normal
// draw -- no time stepping. With antithetic = true, each draw z is paired
// with -z and the pair's payoffs are averaged, which cancels the first-order
// sampling error and measurably tightens the standard error.
inline McResult mc_price(const BsInputs& in, std::uint64_t num_draws,
                         std::uint64_t seed, bool antithetic) {
    SplitMix64 rng(seed);
    double drift = (in.rate - in.div_yield - 0.5 * in.vol * in.vol) * in.time_years;
    double diffusion = in.vol * std::sqrt(in.time_years);
    double disc = std::exp(-in.rate * in.time_years);

    double sum = 0.0;
    double sum_sq = 0.0;
    std::uint64_t n = 0;  // number of independent samples averaged

    auto payoff = [&](double z) {
        double st = in.spot * std::exp(drift + diffusion * z);
        double p = (in.type == OptionType::Call) ? (st - in.strike) : (in.strike - st);
        return disc * (p > 0.0 ? p : 0.0);
    };

    if (antithetic) {
        for (std::uint64_t i = 0; i < num_draws; ++i) {
            double z = rng.next_normal();
            double sample = 0.5 * (payoff(z) + payoff(-z));
            sum += sample;
            sum_sq += sample * sample;
            ++n;
        }
    } else {
        for (std::uint64_t i = 0; i < num_draws; ++i) {
            double sample = payoff(rng.next_normal());
            sum += sample;
            sum_sq += sample * sample;
            ++n;
        }
    }

    double mean = sum / static_cast<double>(n);
    double var = (sum_sq / static_cast<double>(n) - mean * mean);
    if (var < 0.0) var = 0.0;
    McResult r;
    r.price = mean;
    r.std_error = std::sqrt(var / static_cast<double>(n));
    r.paths = antithetic ? 2 * n : n;
    return r;
}

}  // namespace ope
