#pragma once

#include <cmath>

namespace ope {

enum class OptionType { Call, Put };

// Black-Scholes-Merton inputs. Rates and volatility are annualized; T is in
// years. q is a continuous dividend yield. Prices are in spot currency units.
struct BsInputs {
    double spot;
    double strike;
    double time_years;
    double rate;       // continuously compounded risk-free rate
    double div_yield;  // continuous dividend yield
    double vol;        // annualized volatility
    OptionType type;
};

inline double norm_cdf(double x) {
    // N(x) via the complementary error function: accurate in both tails.
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

inline double norm_pdf(double x) {
    return std::exp(-0.5 * x * x) / std::sqrt(2.0 * 3.14159265358979323846);
}

struct D1D2 {
    double d1;
    double d2;
};

inline D1D2 d1_d2(const BsInputs& in) {
    double vol_sqrt_t = in.vol * std::sqrt(in.time_years);
    double d1 = (std::log(in.spot / in.strike) +
                 (in.rate - in.div_yield + 0.5 * in.vol * in.vol) * in.time_years) /
                vol_sqrt_t;
    return {d1, d1 - vol_sqrt_t};
}

// Closed-form European price.
inline double bs_price(const BsInputs& in) {
    if (in.time_years <= 0.0) {
        double intrinsic = (in.type == OptionType::Call) ? (in.spot - in.strike)
                                                         : (in.strike - in.spot);
        return intrinsic > 0.0 ? intrinsic : 0.0;
    }
    if (in.vol <= 0.0) {
        // Deterministic forward: discount the certain payoff.
        double fwd = in.spot * std::exp((in.rate - in.div_yield) * in.time_years);
        double intrinsic = (in.type == OptionType::Call) ? (fwd - in.strike)
                                                         : (in.strike - fwd);
        return intrinsic > 0.0 ? intrinsic * std::exp(-in.rate * in.time_years) : 0.0;
    }
    auto [d1, d2] = d1_d2(in);
    double disc = std::exp(-in.rate * in.time_years);
    double div_disc = std::exp(-in.div_yield * in.time_years);
    if (in.type == OptionType::Call) {
        return in.spot * div_disc * norm_cdf(d1) - in.strike * disc * norm_cdf(d2);
    }
    return in.strike * disc * norm_cdf(-d2) - in.spot * div_disc * norm_cdf(-d1);
}

struct Greeks {
    double delta;  // per 1.00 move in spot
    double gamma;  // per 1.00 move in spot, of delta
    double vega;   // per 1.00 (100 vol points) move in volatility
    double theta;  // per year of time passing
    double rho;    // per 1.00 (100 points) move in the risk-free rate
};

inline Greeks bs_greeks(const BsInputs& in) {
    auto [d1, d2] = d1_d2(in);
    double t = in.time_years;
    double sqrt_t = std::sqrt(t);
    double disc = std::exp(-in.rate * t);
    double div_disc = std::exp(-in.div_yield * t);
    double pdf = norm_pdf(d1);
    bool call = (in.type == OptionType::Call);

    Greeks g;
    g.delta = div_disc * (call ? norm_cdf(d1) : norm_cdf(d1) - 1.0);
    g.gamma = div_disc * pdf / (in.spot * in.vol * sqrt_t);
    g.vega = in.spot * div_disc * pdf * sqrt_t;
    double theta_common = -(in.spot * div_disc * pdf * in.vol) / (2.0 * sqrt_t);
    if (call) {
        g.theta = theta_common - in.rate * in.strike * disc * norm_cdf(d2) +
                  in.div_yield * in.spot * div_disc * norm_cdf(d1);
        g.rho = in.strike * t * disc * norm_cdf(d2);
    } else {
        g.theta = theta_common + in.rate * in.strike * disc * norm_cdf(-d2) -
                  in.div_yield * in.spot * div_disc * norm_cdf(-d1);
        g.rho = -in.strike * t * disc * norm_cdf(-d2);
    }
    return g;
}

}  // namespace ope
