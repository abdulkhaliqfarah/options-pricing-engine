// Dumps a fixed table of pricing cases as CSV for the independent Python
// reference model (scripts/model.py) to verify: analytic price and Greeks,
// Monte Carlo price and standard error, and the implied vol recovered from
// the analytic price.
#include <cstdio>
#include <initializer_list>

#include "ope/black_scholes.hpp"
#include "ope/implied_vol.hpp"
#include "ope/monte_carlo.hpp"

using namespace ope;

int main() {
    struct Case {
        double s, k, t, r, q, vol;
    };
    const Case cases[] = {
        {100, 100, 1.0, 0.05, 0.0, 0.20},  {42, 40, 0.5, 0.10, 0.0, 0.20},
        {100, 90, 0.25, 0.03, 0.01, 0.35},  {50, 60, 2.0, 0.02, 0.03, 0.15},
        {200, 210, 0.75, 0.08, 0.0, 0.50},  {80, 120, 0.5, 0.03, 0.02, 0.30},
        {150, 140, 1.5, 0.04, 0.01, 0.25},  {60, 75, 0.75, 0.04, 0.0, 0.45},
    };
    std::printf("type,spot,strike,t,rate,div,vol,price,delta,gamma,vega,theta,"
                "rho,mc_price,mc_se,iv\n");
    std::uint64_t seed = 20261110;
    for (const auto& c : cases) {
        for (OptionType type : {OptionType::Call, OptionType::Put}) {
            BsInputs in{c.s, c.k, c.t, c.r, c.q, c.vol, type};
            Greeks g = bs_greeks(in);
            McResult mc = mc_price(in, 100000, seed, true);
            double iv = implied_vol(bs_price(in), in);
            std::printf("%s,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,"
                        "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                        type == OptionType::Call ? "call" : "put", c.s, c.k, c.t,
                        c.r, c.q, c.vol, bs_price(in), g.delta, g.gamma, g.vega,
                        g.theta, g.rho, mc.price, mc.std_error, iv);
            seed += 131;
        }
    }
    return 0;
}
