// Scenario checks for the options pricing engine. The checks lean on
// identities that must hold for ANY correct implementation -- put-call
// parity, finite-difference Greeks, Monte Carlo convergence to the closed
// form, and implied-vol round-trips -- plus one textbook benchmark.
// Run via ctest or directly; exits non-zero if any check fails.
#include <cmath>
#include <initializer_list>
#include <cstdio>

#include "ope/black_scholes.hpp"
#include "ope/implied_vol.hpp"
#include "ope/monte_carlo.hpp"

using namespace ope;

static int failures = 0;
static int checks = 0;

#define CHECK(cond)                                     \
    do {                                                \
        ++checks;                                       \
        if (!(cond)) {                                  \
            std::printf("FAIL %s:%d  %s\n", __FILE__,   \
                        __LINE__, #cond);               \
            ++failures;                                 \
        }                                               \
    } while (0)

#define CHECK_NEAR(a, b, tol) CHECK(std::fabs((a) - (b)) <= (tol))

static BsInputs base(double spot, double strike, double t, double r, double q,
                     double vol, OptionType type) {
    return {spot, strike, t, r, q, vol, type};
}

static void test_textbook_benchmark() {
    // Hull's classic example: S=42, K=40, T=0.5, r=10%, sigma=20%.
    // Published values: call ~ 4.76, put ~ 0.81.
    auto call = base(42, 40, 0.5, 0.10, 0.0, 0.20, OptionType::Call);
    auto put = base(42, 40, 0.5, 0.10, 0.0, 0.20, OptionType::Put);
    CHECK_NEAR(bs_price(call), 4.76, 0.01);
    CHECK_NEAR(bs_price(put), 0.81, 0.01);
}

static void test_put_call_parity() {
    // C - P = S e^{-qT} - K e^{-rT} must hold exactly for every parameter set.
    const double cases[][6] = {
        {100, 100, 1.0, 0.05, 0.0, 0.20}, {100, 90, 0.25, 0.03, 0.01, 0.35},
        {50, 60, 2.0, 0.02, 0.03, 0.15},  {200, 210, 0.75, 0.08, 0.0, 0.50},
        {42, 40, 0.5, 0.10, 0.0, 0.20},   {100, 100, 0.05, 0.01, 0.0, 0.80},
    };
    for (auto& c : cases) {
        auto call = base(c[0], c[1], c[2], c[3], c[4], c[5], OptionType::Call);
        auto put = base(c[0], c[1], c[2], c[3], c[4], c[5], OptionType::Put);
        double rhs = c[0] * std::exp(-c[4] * c[2]) - c[1] * std::exp(-c[3] * c[2]);
        CHECK_NEAR(bs_price(call) - bs_price(put), rhs, 1e-10);
    }
}

static void test_price_limits() {
    // Deep in the money, a call approaches discounted intrinsic value.
    auto deep_call = base(200, 100, 1.0, 0.05, 0.0, 0.20, OptionType::Call);
    CHECK_NEAR(bs_price(deep_call), 200 - 100 * std::exp(-0.05), 0.5);
    // Deep out of the money, both prices approach zero but stay non-negative.
    auto far_call = base(50, 200, 0.25, 0.03, 0.0, 0.20, OptionType::Call);
    auto far_put = base(200, 50, 0.25, 0.03, 0.0, 0.20, OptionType::Put);
    CHECK(bs_price(far_call) >= 0.0 && bs_price(far_call) < 0.01);
    CHECK(bs_price(far_put) >= 0.0 && bs_price(far_put) < 0.01);
    // At expiry the pricer returns intrinsic value.
    auto expiry = base(120, 100, 0.0, 0.05, 0.0, 0.20, OptionType::Call);
    CHECK_NEAR(bs_price(expiry), 20.0, 1e-12);
    // Higher volatility never lowers a European price.
    auto lo = base(100, 100, 1.0, 0.05, 0.0, 0.10, OptionType::Call);
    auto hi = base(100, 100, 1.0, 0.05, 0.0, 0.40, OptionType::Call);
    CHECK(bs_price(hi) > bs_price(lo));
}

static void test_greeks_finite_difference() {
    // Every analytic Greek must match a central finite difference of the
    // closed-form pricer itself. This catches sign and scaling errors that
    // parity alone cannot see.
    const double cases[][6] = {
        {100, 100, 1.0, 0.05, 0.0, 0.20}, {80, 120, 0.5, 0.03, 0.02, 0.30},
        {150, 140, 1.5, 0.04, 0.01, 0.25},
    };
    for (auto& c : cases) {
        for (OptionType type : {OptionType::Call, OptionType::Put}) {
            BsInputs in = base(c[0], c[1], c[2], c[3], c[4], c[5], type);
            Greeks g = bs_greeks(in);
            auto shifted = in;
            double h = 1e-4 * c[0];
            shifted.spot += h;
            double up = bs_price(shifted);
            shifted.spot -= 2 * h;
            double dn = bs_price(shifted);
            double fd_delta = (up - dn) / (2 * h);
            double fd_gamma = (up - 2 * bs_price(in) + dn) / (h * h);
            CHECK_NEAR(g.delta, fd_delta, 1e-5);
            CHECK_NEAR(g.gamma, fd_gamma, 1e-5);
            shifted = in;
            shifted.vol += 1e-4;
            double vup = bs_price(shifted);
            shifted.vol -= 2e-4;
            double vdn = bs_price(shifted);
            CHECK_NEAR(g.vega, (vup - vdn) / 2e-4, 1e-4);
            shifted = in;
            shifted.time_years += 1e-4;
            double tup = bs_price(shifted);
            shifted.time_years -= 2e-4;
            double tdn = bs_price(shifted);
            // Theta is the value change per year PASSING, i.e. -dV/dT.
            CHECK_NEAR(g.theta, -(tup - tdn) / 2e-4, 1e-3);
            shifted = in;
            shifted.rate += 1e-4;
            double rup = bs_price(shifted);
            shifted.rate -= 2e-4;
            double rdn = bs_price(shifted);
            CHECK_NEAR(g.rho, (rup - rdn) / 2e-4, 1e-4);
        }
    }
}

static void test_greek_signs_and_ranges() {
    auto call = base(100, 100, 1.0, 0.05, 0.0, 0.20, OptionType::Call);
    auto put = base(100, 100, 1.0, 0.05, 0.0, 0.20, OptionType::Put);
    Greeks gc = bs_greeks(call);
    Greeks gp = bs_greeks(put);
    CHECK(gc.delta > 0.0 && gc.delta < 1.0);
    CHECK(gp.delta < 0.0 && gp.delta > -1.0);
    CHECK_NEAR(gc.delta - gp.delta, 1.0, 1e-12);  // q = 0
    CHECK(gc.gamma > 0.0 && gp.gamma > 0.0);
    CHECK_NEAR(gc.gamma, gp.gamma, 1e-15);
    CHECK(gc.vega > 0.0);
    CHECK_NEAR(gc.vega, gp.vega, 1e-12);
    CHECK(gc.rho > 0.0 && gp.rho < 0.0);
}

static void test_monte_carlo_convergence() {
    // The Monte Carlo estimate must land within 4 standard errors of the
    // closed form across diverse parameter sets -- with a fixed seed this is
    // a deterministic check, not a probabilistic one.
    const double cases[][6] = {
        {100, 100, 1.0, 0.05, 0.0, 0.20}, {42, 40, 0.5, 0.10, 0.0, 0.20},
        {100, 80, 2.0, 0.03, 0.02, 0.35},  {60, 75, 0.75, 0.04, 0.0, 0.45},
    };
    std::uint64_t seed = 12345;
    for (auto& c : cases) {
        for (OptionType type : {OptionType::Call, OptionType::Put}) {
            BsInputs in = base(c[0], c[1], c[2], c[3], c[4], c[5], type);
            McResult mc = mc_price(in, 100000, seed, true);
            double exact = bs_price(in);
            CHECK(std::fabs(mc.price - exact) < 4.0 * mc.std_error);
            CHECK(mc.std_error > 0.0);
            CHECK(mc.paths == 200000);
            seed += 777;
        }
    }
}

static void test_antithetic_reduces_error() {
    // Same number of draws: the antithetic estimator's standard error must
    // not exceed the plain estimator's (they use different seeds, so allow a
    // small margin, but antithetic pairing should clearly help here).
    auto in = base(100, 100, 1.0, 0.05, 0.0, 0.25, OptionType::Call);
    McResult plain = mc_price(in, 50000, 999, false);
    McResult anti = mc_price(in, 50000, 4242, true);
    CHECK(anti.std_error < plain.std_error * 1.05);
    CHECK(std::fabs(plain.price - bs_price(in)) < 4.0 * plain.std_error);
}

static void test_implied_vol_round_trip() {
    // Price with a known vol, then solve: the solver must recover the vol.
    const double vols[] = {0.05, 0.12, 0.20, 0.35, 0.60, 1.00, 1.50};
    for (double v : vols) {
        for (OptionType type : {OptionType::Call, OptionType::Put}) {
            BsInputs in = base(100, 95, 0.75, 0.04, 0.01, v, type);
            double price = bs_price(in);
            double recovered = implied_vol(price, in);
            CHECK_NEAR(recovered, v, 1e-8);
        }
    }
    // Moneyness sweep at a fixed vol.
    for (double k : {70.0, 85.0, 100.0, 115.0, 130.0}) {
        BsInputs in = base(100, k, 1.0, 0.05, 0.0, 0.28, OptionType::Call);
        double recovered = implied_vol(bs_price(in), in);
        CHECK_NEAR(recovered, 0.28, 1e-8);
    }
    // A price below the no-arbitrage floor has no implied vol: NaN.
    BsInputs in = base(100, 100, 1.0, 0.05, 0.0, 0.20, OptionType::Call);
    CHECK(std::isnan(implied_vol(0.001, in)));
    CHECK(std::isnan(implied_vol(1000.0, in)));
}

int main() {
    test_textbook_benchmark();
    test_put_call_parity();
    test_price_limits();
    test_greeks_finite_difference();
    test_greek_signs_and_ranges();
    test_monte_carlo_convergence();
    test_antithetic_reduces_error();
    test_implied_vol_round_trip();
    std::printf("%d checks, %d failures\n", checks, failures);
    return failures == 0 ? 0 : 1;
}
