#pragma once

#include <cmath>
#include <cstdint>

namespace ope {

// splitmix64: a small, fast, fully deterministic PRNG. Same seed -> same
// stream on every platform, which is what makes the Monte Carlo results in
// this repo reproducible bit-for-bit.
class SplitMix64 {
public:
    explicit SplitMix64(std::uint64_t seed) : state_(seed) {}

    std::uint64_t next_u64() {
        std::uint64_t z = (state_ += 0x9E3779B97F4A7C15ull);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
        return z ^ (z >> 31);
    }

    // Uniform on [0, 1).
    double next_uniform() {
        return static_cast<double>(next_u64() >> 11) * (1.0 / 9007199254740992.0);
    }

    // Standard normal via Box-Muller; the second variate is cached.
    double next_normal() {
        if (has_spare_) {
            has_spare_ = false;
            return spare_;
        }
        double u1;
        do {
            u1 = next_uniform();
        } while (u1 <= 0.0);
        double u2 = next_uniform();
        double radius = std::sqrt(-2.0 * std::log(u1));
        double angle = 2.0 * kPi * u2;
        spare_ = radius * std::sin(angle);
        has_spare_ = true;
        return radius * std::cos(angle);
    }

private:
    static constexpr double kPi = 3.14159265358979323846;
    std::uint64_t state_;
    double spare_ = 0.0;
    bool has_spare_ = false;
};

}  // namespace ope
