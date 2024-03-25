#pragma once 
#include <cassert>
#include <cmath>
#include <complex>
#include <cstdint>
#include <random> 
#include <vector>

template<class T>
T fastPower(T base, uint32_t exp)
{
    T result = 1.0;
    while (exp > 0)
    {
        if (exp & 1)
            result *= base;
        exp >>= 1;
        base *= base;
    }
    return result;
}


/**
 * @brief Sample by truncacted geometric distribution
 *
 */
class Geom
{
private:
    using cd = std::complex<double>;
    static constexpr double PI = 3.14159265358979323846;
    std::mt19937 gen;
    std::discrete_distribution<uint32_t> d;

protected:
    static void fft(std::vector<cd>& a, bool invert)
    {
        int n = a.size();

        for (int i = 1, j = 0; i < n; i++)
        {
            int bit = n >> 1;
            for (; j & bit; bit >>= 1)
                j ^= bit;
            j ^= bit;

            if (i < j)
                std::swap(a[i], a[j]);
        }

        for (int len = 2; len <= n; len <<= 1)
        {
            double ang = 2 * PI / len * (invert ? -1 : 1);
            cd wlen(cos(ang), sin(ang));
            for (int i = 0; i < n; i += len)
            {
                cd w(1);
                for (int j = 0; j < len / 2; j++)
                {
                    cd u = a[i + j], v = a[i + j + len / 2] * w;
                    a[i + j] = u + v;
                    a[i + j + len / 2] = u - v;
                    w *= wlen;
                }
            }
        }

        if (invert)
        {
            for (cd& x : a)
                x /= n;
        }
    }
    static std::vector<long double> multiply(std::vector<long double> const& a, std::vector<long double> const& b)
    {
        std::vector<long double> res(a.size() + b.size() - 1, 0);
        for (size_t i = 0; i < a.size(); i++)
            for (size_t j = 0; j < b.size(); j++)
                res[i + j] += a[i] * b[j];
        return res;
        // std::vector<cd> fa(a.begin(), a.end()), fb(b.begin(), b.end());
        // int n = 1;
        // while (n < a.size() + b.size())
        //     n <<= 1;
        // fa.resize(n);
        // fb.resize(n);

        // fft(fa, false);
        // fft(fb, false);
        // for (int i = 0; i < n; i++)
        //     fa[i] *= fb[i];
        // fft(fa, true);

        // std::vector<double> result(n);
        // for (int i = 0; i < n; i++)
        //     result[i] = fa[i].real();
        // result.resize(a.size() + b.size() - 1);
        // return result;
    }

public:
    static std::vector<long double> computeWeight(uint32_t Z_, double epsilon)
    {
        const double alpha = std::exp(-epsilon);
        const uint32_t Z = Z_ / 2;
        std::vector<long double> weights(Z_ + 1, 0.0);
        weights[Z] = (1.0 - alpha) / (1.0 + alpha - 2.0 * fastPower<double>(alpha, Z + 1));

        double multiplier = alpha;

        for (uint32_t i = 1; i <= Z; ++i)
        {
            weights[Z + i] = weights[Z - i] = weights[Z + i - 1] * multiplier;
            multiplier *= alpha;
        }
        return weights;
    }

    static std::vector<long double> geomConv(const uint32_t binCapacity,
        const double epsilon, const uint32_t binNum)
    {
        // compute convolutional distribution
        const std::vector<long double> geom = computeWeight(binCapacity, epsilon);
        std::vector<long double> currentDist(geom);

        for (uint32_t b = 1; b < binNum; ++b)
        {
            currentDist = std::move(multiply(currentDist, geom));
        }
        assert(currentDist.size() == binCapacity * binNum + 1);
        return currentDist;
    }

    Geom(uint32_t Z_, double epsilon_, uint32_t seed_) : gen(seed_)
    {
        assert(Z_ % 2 == 0);
        std::vector<long double> weights = computeWeight(Z_, epsilon_);
        d = std::discrete_distribution<uint32_t>(weights.begin(), weights.end());
    }

    Geom() {}

    uint32_t sample()
    {
        return d(gen);
    }
};


/**
 * @brief Sample by truncacted geometric distribution
 *
 */
class Laplace
{
public:
    Laplace(uint32_t seed_) : gen(seed_), unif(0.0, 1.0) {}

    double sample(double b)
    {
        double cdf = unif(gen);
        return -b * (cdf > 0.5 ? 1.0 : -1.0) * log(1.0 - 2.0 * fabs(cdf - 0.5));
    }
private:
    std::mt19937 gen;
    std::uniform_real_distribution<double> unif;
};
