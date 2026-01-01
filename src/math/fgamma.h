#pragma once

#include <cmath>
#include <limits>

namespace eri::math {

namespace detail {

// numeric constants similar to numpy float64
inline constexpr double EPS   = std::numeric_limits<double>::epsilon();
inline constexpr double MINN  = std::numeric_limits<double>::min();   // smallest positive normal
inline constexpr double FPMIN = MINN / EPS;

inline double gammln(double xx) {
    static constexpr double cof[] = {
        57.1562356658629235,  -59.5979603554754912,
        14.1360979747417471,  -0.491913816097620199,
        0.339946499848118887e-4,  0.465236289270485756e-4,
        -0.983744753048795646e-4, 0.158088703224912494e-3,
        -0.210264441724104883e-3, 0.217439618115212643e-3,
        -0.164318106536763890e-3, 0.844182239838527433e-4,
        -0.261908384015814087e-4, 0.368991826595316234e-5
    };

    if (xx <= 0.0)
        return std::numeric_limits<double>::quiet_NaN();

    double x = xx;
    double y = xx;

    double tmp = x + 5.24218750000000000;
    tmp = (x + 0.5) * std::log(tmp) - tmp;
    double ser = 0.999999999999997092;

    for (double c : cof) {
        y += 1.0;
        ser += c / y;
    }

    return tmp + std::log(2.5066282746310005 * ser / x);
}

inline double gser(double a, double x) {
    const double gln = gammln(a);
    double ap = a;

    double d = 1.0 / a;
    double s = d;

    while (true) {
        ap += 1.0;
        d *= x / ap;
        s += d;
        if (std::fabs(d) < std::fabs(s) * EPS) {
            return s * std::exp(-x + a * std::log(x) - gln);
        }
    }
}

inline double gcf(double a, double x) {
    const double gln = gammln(a);

    double b = x + 1.0 - a;
    double c = 1.0 / FPMIN;
    double d = 1.0 / b;
    double h = d;

    double i = 1.0;
    while (true) {
        const double an = -i * (i - a);
        b += 2.0;

        d = an * d + b;
        if (std::fabs(d) < FPMIN) d = FPMIN;

        c = b + an / c;
        if (std::fabs(c) < FPMIN) c = FPMIN;

        d = 1.0 / d;
        const double dd = d * c;
        h *= dd;

        if (std::fabs(dd - 1.0) <= EPS) break;
        i += 1.0;
    }

    return std::exp(-x + a * std::log(x) - gln) * h;
}

inline double gammpapprox(double a, double x, int psig) {
    static constexpr double y[] = {
        0.0021695375159141994, 0.011413521097787704, 0.027972308950302116,
        0.051727015600492421,  0.082502225484340941, 0.12007019910960293,
        0.16415283300752470,   0.21442376986779355,  0.27051082840644336,
        0.33199876341447887,   0.39843234186401943,  0.46931971407375483,
        0.54413605556657973,   0.62232745288031077,  0.70331500465597174,
        0.78649910768313447,   0.87126389619061517,  0.95698180152629142
    };

    static constexpr double w[] = {
        0.0055657196642445571, 0.012915947284065419, 0.020181515297735382,
        0.027298621498568734,  0.034213810770299537, 0.040875750923643261,
        0.047235083490265582,  0.053244713977759692, 0.058860144245324798,
        0.064039797355015485,  0.068745323835736408, 0.072941885005653087,
        0.076598410645870640,  0.079687828912071670, 0.082187266704339706,
        0.084078218979661945,  0.085346685739338721, 0.085983275670394821
    };

    const double a1 = a - 1.0;
    if (a1 <= 0.0)
        return std::numeric_limits<double>::quiet_NaN();

    const double lna1 = std::log(a1);
    const double sqrta1 = std::sqrt(a1);
    const double gln = gammln(a);

    double xu;
    if (x > a1) xu = std::max(a1 + 11.5 * sqrta1, x + 6.0 * sqrta1);
    else        xu = std::max(0.0, std::min(a1 - 7.5 * sqrta1, x - 5.0 * sqrta1));

    double s = 0.0;
    for (int j = 0; j < (int)(sizeof(y)/sizeof(y[0])); ++j) {
        const double t = x + (xu - x) * y[j];
        s += w[j] * std::exp(-(t - a1) + a1 * (std::log(t) - lna1));
    }

    const double ans = s * (xu - x) * std::exp(a1 * (lna1 - 1.0) - gln);

    if (psig != 0) { // return P(a,x)
        return (ans > 0.0) ? (1.0 - ans) : (-ans);
    } else {         // return Q(a,x)
        return (ans >= 0.0) ? ans : (1.0 + ans);
    }
}

inline double gammp(double a, double x) {
    constexpr double ASWITCH = 100.0;

    if (x < 0.0 || a <= 0.0)
        return std::numeric_limits<double>::quiet_NaN();
    if (x == 0.0) return 0.0;

    if (a >= ASWITCH) return gammpapprox(a, x, 1);
    if (x < a + 1.0)  return gser(a, x);
    return 1.0 - gcf(a, x);
}

inline double gamm_inc(double a, double x) {
    const double p = gammp(a, x);
    const double gln = gammln(a);
    return std::exp(gln) * p; // lower incomplete gamma
}

} // namespace detail

// Boys function F_a(x) using your definition
inline double Fgamma(int a, double x) {
    x = std::fabs(x);
    if (x < 1e-7) x = 1e-7;

    const double val = detail::gamm_inc(a + 0.5, x);
    if (std::isnan(val)) return val;

    return 0.5 * std::pow(x, -a - 0.5) * val;
}

} // namespace eri::math
