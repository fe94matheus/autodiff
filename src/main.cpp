#include <iostream>
#include <string>
#include <cmath>
#include <numbers>
#include "autodiff.hpp"
#include "newton_function.hpp"
#include "interval.hpp"
#include "ap_number.hpp"

const uint precison = 1024;

namespace flib
{
    double newton(newton_function &f, double x)
    {
        for (int i = 0; i < 10; ++i)
        {
            fdh fdhx = f(x);
            x = x - fdhx.f / fdhx.d;
        }

        return x;
    }

    void test()
    {
        // double x = std::numbers::pi * 0.25;
        double x = 6;
        my_f f;

        // double a = newton(f, x);

        fdh r = my_function(fdh{2, 1, 0});
    }

    interval<mpfr_t, precison> f(interval<mpfr_t, precison> x)
    {
        ArbitraryPrecision a("2", precison);
        interval<mpfr_t, precison> y(a, a);

        return interval<mpfr_t, precison>::exp(x) - y;
    }

} // namespace flib

int main()
{
    uint number_digits = std::ceil(precison * std::log10(2));

    using interval = flib::interval<mpfr_t, precison>;
    using ap = flib::ArbitraryPrecision;

    ap num1("0", precison);
    ap num2("1", precison);

    interval a(num1, num2);

    for (uint i = 0; i < 20; i++)
    {
        a = interval::intersection(a, (a.mid() - (flib::f(a.mid()) / interval::exp(a))));
    }
    std::cout << a << std::endl;

    mpfr_free_cache(); // free the cache for constants like pi

    return 0;
}

