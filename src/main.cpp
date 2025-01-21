#include <iostream>
#include <string>
#include <cmath>
#include <numbers>
#include "autodiff.hpp"
#include "newton_function.hpp"
#include "interval.hpp"
#include "ap_number.hpp"

const uint precison = 53;

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

    interval<mpfr_t, precison> g(interval<mpfr_t, precison> x)
    {
        ArbitraryPrecision a("2", precison);
        ArbitraryPrecision b("1e-1024", precison);
        interval<mpfr_t, precison> y(a, a);
        interval<mpfr_t, precison> w(b, b);

        return (x - y) * (x - y) - w;
    }

    interval<mpfr_t, precison> dg(interval<mpfr_t, precison> x)
    {
        ArbitraryPrecision a("2", precison);
        ArbitraryPrecision b("4", precison);
        interval<mpfr_t, precison> y(a, a);
        interval<mpfr_t, precison> w(b, b);

        return y * x - w;
    }


} // namespace flib

int main()
{
    uint number_digits = std::ceil(precison * std::log10(2));

    using interval = flib::interval<mpfr_t, precison>;
    using ap = flib::ArbitraryPrecision;


    ap num1("2.1", precison);
    ap num2("3", precison);

    interval a(num1, num2);

    for (uint i = 0; i < 20; i++)
    {
        std::cout << a << std::endl;
        std::cout << (a.mid() - (flib::g(a.mid()) / flib::dg(a))) << std::endl;

        a = interval::intersection(a, (a.mid() - (flib::g(a.mid()) / flib::dg(a))));
    }
    std::cout << a << std::endl;

    mpfr_free_cache(); // free the cache for constants like pi

    return 0;
}
