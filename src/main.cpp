#include <iostream>
#include <string>
#include <cmath>
#include <numbers>
#include "autodiff.hpp"
#include "newton_function.hpp"
#include "interval.hpp"

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

    void teste_interval()
    {
        mpfr_t value1d;
        mpfr_t value1u;
        mpfr_t result1;

        mpfr_t x;

        int p = 128;

        mpfr_init2(value1d, p);
        mpfr_init2(value1u, p);
        mpfr_init2(x, p);
        mpfr_init2(result1, p);

        mpfr_const_pi(x, MPFR_RNDD);

        mpfr_pow_si(value1d, x, 25, MPFR_RNDD);
        mpfr_pow_si(value1u, x, 25, MPFR_RNDU);

        mpfr_printf("value1d: %.10RNa\n", value1d);
        mpfr_printf("value1u: %.10RNa\n", value1u);

        mpfr_sub(result1, value1u, value1d, MPFR_RNDN);

        mpfr_printf("result1: %.10Re\n", result1);

        mpfr_clear(x);
        mpfr_clear(value1d);
        mpfr_clear(value1u);
        mpfr_clear(result1);
    }

} // namespace flib

int main()
{
    mpfr_t x;
    mpfr_t y, z, w;

    const uint precison = 1024;
    uint number_digits = std::ceil(precison * std::log10(2));

    char buffer[number_digits];

    mpfr_init2(x, precison);
    mpfr_init2(y, precison);
    mpfr_init2(z, precison);
    mpfr_init2(w, precison);

    mpfr_set_str(x, "-1",10, MPFR_RNDN);
    mpfr_set_str(y, "2", 10, MPFR_RNDN);
    mpfr_set_str(z, "2", 10, MPFR_RNDN);
    mpfr_set_str(w, "4", 10, MPFR_RNDN);

    using interval = flib::interval<mpfr_t, precison>;

    interval a(x, y);
    interval b(z, w);

    std::cout << a << std::endl;
    std::cout << b << std::endl;

    std::cout << interval::intersection(a, b).width() << std::endl;

    mpfr_clear(x);
    mpfr_clear(y);
    mpfr_clear(z);
    mpfr_clear(w);
    mpfr_free_cache(); // free the cache for constants like pi

    return 0;
}
