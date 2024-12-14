#pragma once
#include "autodiff.hpp"

//   h(x) = g( f(x) )
//
//  h'(x) = g'( f(x) ) f'(x)
//
// h''(x) = g''( f(x) ) f'(x) f'(x) + g'( f(x) ) f''(x)

namespace flib
{
    inline fdh cos(fdh x)
    {
        return fdh{std::cos(x.f), -std::sin(x.f) * x.d, - std::cos(x.f) * x.d * x.d - std::sin(x.f) * x.h};
    }

    inline fdh sin(fdh x)
    {
        return fdh{std::sin(x.f), std::cos(x.f) * x.d, std::cos(x.f) * x.h - std::sin(x.f) * x.d * x.d};
    }

    inline fdh exp(fdh x)
    {
        double e = std::exp(x.f);
        return fdh{e, e * x.d, e * x.d * x.d  + e * x.h};
    }

    inline fdh log(fdh x)
    {
        return fdh { std::log(x.f), x.d / x.f, (x.h * x.f - x.d * x.d) / ( x.d * x.d ) };
    }

    inline fdh x_pwr_k(fdh x, int k)
    {
        return fdh { std::pow(x.f, k),   k * std::pow(x.f, k-1) * x.d, k * (k - 1) * std::pow(x.f, k - 2) * x.d + k * std::pow(x.f, k - 1) * x.h};
    }

    // template<class T>
    // T my_function(T x)
    // {
    //     return (x_pwr_k(x, 2) + sin(exp(x)))/cos(x) - log(x);
    // }

    template<class T>
    T my_function(T x)
    {
        return sin( x ) / exp(x) ;
    }
} // namespace flib
