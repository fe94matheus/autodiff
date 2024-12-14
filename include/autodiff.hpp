#pragma once
#include <cassert>
#include <cmath>

//----------------------------------------------------------------------------------------
// diferenciação automática forward.
//----------------------------------------------------------------------------------------

namespace flib
{
    struct fdh
    {
        double f;
        double d;
        double h;

        friend std::ostream &operator<<(std::ostream &os, const fdh &obj)
        {
            os << "fdh(f: " << obj.f << ", d: " << obj.d << ", h: " << obj.h << ")";
            return os;
        }

        //--------------------
        // f(x) * a
        //--------------------

        fdh operator*(double x) const
        {
            return fdh{x * f, x * d, x * h};
        }

        //--------------------
        // a * f(x)
        //--------------------

        friend fdh operator*(double x, fdh u)
        {
            return fdh{x * u.f, x * u.d, x * u.h};
        }

        //--------------------
        // f(x) / a
        //--------------------

        fdh operator/(double x) const
        {
            return fdh{f / x, d / x, h / x};
        }

        //--------------------
        // a / f(x)
        //--------------------

        friend fdh operator/(double x, fdh u)
        {
            return fdh{x / u.f, -(x / (u.f * u.f)) * u.d, ((2 * x * u.d * u.d) - x * u.h * u.d * u.d) / (u.f * u.f * u.f * u.f)};
        }

        //--------------------
        // f(x) + a
        //--------------------

        fdh operator+(double x) const
        {
            return fdh{x + f, d, h};
        }

        //--------------------
        // a + f(x)
        //--------------------

        friend fdh operator+(double x, fdh u)
        {
            return fdh{x + u.f, u.d, u.h};
        }

        //--------------------
        // f(x) - a
        //--------------------

        fdh operator-(double x) const
        {
            return fdh{f - x, d, h};
        }

        //--------------------
        // a - f(x)
        //--------------------

        friend fdh operator-(double x, fdh u)
        {
            return fdh{x - u.f, -u.d, -u.h};
        }

        //--------------------
        // f(x) / g(x)
        //--------------------

        fdh operator/(fdh x) const
        {
            return fdh{
                f / x.f,
                (d * x.f - f * x.d) / (x.f * x.f),
                (h * x.f *  x.f * x.f - f * x.h * x.f * x.f - 2 * (d * x.f - f * x.d) * x.f * x.d) / (x.f * x.f * x.f * x.f)
            };
        }

        //--------------------
        // f(x) * g(x)
        //--------------------

        fdh operator*(fdh x) const
        {
            return fdh{f * x.f, f * x.d + d * x.f, d * x.d + f * x.h + h * x.f + d * x.d};
        }

        //--------------------
        // f(x) + g(x)
        //--------------------

        fdh operator+(fdh x) const
        {
            return fdh{f + x.f, d + x.d, h + x.h};
        }

        //--------------------
        // -f(x) 
        //--------------------

        fdh operator-() const
        {
            return fdh{-f, -d, -h};
        }

        //--------------------
        // f(x) - g(x)
        //--------------------

        fdh operator-(fdh x) const
        {
            return fdh{f - x.f, d - x.d, h - x.h};
        }
    };

} // namespace flib
