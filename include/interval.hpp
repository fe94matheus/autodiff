#pragma once
#include "mpfr.h"
#include <iostream>

namespace flib
{
    template <class T, size_t Prec>
    class interval
    {
    private:
        mpfr_t l;
        mpfr_t u;

    public:
        interval(T a, T b)
        {
            if (mpfr_cmp(a, b) > 0)
            {
                throw std::invalid_argument("Invalid interval: lower bound is greater than upper bound");
            }
            mpfr_init2(l, Prec);
            mpfr_init2(u, Prec);

            mpfr_set(l, a, MPFR_RNDD);
            mpfr_set(u, b, MPFR_RNDU);
        }

        ~interval()
        {
            mpfr_clear(l);
            mpfr_clear(u);
        }

        //---------------------------------------
        // Lower bound
        //---------------------------------------

        interval lower()
        {
            mpfr_t r;
            mpfr_init2(r, Prec);

            interval<T, Prec> result(this->l, this->l);

            mpfr_clear(r);

            return result;
        }

        //---------------------------------------
        // Upper bound
        //---------------------------------------

        interval upper()
        {
            mpfr_t r;
            mpfr_init2(r, Prec);

            interval<T, Prec> result(this->u, this->u);

            mpfr_clear(r);

            return result;
        }

        //---------------------------------------
        // Returns the magnitude of the interval
        //---------------------------------------

        interval width()
        {
            mpfr_t r;
            mpfr_init2(r, Prec);

            mpfr_sub(r, u, l, MPFR_RNDN);

            interval<T, Prec> result(r, r);

            mpfr_clear(r);

            return result;
        }

        //---------------------------------------
        // Returns the norm of the interval
        //---------------------------------------

        interval norm()
        {
            mpfr_t r;
            mpfr_t abs_l;
            mpfr_t abs_u;

            mpfr_inits2(Prec, r, abs_l, abs_u, NULL);

            mpfr_abs(abs_l, l, MPFR_RNDD);
            mpfr_abs(abs_u, u, MPFR_RNDU);

            mpfr_max(r, abs_l, abs_u, MPFR_RNDU);
            interval<T, Prec> result(r, r);

            mpfr_clears(r, abs_l, abs_u, NULL);

            return result;
        }

        //---------------------------------------
        // Returns the mid point of the interval
        //---------------------------------------

        interval mid()
        {
            mpfr_t r;

            mpfr_init2(r, Prec);

            mpfr_add(r, this->u, this->l, MPFR_RNDU);
            mpfr_div_ui(r, r, 2, MPFR_RNDN);

            interval<T, Prec> result(r, r);

            mpfr_clear(r);

            return result;
        }

        //---------------------------------------
        // [a , b] + [c , d] = [a+c, b+d]
        //---------------------------------------

        interval operator+(const interval &iv) const
        {
            mpfr_t result_l;
            mpfr_t result_u;

            mpfr_init2(result_l, Prec);
            mpfr_init2(result_u, Prec);

            mpfr_add(result_l, l, iv.l, MPFR_RNDD);
            mpfr_add(result_u, u, iv.u, MPFR_RNDU);

            interval<T, Prec> r(result_l, result_u);

            mpfr_clear(result_l);
            mpfr_clear(result_u);

            return r;
        }

        //---------------------------------------
        // [a , b] + alpha = [a+alpha, b+alpha]
        //---------------------------------------

        interval operator+(const T &a) const
        {
            mpfr_t result_l;
            mpfr_t result_u;

            mpfr_init2(result_l, Prec);
            mpfr_init2(result_u, Prec);

            mpfr_add(result_l, l, a, MPFR_RNDD);
            mpfr_add(result_u, u, a, MPFR_RNDU);

            interval<T, Prec> r(result_l, result_u);

            mpfr_clear(result_l);
            mpfr_clear(result_u);

            return r;
        }

        //-------------------------------------------
        // alpha + [a , b]  = [alpha + a, alpha + b]
        //-------------------------------------------

        friend interval operator+(const T &a, const interval &iv)
        {
            mpfr_t result_l;
            mpfr_t result_u;

            mpfr_init2(result_l, Prec);
            mpfr_init2(result_u, Prec);

            mpfr_add(result_l, a, iv.l, MPFR_RNDD);
            mpfr_add(result_u, a, iv.u, MPFR_RNDU);

            interval<T, Prec> r(result_l, result_u);

            mpfr_clear(result_l);
            mpfr_clear(result_u);

            return r;
        }

        //---------------------------------------
        // [a , b] - alpha = [ a - alpha, b - alpha ]
        //---------------------------------------

        interval operator-(const T &a) const
        {
            mpfr_t result_l;
            mpfr_t result_u;

            mpfr_init2(result_l, Prec);
            mpfr_init2(result_u, Prec);

            mpfr_sub(result_l, l, a, MPFR_RNDD);
            mpfr_sub(result_u, u, a, MPFR_RNDU);

            interval<T, Prec> r(result_l, result_u);

            mpfr_clear(result_l);
            mpfr_clear(result_u);

            return r;
        }

        //--------------------------------------------
        // alpha - [a , b]  = [alpha - a, alpha - b]
        //--------------------------------------------

        friend interval operator-(const T &a, const interval &iv)
        {
            mpfr_t result_l;
            mpfr_t result_u;

            mpfr_init2(result_l, Prec);
            mpfr_init2(result_u, Prec);

            mpfr_add(result_l, a, (-iv).l, MPFR_RNDD);
            mpfr_add(result_u, a, (-iv).u, MPFR_RNDU);

            interval<T, Prec> r(result_l, result_u);

            mpfr_clear(result_l);
            mpfr_clear(result_u);

            return r;
        }

        //---------------------------------------
        // [a , b] - [c , d] = [a-d, b-c]
        //---------------------------------------

        interval operator-(const interval &iv) const
        {
            mpfr_t result_l;
            mpfr_t result_u;

            mpfr_init2(result_l, Prec);
            mpfr_init2(result_u, Prec);

            mpfr_sub(result_l, l, iv.u, MPFR_RNDD);
            mpfr_sub(result_u, u, iv.l, MPFR_RNDU);

            interval<T, Prec> r(result_l, result_u);

            mpfr_clear(result_l);
            mpfr_clear(result_u);

            return r;
        }

        //---------------------------------------
        // - [c , d] = [-d, -c]
        //---------------------------------------

        interval operator-() const
        {
            mpfr_t result_l;
            mpfr_t result_u;

            mpfr_init2(result_l, Prec);
            mpfr_init2(result_u, Prec);

            mpfr_neg(result_l, u, MPFR_RNDD);
            mpfr_neg(result_u, l, MPFR_RNDU);

            interval<T, Prec> r(result_l, result_u);

            mpfr_clear(result_l);
            mpfr_clear(result_u);

            return r;
        }

        //----------------------------------------------------------
        // [a , b] * [c , d] = [min A, max A], A = {ac, ad, bc, bd}
        //----------------------------------------------------------

        interval operator*(const interval &iv) const
        {
            mpfr_t p1, p2, p3, p4;
            mpfr_t result_l;
            mpfr_t result_u;
            mpfr_inits2(Prec, p1, p2, p3, p4, result_l, result_u, NULL);

            mpfr_mul(p1, l, iv.l, MPFR_RNDD);
            mpfr_mul(p2, l, iv.u, MPFR_RNDD);
            mpfr_mul(p3, u, iv.l, MPFR_RNDU);
            mpfr_mul(p4, u, iv.u, MPFR_RNDU);

            mpfr_min(result_l, p1, p2, MPFR_RNDD);
            mpfr_min(result_l, result_l, p3, MPFR_RNDD);
            mpfr_min(result_l, result_l, p4, MPFR_RNDD);

            mpfr_max(result_u, p1, p2, MPFR_RNDU);
            mpfr_max(result_u, result_u, p3, MPFR_RNDU);
            mpfr_max(result_u, result_u, p4, MPFR_RNDU);

            interval<T, Prec> r(result_l, result_u);

            mpfr_clears(p1, p2, p3, p4, result_l, result_u, NULL);

            return r;
        }

        //----------------------------------------------------------
        //  [a , b] * alpha  = [alpha * a , alpha * b]
        //----------------------------------------------------------

        interval operator*(const T &a) const
        {
            mpfr_t result_l;
            mpfr_t result_u;
            mpfr_inits2(Prec, result_l, result_u, NULL);

            mpfr_mul(result_l, l, a, MPFR_RNDD);
            mpfr_mul(result_u, u, a, MPFR_RNDU);

            interval<T, Prec> r(result_l, result_u);

            mpfr_clears(result_l, result_u, NULL);

            return r;
        }

        friend interval operator*(const T &a, const interval &iv)
        {
            mpfr_t result_l;
            mpfr_t result_u;
            mpfr_inits2(Prec, result_l, result_u, NULL);

            mpfr_mul(result_l, iv.l, a, MPFR_RNDD);
            mpfr_mul(result_u, iv.u, a, MPFR_RNDU);

            interval<T, Prec> r(result_l, result_u);

            mpfr_clears(result_l, result_u, NULL);

            return r;
        }

        //----------------------------------------------------------
        //  [a , b] / alpha  = [ a / alpha, b / alpha ]
        //----------------------------------------------------------

        interval operator/(const T &a) const
        {

            if (a == 0)
            {
                throw std::invalid_argument("Division by zero is undefined");
            }

            mpfr_t result_l;
            mpfr_t result_u;
            mpfr_inits2(Prec, result_l, result_u, NULL);

            mpfr_div(result_l, l, a, MPFR_RNDD);
            mpfr_div(result_u, u, a, MPFR_RNDU);

            interval<T, Prec> r(result_l, result_u);

            mpfr_clears(result_l, result_u, NULL);

            return r;
        }

        //----------------------------------------------------------
        //  alpha / [a , b] = [ alpha / b , alpha / a ]
        //----------------------------------------------------------

        friend interval operator/(const T &a, const interval &iv)
        {
            int sign_l = mpfr_sgn(iv.l);
            int sign_u = mpfr_sgn(iv.u);

            if (sign_l * sign_u <= 0)
            {
                throw std::domain_error("Division by an interval containing zero is undefined");
            }

            mpfr_t result_l;
            mpfr_t result_u;
            mpfr_inits2(Prec, result_l, result_u, NULL);

            mpfr_div(result_l, a, iv.u, MPFR_RNDD);
            mpfr_div(result_u, a, iv.l, MPFR_RNDU);

            interval<T, Prec> r(result_l, result_u);

            mpfr_clears(result_l, result_u, NULL);

            return r;
        }

        //----------------------------------------------------------
        // [a , b] / [c, d] =  [a , b] * ( 1 / [c, d] )
        //----------------------------------------------------------

        interval operator/(const interval &iv) const
        {

            int sign_l = mpfr_sgn(iv.l);
            int sign_u = mpfr_sgn(iv.u);

            if (sign_l * sign_u <= 0)
            {
                throw std::domain_error("Division by an interval containing zero is undefined");
            }

            mpfr_t p1, p2, p3, p4;
            mpfr_t result_l;
            mpfr_t result_u;
            mpfr_t recip_l, recip_u;

            mpfr_inits2(Prec, p1, p2, p3, p4, result_l, result_u, recip_l, recip_u, NULL);

            mpfr_ui_div(recip_l, 1, iv.u, MPFR_RNDD);
            mpfr_ui_div(recip_u, 1, iv.l, MPFR_RNDU);

            mpfr_mul(p1, l, recip_l, MPFR_RNDD);
            mpfr_mul(p2, l, recip_u, MPFR_RNDD);
            mpfr_mul(p3, u, recip_l, MPFR_RNDU);
            mpfr_mul(p4, u, recip_u, MPFR_RNDU);

            mpfr_min(result_l, p1, p2, MPFR_RNDD);
            mpfr_min(result_l, result_l, p3, MPFR_RNDD);
            mpfr_min(result_l, result_l, p4, MPFR_RNDD);

            mpfr_max(result_u, p1, p2, MPFR_RNDU);
            mpfr_max(result_u, result_u, p3, MPFR_RNDU);
            mpfr_max(result_u, result_u, p4, MPFR_RNDU);

            interval<T, Prec> r(result_l, result_u);

            mpfr_clears(p1, p2, p3, p4, result_l, result_u, recip_l, recip_u, NULL);

            return r;
        }

        friend std::ostream &operator<<(std::ostream &os, const interval &iv)
        {
            char lb[Prec], ub[Prec];
            mpfr_sprintf(lb, "%.20Re", iv.l);
            mpfr_sprintf(ub, "%.20Re", iv.u);
            os << "[ " << lb << " , " << ub << " ]";
            return os;
        }

        //------------------------------------------------
        //Intersection
        //------------------------------------------------

        static interval intersection(const interval &a, const interval &b)  
        {
            if (mpfr_cmp(b.u, a.l) < 0 || mpfr_cmp(a.u, b.l) < 0)
            {
                throw std::domain_error("Empty intersection");
            }
            
            mpfr_t min, max;
            mpfr_inits2(Prec, min, max, NULL);

            mpfr_max(max, a.l, b.l, MPFR_RNDU);
            mpfr_min(min, a.u, b.u, MPFR_RNDD);

            interval<T, Prec> r(max, min);

            return r;
        }
    };

} // namespace flib
