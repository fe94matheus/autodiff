#include <mpfr.h>
#include <iostream>
#include <string>

namespace flib
{
    class ArbitraryPrecision
    {
    private:
        mpfr_t value;

    public:
        ArbitraryPrecision(mpfr_prec_t prec = 256)
        {
            mpfr_init2(value, prec);
            mpfr_set_d(value, 0.0, MPFR_RNDZ);
        }

        ArbitraryPrecision(const std::string &str, mpfr_prec_t prec = 53, mpfr_rnd_t rnd = MPFR_RNDZ)
        {
            mpfr_init2(value, prec);
            if (mpfr_set_str(value, str.c_str(), 10, rnd) != 0)
            {
                throw std::invalid_argument("Invalid string for MPFR number.");
            }
        }

        ArbitraryPrecision(ArbitraryPrecision const &x, mpfr_prec_t prec = 256)
        {
            mpfr_init2(value, prec);
            mpfr_set(value, x.value, MPFR_RNDZ);
        }

        ~ArbitraryPrecision()
        {
            mpfr_clear(value);
        }

        ArbitraryPrecision(long int i, mpfr_prec_t prec = 53, mpfr_rnd_t rnd = MPFR_RNDZ)
        {
            mpfr_init2(value, prec);
            mpfr_set_si(value, i, rnd);
        }

        ArbitraryPrecision(const ArbitraryPrecision &other)
        {
            mpfr_init2(value, mpfr_get_prec(other.value));
            mpfr_set(value, other.value, MPFR_RNDZ);
        }

        ArbitraryPrecision &operator=(const ArbitraryPrecision &other)
        {
            if (this != &other)
            {
                mpfr_set_prec(value, mpfr_get_prec(other.value));
                mpfr_set(value, other.value, MPFR_RNDZ);
            }
            return *this;
        }

        operator mpfr_ptr()
        {
            return value;
        }

        operator mpfr_srcptr() const
        {
            return value;
        }
    };
} // namespace flib
