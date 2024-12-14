#pragma once
#include "autodiff.hpp"
#include "elementary_functions.hpp"

//----------------------------------------------------------------------------------------
// método de Newton com heranca
//----------------------------------------------------------------------------------------

namespace flib
{
    class newton_function
    {
    public:
        virtual ~newton_function() {}
        virtual fdh operator()(double x) = 0;
    };

    class my_f : public newton_function
    {
    public:
        fdh operator()(double x) override
        {
            fdh fdhx{x, 1, 0};
            fdh r = my_function(fdhx);
            return r;
        }
    };
} // namespace flib
