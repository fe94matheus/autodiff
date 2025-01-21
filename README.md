# Interval Arithmetic Library

A C++ template library for performing interval arithmetic with arbitrary precision using MPFR. This library provides robust interval arithmetic operations with directed rounding for reliable numerical computations.

## Features

* Arbitrary precision interval arithmetic using MPFR
* Template-based implementation for flexibility
* Comprehensive set of arithmetic operations
* Interval analysis operations (width, norm, midpoint)
* Exception handling for undefined operations
* Directed rounding for guaranteed containment
* Stream output support for intervals

## Dependencies

* MPFR (GNU Multiple Precision Floating-Point Reliable Library)
* C++ Standard Library

## Usage

### Basic Interval Operations

```cpp
#include "interval.hpp"
using namespace flib;

// Create intervals with specified precision
interval<mpfr_t, 256> a(1.0, 2.0);  // Interval [1.0, 2.0] with 256-bit precision
interval<mpfr_t, 256> b(2.0, 3.0);  // Interval [2.0, 3.0]

// Basic arithmetic operations
auto sum = a + b;        // Addition
auto diff = a - b;       // Subtraction
auto prod = a * b;       // Multiplication
auto quot = a / b;       // Division
```

### Analysis Operations

```cpp
// Get interval properties
auto width = a.width();   // Width of interval
auto mid = a.mid();      // Midpoint
auto norm = a.norm();    // Norm (maximum absolute value)
auto lower = a.lower();  // Lower bound
auto upper = a.upper();  // Upper bound
```

### Mathematical Functions

```cpp
// Mathematical operations
auto exp_result = interval<mpfr_t, 256>::exp(a);   // Exponential
auto sqrt_result = interval<mpfr_t, 256>::sqrt(a);  // Square root
```

### Set Operations

```cpp
// Set operations
auto intersect = interval<mpfr_t, 256>::intersection(a, b);  // Intersection of intervals
```

## Error Handling

The library includes robust error handling for undefined operations:
* Division by an interval containing zero
* Square root of negative intervals
* Invalid interval bounds (lower > upper)
* Empty intersection of intervals

## Implementation Details

* Uses directed rounding for reliable results
* Lower bound rounds down (MPFR_RNDD)
* Upper bound rounds up (MPFR_RNDU)
* Template parameter `Prec` determines precision in bits
* RAII design for proper MPFR resource management

## Contributing

When contributing to this project, please:
* Follow the existing code style
* Add tests for new functionality
* Document any changes in behavior
* Update the README for significant changes

## License

This project is licensed under the MIT License - see the LICENSE file for details.