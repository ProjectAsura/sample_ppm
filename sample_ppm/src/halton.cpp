﻿//-----------------------------------------------------------------------------
// File : halton.cpp
// Desc : Halto Sequence.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes
//-----------------------------------------------------------------------------
#include <math_util.h>


namespace /* anonymous */ {

// Halton sequence with reverse permutation
int primes[61] = {
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79,
    83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181,
    191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283
};

inline int rev( const int i, const int p )
{
    if (i == 0)
        return i;

    return p - i;
}

} // namespace /*anonymous */


double halton( const int b, int j )
{
    const auto p = primes[b];
    auto h = 0.0;
    auto f = 1.0 / (double)p;
    auto fct = f;
    while ( j > 0 ) 
    {
        h += rev( j % p, p ) * fct;
        j /= p;
        fct *= f;
    }
    return h;
}