#pragma once

#include <utilities/types.hpp>

namespace Erwin
{
namespace math
{
    // Constants:

    template <typename T>
    constexpr T PI = static_cast<T>(
        3.1415926535897932384626433832795028841971693993751058209 );
    template <typename T>
    constexpr T C = static_cast<T>( 137.03599907444 );
    template <typename T>
    constexpr T ALPHA = static_cast<T>( 1. / 137.03599907444 );

    template <typename T>
    inline constexpr int signum( T x, std::false_type /*is_signed*/ )
    {
        return T( 0 ) < x;
    }

    template <typename T>
    inline constexpr int signum( T x, std::true_type /*is_signed*/ )
    {
        return ( T( 0 ) < x ) - ( x < T( 0 ) );
    }

    template <typename T>
    inline constexpr int signum( T x )
    {
        return signum( x, std::is_signed<T>() );
    }

    std::complex<double> cg_coefficient( const Angular& init,
                                         const Angular& fin );

    std::vector<std::complex<double>> make_ecs_grid( size_t grid_size,
                                                     double rmax,
                                                     double exterior_percent,
                                                     double alpha );

    std::vector<double> make_equally_spaced_grid( size_t grid_size,
                                                  double rmax );
}
}
