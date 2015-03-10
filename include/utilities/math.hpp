#pragma once

#include <utilities/types.hpp>

namespace Erwin
{
namespace math
{
    // Constants:

    template <typename T>
    const T PI = T( std::atan( 1.0 ) * 4.0 );
    template <typename T>
    const T C = T( 137.03599907444 );
    template <typename T>
    const T ALPHA = T( 1. / 137.03599907444 );

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

    std::vector<std::complex<double>> make_ecs_grid( int grid_size,
                                                     double rmax,
                                                     double exterior_percent,
                                                     double alpha );

    std::vector<double> make_equally_spaced_grid( int grid_size, double rmax );
}
}
