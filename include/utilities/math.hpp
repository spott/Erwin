#pragma once

// gsl
#include <gsl/gsl_sf_coupling.h>

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

    struct Angular {
        unsigned int l;
        int m;
    };

    template <typename scalar>
    scalar CGCoefficient( const Angular& init, const Angular& fin )
    {
        scalar out = gsl_sf_coupling_3j( init.l * 2, 2, fin.l * 2, 0, 0, 0 );
        out *= out;
        out *= std::sqrt( ( 2 * init.l + 1 ) * ( 2 * fin.l + 1 ) );
        return out;
    }

    std::vector<std::complex<double>> make_ecs_grid( int grid_size,
                                                     double rmax,
                                                     double exterior_percent,
                                                     double alpha )
    {
        typedef std::complex<double> complex;
        using std::vector;
        std::vector<complex> grid;
        grid.reserve( grid_size );

        auto dr = rmax / grid_size;

        auto boundary = rmax - rmax * exterior_percent;

        while ( grid.size() <= grid_size ) {
            auto r = dr * ( grid.size() + 1 );
            if ( r <= boundary )
                grid.emplace_back( r );
            else
                grid.emplace_back(
                    boundary + ( r - boundary ) * exp( complex( 0, alpha ) ) );
        }

        return grid;
    }

    std::vector<double> make_equally_spaced_grid( int grid_size, double rmax )
    {
        using std::vector;
        vector<double> grid;
        grid.reserve( grid_size );

        auto dr = rmax / grid_size;

        while ( grid.size() <= grid_size ) {
            auto r = dr * ( grid.size() + 1 );
            grid.emplace_back( r );
        }
        return grid;
    }
}
}
