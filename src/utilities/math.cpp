#include <utilities/math.hpp>
#include <complex>
#include <vector>
// gsl
#include <gsl/gsl_sf_coupling.h>

namespace Erwin
{
namespace math
{

    std::complex<double> cg_coefficient( const Angular& init,
                                         const Angular& fin )
    {
        auto l_init = static_cast<int>( init.l );
        auto l_fin = static_cast<int>( fin.l );
        std::complex<double> out =
            gsl_sf_coupling_3j( l_init * 2, 2, l_fin * 2, 0, 0, 0 );
        out *= out;
        out *= std::sqrt( ( 2 * l_fin + 1 ) * ( 2 * l_fin + 1 ) );
        return out;
    }

    std::vector<std::complex<double>> make_ecs_grid( size_t grid_size,
                                                     double rmax,
                                                     double exterior_percent,
                                                     double alpha )
    {
        typedef std::complex<double> complex;
        using std::vector;
        std::vector<complex> grid;
        grid.reserve( grid_size + 1 );

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

    std::vector<double> make_equally_spaced_grid( size_t grid_size,
                                                  double rmax )
    {
        using std::vector;
        vector<double> grid;
        grid.reserve( grid_size + 1 );

        auto dr = rmax / grid_size;

        while ( grid.size() <= grid_size ) {
            auto r = dr * ( grid.size() + 1 );
            grid.emplace_back( r );
        }
        return grid;
    }
}
}
