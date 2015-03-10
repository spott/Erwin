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
        std::complex<double> out =
            gsl_sf_coupling_3j( init.l * 2, 2, fin.l * 2, 0, 0, 0 );
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

    std::vector<double> make_equally_spaced_grid( int grid_size, double rmax )
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
