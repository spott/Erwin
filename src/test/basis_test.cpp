#include <time_independent/Basis.hpp>
#include <utilities/math.hpp>
#include <parameters/basis.hpp>

int main( int argc, const char** argv )
{
    // using namespace petsc;
    using namespace std;
    using namespace Erwin;

    PetscContext pc( argc, argv );

    auto parameters = make_BasisParameters( argc, argv );
    if ( !pc.rank() ) cout << parameters.print();
    if ( !pc.rank() ) parameters.write();

    if ( parameters.ecs_percent ) {
        auto grid = Erwin::math::make_ecs_grid(
            parameters.points, parameters.rmax, *( parameters.ecs_percent ),
            *( parameters.ecs_alpha ) );

        auto H = make_SphericalHamiltonian(
            grid, []( auto r ) { return -1. / r; }, 0 );
        auto B = make_Basis( H, parameters.nmax );
        vector<BasisID> prototype;

        Vector a = H.H.get_right_vector();
        Vector b = H.H.get_right_vector();

        for ( auto l = 0u; l <= parameters.lmax; ++l ) {
            if ( !pc.rank() ) cout << "l: " << l;
            H.l( l );
            B.nstates = parameters.nmax - l;
            auto gs = B.find();
            if ( !pc.rank() ) cout << " gs: " << gs << endl;
            // B.el.print();
            // B.er.print();
            // for ( auto i = 0u; i < std::max(parameters.nmax - l - 1, 40); ++i
            // ) {must be exhaustive
            //     auto left = B.el.get_eigenvector( i );
            //     auto right = B.er.get_eigenvector( i );
            //     Vector::draw( {{left, right} );
            //     }
            //}

            B.save_basis( parameters.l_filename_left( l ),
                          parameters.l_filename_right( l ) );
            B.add_evalues( prototype );
        }
        if ( !pc.rank() ) {
            parameters.write_grid( grid );
            parameters.write_prototype( prototype );
        }
    } else {
        auto grid = Erwin::math::make_equally_spaced_grid(
            parameters.points, parameters.rmin, parameters.rmax );

        // for (auto&& a: grid)
        //     cout << a << ", ";

        auto H = make_SphericalHamiltonian(
            grid, []( auto r ) { return -1. / r; }, 0 );

        auto B = make_Basis( H, parameters.nmax );
        vector<BasisID> prototype;
        Vector a = H.H.get_right_vector();
        Vector b = H.H.get_right_vector();

        for ( auto l = 0u; l <= parameters.lmax; ++l ) {
            if ( !pc.rank() ) cout << "l: " << l;
            cout.flush();
            H.l( l );
            B.nstates = parameters.nmax - l;
            B.find();
            for ( auto i = 0u; i < parameters.nmax - l - 1; ++i ) {
                auto left = B.e.get_eigenpair( i );
                if ( i == 0 ) {
                    left.evector.draw();
                    if ( !pc.rank() ) cout << " gs: " << left.evalue << endl;
                }
                if ( i + l + 1 == 2 && l == 0 ) a = left.evector;
                if ( i + l + 1 == 2 && l == 1 ) b = left.evector;
            }

            B.save_basis( parameters.l_filename( l ) );
            B.add_evalues( prototype );
        }
        if ( !pc.rank() ) {
            parameters.write_grid( grid );
            parameters.write_prototype( prototype );
        }
    }
}
