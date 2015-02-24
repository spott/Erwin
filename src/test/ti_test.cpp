#include <time_independent/time_independent.hpp>
#include <utilities/math.hpp>
#include <parameters/basis.hpp>

int main( int argc, const char** argv )
{
    using namespace petsc;
    using namespace std;
    using namespace Erwin;

    PetscContext pc( argc, argv );

    auto parameters = make_BasisParameters( argc, argv );
    if ( !pc.rank() ) cout << parameters.print();
    if ( !pc.rank() ) parameters.write();

    auto grid =
        math::make_ecs_grid( parameters.points, parameters.rmax,
                             parameters.ecs_percent, parameters.ecs_alpha );

    // auto potential = ;
    auto H = make_SphericalHamiltonian(
        grid, []( complex<double> r ) { return -1. / r; }, 0 );
    auto B = make_Basis( H, parameters.nmax );
    vector<BasisID> prototype;

    for ( int l = 0; l <= parameters.lmax; ++l ) {
        cout << "l: " << l;
        H.l( l );
        auto gs = B.find();
        if ( !pc.rank() ) cout << " gs: " << gs << endl;
        B.er.print();
        B.el.print();
        B.save_basis( parameters.l_filename_left( l ),
                      parameters.l_filename_right( l ) );
        B.add_evalues( prototype );
    }
    B.save_grid( parameters.grid_filename() );
    io::export_vector_binary( parameters.prototype_filename(), prototype );
}
