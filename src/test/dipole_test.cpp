
#include <petsc_cpp/Petsc.hpp>
//#include <time_independent/time_independent.hpp>
#include <parameters/hamiltonian.hpp>
#include <time_independent/BasisLoader.hpp>
#include <time_independent/dipole_matrix.hpp>
#include <utilities/io.hpp>


int main( int argc, const char** argv )
{
    using namespace petsc;
    using namespace std;
    using namespace Erwin;

    PetscContext pc( argc, argv );

    auto parameters = make_HamiltonianParameters( argc, argv );
    if ( !pc.rank() ) cout << parameters.print();
    if ( !pc.rank() ) parameters.write();

    if ( !parameters.basis ) {
        if ( !pc.rank() ) cerr << "Need a basis!";
        return -1;
    }
    auto old_prototype = parameters.basis->read_prototype();
    auto new_prototype =
        shrink_prototype( old_prototype, parameters.max_basis() );
    parameters.write_prototype( new_prototype );
    if ( !pc.rank() )
        for ( auto& a : new_prototype ) cout << a << endl;

    auto H = make_field_free( new_prototype );
    parameters.write_field_free( H );
    H.print();

    if ( !parameters.basis->ecs_percent ) {
        auto D = make_dipole_matrix<BasisID, double>( *( parameters.basis ),
                                                      new_prototype );

        D.print();

        parameters.write_dipole( D );
    } else {
        auto D = make_dipole_matrix<BasisID, complex<double>>(
            *( parameters.basis ), new_prototype );

        D.print();

        parameters.write_dipole( D );
    }
}
