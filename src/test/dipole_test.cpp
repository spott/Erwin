
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

    auto old_prototype = parameters.basis.read_prototype();
    auto new_prototype =
        shrink_prototype( old_prototype, parameters.max_basis() );
    parameters.write_prototype( new_prototype );

    Matrix H = make_field_free( new_prototype );
    Matrix D = make_dipole_matrix(
        parameters.basis, new_prototype);
}
