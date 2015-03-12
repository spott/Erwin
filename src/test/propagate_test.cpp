#include <parameters/hamiltonian.hpp>
#include <parameters/laser.hpp>
#include <petsc_cpp/Petsc.hpp>

int main( int argc, const char** argv )
{
    using namespace petsc;
    using namespace Erwin;

    PetscContext pc( argc, argv );

    auto laser = make_LaserParameters( argc, argv );
    cout << laser.print();
    auto hamiltonian = make_HamiltonianParameters( argc, argv );

    auto H0 = hamiltonian.read_field_free();
    auto D = hamiltonian.read_dipole();

    auto E = laser.efield();
    Matrix H( D );

    TimeStepper ts( [&H0, &D, &E]( Vector&&, Matrix&& A, Matrix&& B,
                                   TimeStepper&&, double t ) {
                        A = D;
                        A *= E( t );
                        A += H0;
                        A *= std::complex<double>( 0, -1 );
                        B.shallow_copy( A );
                    },
                    H, {0, 10, .01}, TSCN );

    auto psi = H.get_right_vector();
    psi.set_value( 0, 1. );
    ts.solve( psi );
}
