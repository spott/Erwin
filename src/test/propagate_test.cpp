#define SLEPC

#include <parameters/hamiltonian.hpp>
#include <parameters/laser.hpp>
#include <parameters/propagate.hpp>
#include <parameters/absorber.hpp>
#include <petsc_cpp/Petsc.hpp>
#include <time_dependent/propagator.hpp>
#include <parameters/dipole.hpp>
#include <parameters/eigenstates.hpp>

int main( int argc, const char** argv )
{
    using namespace petsc;
    using namespace Erwin;

    PetscContext pc( argc, argv );

    auto laser = make_LaserParameters( argc, argv );
    auto hamiltonian = make_HamiltonianParameters( argc, argv );
    auto absorber = make_AbsorberParameters( argc, argv );
    auto propagation = make_PropagationParameters( argc, argv );
    auto dipole = make_DipoleParameters( argc, argv );
    cout << laser.print();
    cout << hamiltonian.print();
    cout << absorber.print();
    cout << propagation.print();
    cout << dipole.print();

    laser.write();
    hamiltonian.write();
    absorber.write();
    propagation.write();
    dipole.write();

    auto prototype = hamiltonian.read_prototype();
    vector<double> energy_values;
    transform( prototype.begin(), prototype.end(),
               back_inserter( energy_values ),
               []( auto a ) { return a.e.real(); } );

    auto H0 = hamiltonian.read_field_free();
    auto D = hamiltonian.read_dipole();

    Propagator p( D, H0, propagation.time(), laser.efield() );

    p.register_observable( laser.get_observer() );
    p.register_observable(
        absorber.get_observer( prototype, p.get_operator_vector() ) );
    p.register_observable( dipole.get_observer( D, prototype ) );
    p.register_observable( std::unique_ptr<Observable>(
        new EigenstateObserver( p.H, 10, "./", 10, prototype.front() ) ) );

    cout << p.observables_names() << endl;

    auto psi = p.get_operator_vector();
    psi.set_value( 0, 1. );
    psi.assemble();
    p.run( psi );
    // H0.conjugate();
    // H0.print();

    // auto E = []( double ) { return 0; }; // laser.efield();
    // Matrix H = D.get_empty_matrix();
    // H += H0;
    // H.assemble();

    // TimeStepper ts(
    //     [&H0, &D, E]( Vector&, Matrix& A, Matrix& B, TimeStepper&, double t )
    //     {
    //         A = D;
    //         A *= E( t );
    //         A += H0;
    //         A *= std::complex<double>( 0, -1 );
    //         B.shallow_copy( A );
    //     },
    //     H, {0, 500, .01}, TSCN );

    // ts.set_monitor(
    //     [&energy_values, E]( TimeStepper&, int step, double t, Vector& U ) {
    //         cout << "t: " << t << " e_field: " << E( t ) << " step: " << step
    //              << " norm-1: " << U.norm() - 1 << endl;
    //         abs_square( U ).draw(
    //             []( double d ) { return d == 0.0 ? -20 : std::log10( d ); }
    //             );
    //     } );
    // auto psi = H.get_right_vector();

    // psi.set_all( 1. );
    // psi.assemble();
    // ts.solve( psi );
}
