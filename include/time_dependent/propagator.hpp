#pragma once

#include <petsc_cpp/Petsc.hpp>
#include <parameters/propagate.hpp>
#include <time_dependent/observables.hpp>

namespace Erwin
{

struct Propagator {
    template <typename E>
    Propagator( petsc::Matrix& Dipole_,
                petsc::Vector& H0_,
                petsc::TimeStepper::times tt,
                E& ef )
        : Dipole( Dipole_ ), H( Dipole.get_empty_matrix() ), H0( H0_ ),
          observables(), ts( [this, &ef]( petsc::Vector& U,
                                          petsc::Matrix& A,
                                          petsc::Matrix& B,
                                          petsc::TimeStepper& T,
                                          double t ) {
                                 A = this->Dipole;
                                 A *= ef( t );
                                 A += this->H0;
                                 A *= std::complex<double>( 0, -1 );
                                 B.shallow_copy(A);

                                 for ( auto& o : this->observables )
                                     o->modify( A, U, U, T );
                             },
                             H,
                             tt,
                             TSTHETA )
    {
        H.assemble();
        H = this->Dipole;
        H *= ef( tt.ti );
        H += this->H0;
        H *= std::complex<double>( 0, -1 );
        ts.set_monitor( [ over = this, tt ]( petsc::TimeStepper & T, int step,
                                             double t,
                                             const petsc::Vector& U ) {
            static petsc::Draw draw = [&]() {
                petsc::Draw dr( 0, U.size(), -100, 0, 1, U.comm() );
                dr.set_function( []( PetscScalar d, unsigned ) {
                    auto x = std::pow( std::abs( d ), 2 );
                    return x == 0.0 ? -30 : std::log10( x );
                } );
                dr.set_title( "wavefunction propability" );
                return dr;
            }();
            auto norm = U.norm() - 1;
            if ( !draw.rank() )
                std::printf( "t: %8.3f step: %8i norm-1: %8.3e ", t, step,
                             norm );
            for ( auto& o : over->observables ) ( *o )( over->H, U, T );
            for ( auto& o : over->observables )
                if ( !draw.rank() )
                    std::printf( "%s: %s ", o->name().c_str(),
                                 o->last().c_str() );
            if ( !draw.rank() ) cout << endl;
            draw.draw_vector( U );
        } );
    }

    void register_observable( std::unique_ptr<Observable>&& O )
    {
        if ( O != nullptr ) observables.emplace_back( std::move( O ) );
    }

    std::string observables_names()
    {
        using namespace std;
        stringstream ss;
        cout << observables.size() << std::endl;
        for ( auto& o : observables ) ss << o->name() << ", ";
        cout << observables.size() << std::endl;
        return ss.str();
    }

    void run( petsc::Vector psi ) { ts.solve( psi ); }

    petsc::Vector get_operator_vector() { return H.get_right_vector(); }

    petsc::Matrix& Dipole;
    petsc::Matrix H;
    petsc::Vector& H0;
    std::vector<std::unique_ptr<Observable>> observables;
    petsc::TimeStepper ts;
};
}
