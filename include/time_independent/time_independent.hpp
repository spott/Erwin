#pragma once

#define SLEPC
#include <petsc_cpp/Petsc.hpp>

#include <utilities/types.hpp>
#include <utilities/io.hpp>
#include <time_independent/Hamiltonian.hpp>

namespace Erwin
{

using namespace std;
using namespace petsc;

template <typename HamiltonianType,
          bool hermitian = HamiltonianTraits<HamiltonianType>::hermitian>
struct Basis;

template <typename HamiltonianType>
struct Basis<HamiltonianType, true> {
    using QuantumNumbers =
        typename HamiltonianTraits<HamiltonianType>::QuantumNumbers;
    using Scalar = typename HamiltonianTraits<HamiltonianType>::Scalar;

    Basis( Hamiltonian<HamiltonianType>& H_, int num_states )
        : nstates( num_states ), H( H_ ),
          e( H.H,
             num_states,
             EigenvalueSolver::Which::smallest_real,
             EigenvalueSolver::Type::hermitian )
    {
        e.dimensions( num_states, max( num_states, 600 ) );
        e.balance( EPS_BALANCE_TWOSIDE, 10 );
        e.inner_product_space( inner_product_space_diag( H.grid ) );
    }

    Basis( Hamiltonian<HamiltonianType>& H_ ) : Basis( H_, 100 ) {}

    Vector inner_product_space_diag( vector<Scalar> grid )
    {
        Vector v = H.H.get_right_vector();
        populate_vector( v, [&grid]( int i ) {
            return i == 0 ? grid[0] : grid[i] - grid[i - 1];
        } );
        return v;
    }

    complex<double> find()
    {
        // If not hermitian, we need to solve the system twice...
        // one for left and one for right
        // not implemented at the moment
        assert( H.hermitian );
        // confirm that we have the right operator
        e.op( H.H );

        // find ground state...
        e.dimensions( 1 );
        e.tolerances( 1e-2, 400 );
        e.shift_invert( -10 );
        e.solve();
        auto ep = e.get_eigenpair( 0 );

        // shift/invert to evalue, and set initial vector to speed up
        // computation.
        e.dimensions( nstates, max( nstates, 600 ) );
        e.tolerances( 1e-16, 400 );
        e.shift_invert( ep.evalue );
        e.set_initial_vector( ep.evector );

        e.solve();
        return ep.evalue;
    }

    void save_basis( string filename )
    {
        // clear file first:
        if ( io::file_exists( filename ) ) io::empty_file( filename );
        e.save_basis( filename, []( Vector& v ) {
            petsc::map( v, []( auto a, auto i ) { return a.real(); } );
        } );
    }

    void save_grid( string filename )
    {
        if ( !H.H.rank() ) io::export_vector_binary( filename, H.grid );
    }

    vector<QuantumNumbers>& add_evalues( vector<QuantumNumbers>& v )
    {
        for ( int i = 0; i < min( nstates, e.num_converged() ); i++ ) {
            v.push_back( H.basis_set_inserter( e.get_eigenvalue( i ) ) );
        }

        return v;
    }

    int nstates;
    Hamiltonian<HamiltonianType>& H;
    EigenvalueSolver e;

  private:
};

template <typename HamiltonianType>
struct Basis<HamiltonianType, false> {
    using QuantumNumbers =
        typename HamiltonianTraits<HamiltonianType>::QuantumNumbers;
    using Scalar = typename HamiltonianTraits<HamiltonianType>::Scalar;


    Basis( Hamiltonian<HamiltonianType>& H_, int num_states )
        : nstates( num_states ), H( H_ ),
          er( H.H,
              nstates,
              EigenvalueSolver::Which::smallest_real,
              EigenvalueSolver::Type::nonhermitian ),
          el( H.HT,
              nstates,
              EigenvalueSolver::Which::smallest_real,
              EigenvalueSolver::Type::nonhermitian )
    {
        er.dimensions( num_states, max( num_states, 600 ) );
        er.balance( EPS_BALANCE_TWOSIDE, 10 );
        er.inner_product_space( inner_product_space_diag( H.grid ) );
        el.dimensions( num_states, max( num_states, 600 ) );
        el.balance( EPS_BALANCE_TWOSIDE, 10 );
        el.inner_product_space(
            move( inner_product_space_diag( H.grid ).conjugate() ) );
    }

    Basis( Hamiltonian<HamiltonianType>& H_ ) : Basis( H_, 100 ) {}

    Vector inner_product_space_diag( vector<Scalar> grid )
    {
        Vector v = H.H.get_right_vector();
        populate_vector( v, [&grid]( int i ) {
            return i == 0 ? grid[0] : grid[i] - grid[i - 1];
        } );
        return v;
    }

    complex<double> find()
    {
        // If not hermitian, we need to solve the system twice...
        // one for left and one for right
        // not implemented at the moment
        // confirm that we have the right operator
        er.op( H.H );

        // find ground state...
        er.dimensions( 1 );
        er.tolerances( 1e-2, 400 );
        er.shift_invert( -10 );
        er.solve();
        auto ep = er.get_eigenpair( 0 );

        // shift/invert to evalue, and set initial vector to speed up
        // computation.
        er.dimensions( nstates, max( nstates, 600 ) );
        er.tolerances( 1e-16, 400 );
        er.shift_invert( ep.evalue );
        er.set_initial_vector( ep.evector );
        er.solve();

        // next version...
        el.op( H.HT );
        el.dimensions( nstates, max( nstates, 600 ) );
        el.tolerances( 1e-16, 400 );
        el.shift_invert( ep.evalue );
        el.set_initial_vector( ep.evector );
        el.solve();


        return ep.evalue;
    }

    void save_basis( string lfilename, string rfilename )
    {
        // clear file first:
        if ( io::file_exists( lfilename ) ) io::empty_file( lfilename );
        if ( io::file_exists( rfilename ) ) io::empty_file( rfilename );
        el.save_basis( lfilename );
        er.save_basis( rfilename );
    }

    void save_grid( string filename )
    {
        if ( !H.H.rank() ) io::export_vector_binary( filename, H.grid );
    }

    vector<QuantumNumbers>& add_evalues( vector<QuantumNumbers>& v )
    {
        for ( int i = 0; i < min( nstates, er.num_converged() ); i++ ) {
            v.push_back( H.basis_set_inserter( er.get_eigenvalue( i ) ) );
        }

        return v;
    }

    int nstates;
    Hamiltonian<HamiltonianType>& H;
    EigenvalueSolver er;
    EigenvalueSolver el;

  private:
};


template <typename T>
Basis<T> make_Basis( Hamiltonian<T>& H, int num_states )
{
    return Basis<T>( H, num_states );
}

template <typename Scalar, typename Funct>
SphericalHamiltonian<Scalar> make_SphericalHamiltonian( vector<Scalar>& grid,
                                                        Funct potential,
                                                        int l_quantum_number )
{
    return SphericalHamiltonian<Scalar>( grid, potential, l_quantum_number );
}
}
