#pragma once

#define SLEPC
#include <petsc_cpp/Petsc.hpp>

#include <utilities/types.hpp>
#include <utilities/io.hpp>

namespace Erwin
{

using namespace std;
using namespace petsc;


template <typename>
struct HamiltonianTraits;


template <typename HamiltonianType>
struct Hamiltonian {
    using Scalar = typename HamiltonianTraits<HamiltonianType>::Scalar;
    using QuantumNumbers =
        typename HamiltonianTraits<HamiltonianType>::QuantumNumbers;

    Hamiltonian( vector<Scalar> grid_ ) : H( grid_.size() ), grid( grid_ ) {}

    QuantumNumbers basis_set_inserter( complex<double> ev )
    {
        return static_cast<HamiltonianType*>( this )
            ->basis_set_insert( ev );
    }

    Matrix H;
    vector<Scalar> grid;
};


template <typename Scalar>
struct SphericalHamiltonian final
    : Hamiltonian<SphericalHamiltonian<Scalar>> {
    using ThisType = SphericalHamiltonian<Scalar>;
    using QuantumNumbers =
        typename HamiltonianTraits<ThisType>::QuantumNumbers;

    SphericalHamiltonian( vector<Scalar> grid,
                          function<Scalar( Scalar )> potential,
                          int quantum_number_l )
        : Hamiltonian<ThisType>( grid ), potential_( potential ),
          ll( quantum_number_l )

    {
        // Make the matrix (for second order accurate):
        auto nonzeros = []( int i, int j ) {
            return i == j || i == j + 1 || i == j - 1;
        };

        this->H.reserve( nonzeros );
        populate_matrix( this->H, nonzeros,
                         [this]( int i, int j ) -> Scalar {
            if ( i != j )
                return this->second_derivative_2( i, j );
            else if ( i == j )
                return this->second_derivative_2( i, j ) +
                       this->potential_( this->grid[i] ) +
                       this->centrifugal_potential( this->grid[i] );
            else
                throw domain_error(
                    "attempting to put a value where one doesn't belong" );
        } );

        this->H.assemble();
    }

    // 2nd order accurate 2nd derivative
    Scalar second_derivative_2( int i, int j )
    {
        auto& grid = Hamiltonian<ThisType>::grid;
        if ( i == j ) {
            auto a = grid[i + 1];
            auto b = grid[i];
            auto c = i > 0 ? grid[i - 1] : 0;
            auto dr2 = ( a - b ) * ( b - c );
            return 1. / ( dr2 );
        } else if ( i > j ) {
            auto a = grid[i];

            auto b = grid[i - 1];
            auto c = i > 1 ? grid[i - 2] : 0;
            auto dr2 = ( a - b ) * ( b - c );
            return -1. / ( 2. * dr2 );
        } else if ( i < j ) {
            auto a = grid[i + 2];
            auto b = grid[i + 1];
            auto c = grid[i];
            auto dr2 = ( a - b ) * ( b - c );
            return -1. / ( 2. * dr2 );
        } else
            throw domain_error( "derivative: attempting to put a "
                                "value where one doesn't belong" );
    }

    Scalar centrifugal_potential( Scalar r )
    {
        return ll * ( ll + 1. ) / ( 2. * r * r );
    }

    int& l() { return ll; }
    // reminder: expensive
    int& l( int l )
    {
        if ( l == ll ) return ll;

        ll = l;
        n = l + 1;
        auto nonzeros = []( int i, int j ) {
            return i == j || i == j + 1 || i == j - 1;
        };
        populate_matrix( this->H, nonzeros,
                         [this]( int i, int j ) -> Scalar {
            if ( i != j )
                return this->second_derivative_2( i, j );
            else if ( i == j )
                return this->second_derivative_2( i, j ) +
                       this->potential_( this->grid[i] ) +
                       this->centrifugal_potential( this->grid[i] );
            else
                throw domain_error(
                    "attempting to put a value where one doesn't belong" );
        } );
        return ll;
    }

    // for inserting into vector
    QuantumNumbers basis_set_insert( complex<double> ev )
    {
        return QuantumNumbers{n++, ll, 0, ev};
    }

  private:
    int n;
    function<Scalar( Scalar )> potential_;
    int ll;
};

template <typename Scalar_>
struct HamiltonianTraits<SphericalHamiltonian<Scalar_>> {
    using Scalar = Scalar_;
    using QuantumNumbers = BasisID;
};

template <typename HamiltonianType>
struct Basis {
    using QuantumNumbers =
        typename HamiltonianTraits<HamiltonianType>::QuantumNumbers;
    using Scalar = typename HamiltonianTraits<HamiltonianType>::Scalar;

    Basis( Hamiltonian<HamiltonianType>& H_ )
        : nstates( 100 ), H( H_ ),
          e( H.H,
             nstates,
             EigenvalueSolver::Which::smallest_real,
             EigenvalueSolver::Type::hermitian )
    {
    }

    Basis( Hamiltonian<HamiltonianType>& H_, int num_states )
        : nstates( num_states ), H( H_ ),
          e( H.H,
             num_states,
             EigenvalueSolver::Which::smallest_real,
             EigenvalueSolver::Type::hermitian )

    {
        e.dimensions( num_states, max( num_states, 600 ) );
        e.balance( EPS_BALANCE_TWOSIDE, 10 );
        //e.shift_invert( -1 );
        e.tolerances( 1e-16, 400 );
    }


    void find() { e.solve(); }
    
    void save_basis( string filename )
    {
        // clear file first:
        if ( io::file_exists( filename ) ) io::empty_file( filename );
        save_basis_overload( filename, Scalar() );
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
    void save_basis_overload( string filename, complex<double> )
    {
        e.save_basis( filename );
    }
    void save_basis_overload( string filename, double )
    {
        e.save_basis( filename, []( Vector& v ) {
            map( v, []( auto a, auto i ) { return a.real(); } );
        } );
    }
};

template <typename T>
Basis<T> make_basis( Hamiltonian<T>& H, int num_states )
{
    return Basis<T>( H, num_states );
}
}
