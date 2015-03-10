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

template <typename HamiltonianType,
          bool hermitian_ = HamiltonianTraits<HamiltonianType>::hermitian>
struct Hamiltonian;

template <typename HamiltonianType>
struct Hamiltonian<HamiltonianType, true> {
    using Scalar = typename HamiltonianTraits<HamiltonianType>::Scalar;
    using QuantumNumbers =
        typename HamiltonianTraits<HamiltonianType>::QuantumNumbers;

    // matrix uses grid_.size() - 1.  grid should probably be split out into its
    // own class eventually.
    Hamiltonian( vector<Scalar>& grid_ ) : H( grid_.size() - 1 ), grid( grid_ )
    {
    }

    QuantumNumbers basis_set_inserter( complex<double> ev )
    {
        return static_cast<HamiltonianType*>( this )->basis_set_insert( ev );
    }

    void assemble() { H.assemble(); }

    Matrix H;
    vector<Scalar> grid;
    static constexpr bool hermitian =
        HamiltonianTraits<HamiltonianType>::hermitian;
};

template <typename HamiltonianType>
struct Hamiltonian<HamiltonianType, false> {
    using Scalar = typename HamiltonianTraits<HamiltonianType>::Scalar;
    using QuantumNumbers =
        typename HamiltonianTraits<HamiltonianType>::QuantumNumbers;

    // matrix uses grid_.size() - 1.  grid should probably be split out into its
    // own class eventually.
    Hamiltonian( vector<Scalar>& grid_ ) : H( grid_.size() - 1 ), grid( grid_ )
    {
    }

    QuantumNumbers basis_set_inserter( complex<double> ev )
    {
        return static_cast<HamiltonianType*>( this )->basis_set_insert( ev );
    }

    void assemble()
    {
        H.assemble();
        HT = hermitian_transpose( H );
    }

    Matrix H;

    // Hermitian transpose of H.
    Matrix HT;
    vector<Scalar> grid;
    static constexpr bool hermitian =
        HamiltonianTraits<HamiltonianType>::hermitian;
};


template <typename Scalar>
struct SphericalHamiltonian final : Hamiltonian<SphericalHamiltonian<Scalar>> {
    using ThisType = SphericalHamiltonian<Scalar>;
    using QuantumNumbers = typename HamiltonianTraits<ThisType>::QuantumNumbers;

    SphericalHamiltonian( vector<Scalar> grid,
                          function<Scalar( Scalar )> potential,
                          int quantum_number_l = 0 )
        : Hamiltonian<ThisType>( grid ), potential_( potential ),
          ll( quantum_number_l )

    {
        // Make the matrix (for second order accurate):
        auto nonzeros =
            []( int i, int j ) { return i == j || i == j + 1 || i == j - 1; };

        this->H.reserve( nonzeros );
        populate_matrix( this->H, nonzeros, [this]( int i, int j ) -> Scalar {
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

        this->assemble();
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
        auto nonzeros =
            []( int i, int j ) { return i == j || i == j + 1 || i == j - 1; };
        populate_matrix( this->H, nonzeros, [this]( int i, int j ) -> Scalar {
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

        this->assemble();
        return ll;
    }

    // for inserting into vector
    QuantumNumbers basis_set_insert( complex<double> ev )
    {
        return QuantumNumbers{n++, ll, 0, ev};
    }

  private:
    int n{1};
    function<Scalar( Scalar )> potential_;
    int ll;
};

template <typename Scalar_>
struct HamiltonianTraits<SphericalHamiltonian<Scalar_>> {
    using Scalar = Scalar_;
    using QuantumNumbers = BasisID;
    static constexpr bool hermitian = is_same<Scalar_, double>::value;
};
}
