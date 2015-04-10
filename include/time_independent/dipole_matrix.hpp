#pragma once

#include <time_independent/BasisLoader.hpp>
#include <utilities/math.hpp>

#include <petsc_cpp/Petsc.hpp>

namespace Erwin
{


template <typename B>
Vector make_field_free( const std::vector<B>& prototype )
{
    using namespace petsc;
    // we make the assumption that all basis types have "e" as the energy
    // parameter, and it is always complex
    Vector m( prototype.size() );

    populate_vector( m, [&prototype]( unsigned i ) { return prototype[i].e; } );

    return m;
}

template <typename B, typename Scalar>
Matrix make_dipole_matrix( BasisParameters bparams, std::vector<B> prototype )
{
    using namespace std;
    using namespace petsc;
    Matrix m( static_cast<unsigned>(prototype.size() ));
    auto dipole_selection_rules = [&]( unsigned i, unsigned j ) {
        return ( abs( static_cast<int>( prototype[i].l - prototype[j].l ) ) ==
                     1 &&
                 abs( prototype[i].m - prototype[j].m ) <= 1 ) ||
               prototype[i] == prototype[j];
    };
    m.reserve( dipole_selection_rules );

    auto ranges = m.get_ownership_rows();
    auto rowstart = ranges[0];
    auto rowend = ranges[1];

    PetscScalar value;

    auto grid = io::import_vector_binary<Scalar>( bparams.grid_filename() );
    // integration is r^3 dr;
    Vector integrator( grid.size() - 1, Vector::type::seq );
    populate_vector( integrator, [&grid]( unsigned i ) -> PetscScalar {
        return grid[i] * ( i == 0 ? grid[i] : grid[i] - grid[i - 1] );
    } );

    Vector dr( grid.size() - 1, Vector::type::seq );
    populate_vector( dr, [&grid]( unsigned i ) -> PetscScalar {
        return ( i == 0 ? grid[i] : grid[i] - grid[i - 1] );
    } );

    BasisLoader<Scalar> bl( bparams );
    for ( PetscInt i_ = rowstart; i_ < rowend; i_++ ) {
        for ( PetscInt j_ = 0; j_ < static_cast<int>( prototype.size() );
              j_++ ) {
            unsigned i = static_cast<unsigned>( i_ );
            unsigned j = static_cast<unsigned>( j_ );
            if ( dipole_selection_rules( i, j ) ) {
                if ( i == j ) {
                    m.set_value( i_, j_, 0. );
                    continue;
                }
                // the meat:
                Vector left = bl.left( prototype[i].n, prototype[i].l );
                Vector right = bl.right( prototype[j].n, prototype[j].l );

                // fix it, see if that changes things:
                auto rpart = inner_product( left, integrator, right );
                auto angularpart =
                    math::cg_coefficient( prototype[i], prototype[j] );

                if ( prototype[i].n == 2 && prototype[j].n == 2 ) {
                    Vector::draw( {{left, right}} );
                    cout << "n = 2: " << rpart << " * " << angularpart << " = "
                         << rpart* angularpart << endl;
                    rpart =
                        inner_product( right, conjugate( integrator ), left );
                    cout << "n = 2 *: " << rpart << " * " << angularpart
                         << " = " << rpart* angularpart << endl;
                    rpart =
                        inner_product( conjugate( left ), integrator, right );
                    cout << "n = 2 **: " << rpart << " * " << angularpart
                         << " = " << rpart* angularpart << endl;
                }
                m.set_value( i_, j_, rpart * angularpart );
            }
        }
    }

    m.assemble();

    return m;
}
}
