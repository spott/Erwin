#pragma once

#include <time_independent/BasisLoader.hpp>

#include <petsc_cpp/Petsc.hpp>

namespace Erwin
{

using namespace std;
using namespace petsc;

template <typename B>
Matrix make_field_free( const vector<B>& prototype )
{
    // we make the assumption that all basis types have "e" as the energy
    // parameter, and it is always complex
    Matrix m( prototype.size() );

    populate_matrix( m, []( int i, int j ) { return i == j; },
                     [&prototype]( int i, int j ) { return prototype[i].e; } );

    return m;
}

template <typename B, typename Scalar>
Matrix make_dipole_matrix( BasisParameters bparams, vector<B> prototype )
{
    Matrix m( prototype.size() );
    auto dipole_selection_rules = [&]( int i, int j ) {
        return abs( prototype[i].l - prototype[j].l ) == 1 &&
               abs( prototype[i].m - prototype[j].m ) <= 1;
    };
    m.reserve( dipole_selection_rules );

    auto ranges = m_.get_ownership_rows();
    auto rowstart = ranges[0];
    auto rowend = ranges[1];

    PetscScalar value;

    // l values that matter:
    const vector<int> l_values = [&]() {
        vector<B> l_value_tmp;
        vector<int> ls;
        unique_copy( prototype.begin() + rowstart, prototype.begin() + rowend,
                     l_value_tmp, []( auto a, auto b ) { return a.l == b.l; } );
        transform( l_value_tmp.begin(), l_value_tmp.end(), ls.begin(),
                   []( auto a ) { return a.l; } );
        return ls;
    }();

    BasisLoader<Scalar> bl( bparams.l_filename( prototype[rowstart].l ),
                            prototype[rowstart].n - prototype[rowstart].l - 1,
                            bparams.l_filename( prototype[rowstart].l - 1 ),
                            0, bparams.points );
    for ( PetscInt i = rowstart; i < rowend; i++ ) {
        for ( PetscInt j = 0u; j < prototype.size(); j++ ) {
            if ( dipole_selection_rules( i, j ) ) {
                // the meat:
                auto left = bl.left_vector();
                auto right = bl.right_vector();

                auto rpart = math::integrate(left, grid, right);
                auto angularpart = math::cg_coefficients(prototype[i], prototype[j]);

                m_.set_value(i,j,rpart * angularpart);
            }
        }
    }
}
}
