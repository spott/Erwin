#pragma once

#include <parameters/basis.hpp>

namespace Erwin
{

using namespace std;

struct HamiltonianParameters {

    HamiltonianParameters( string folder,
                           BasisParameters basis,
                           size_t nmax,
                           size_t lmax,
                           size_t mmax,
                           double emax )
        : nmax( nmax ), lmax( lmax ), mmax( mmax ), emax( emax ),
          folder( folder ), basis( basis )
    {
    }

    void write() const;
    string prototype_filename() const { return folder + "/prototype.dat"; }
    void write_prototype( const vector<BasisID>& rhs ) const
    {
        io::export_vector_binary( prototype_filename(), rhs );
    }
    vector<BasisID> read_prototype() const
    {
        return io::import_vector_binary<BasisID>( prototype_filename() );
    }

    string dipole_filename() const { return folder + "/dipole_matrix.dat"; }
    void write_dipole( const petsc::Matrix& rhs ) const
    {
        rhs.to_file( dipole_filename() );
    }
    petsc::Matrix read_dipole() const
    {
        petsc::Matrix m( petsc::Matrix::type::block_aij );
        petsc::binary_import( m, dipole_filename() );
        return m;
    }

    string field_free_filename() const
    {
        return folder + "/energy_eigenvalues_vector.dat";
    }
    void write_field_free( const petsc::Vector& rhs ) const
    {
        rhs.to_file( field_free_filename() );
    }
    petsc::Vector read_field_free() const
    {
        petsc::Vector m;
        petsc::binary_import( m, dipole_filename() );
        return m;
    }

    BasisID max_basis() const
    {
        return BasisID{
            int( nmax ), int( lmax ), int( mmax ), complex<double>( emax )};
    }
    string print() const;

    size_t nmax;
    size_t lmax;
    size_t mmax;
    double emax;
    string folder;
    BasisParameters basis;
};


vector<BasisID> shrink_prototype( vector<BasisID> rhs,
                                  BasisID largest_inclusive );

const HamiltonianParameters make_HamiltonianParameters( int argc,
                                                        const char** argv );
}
