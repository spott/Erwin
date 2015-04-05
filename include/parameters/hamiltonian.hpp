#pragma once

#include <parameters/basis.hpp>
#include <experimental/optional>

namespace Erwin
{

using namespace std;

struct HamiltonianParameters {

    HamiltonianParameters( string folder_,
                           BasisParameters basis_,
                           unsigned nmax_,
                           unsigned lmax_,
                           unsigned mmax_,
                           double emax_ )
        : nmax( nmax_ ), lmax( lmax_ ), mmax( mmax_ ), emax( emax_ ),
          folder( folder_ ), basis( basis_ )
    {
    }
    HamiltonianParameters( string folder_,
                           unsigned nmax_,
                           unsigned lmax_,
                           unsigned mmax_,
                           double emax_ )
        : nmax( nmax_ ), lmax( lmax_ ), mmax( mmax_ ), emax( emax_ ),
          folder( folder_ )
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
    petsc::Matrix read_dipole( MPI_Comm comm = PETSC_COMM_WORLD ) const
    {
        return petsc::binary_import_matrix( comm, petsc::Matrix::type::aij,
                                            dipole_filename() );
    }

    string field_free_filename() const
    {
        return folder + "/energy_eigenvalues_vector.dat";
    }
    void write_field_free( const petsc::Vector& rhs ) const
    {
        rhs.to_file( field_free_filename() );
    }
    petsc::Vector read_field_free( MPI_Comm comm = PETSC_COMM_WORLD ) const
    {
        return petsc::binary_import_vector( comm, field_free_filename() );
    }

    BasisID max_basis() const
    {
        return BasisID{
            nmax, lmax, static_cast<int>( mmax ), complex<double>( emax )};
    }
    string print() const;

    unsigned nmax;
    unsigned lmax;
    unsigned mmax;
    double emax;
    string folder;
    experimental::optional<BasisParameters> basis;
};


vector<BasisID> shrink_prototype( vector<BasisID> rhs,
                                  BasisID largest_inclusive );

const HamiltonianParameters make_HamiltonianParameters( int argc,
                                                        const char** argv );
}
