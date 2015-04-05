#pragma once

#include <utilities/io.hpp>
#include <utilities/types.hpp>
#include <experimental/optional>

namespace Erwin
{

using namespace std;

struct BasisParameters {

    BasisParameters( string folder_,
                     double rmax_,
                     experimental::optional<double> rmin_,
                     size_t points_,
                     unsigned nmax_,
                     unsigned lmax_,
                     double charge_,
                     string atom_ )
        : rmax( rmax_ ), points( points_ ), nmax( nmax_ ), lmax( lmax_ ),
          charge( charge_ ), folder( folder_ ), atom( atom_ )
    {
        assert( nmax > lmax + 1 );
        if ( rmin_ )
            rmin = *rmin_;
        else
            rmin = rmax / points;
    }

    BasisParameters( string folder_,
                     double rmax_,
                     experimental::optional<double> rmin_,
                     size_t points_,
                     unsigned nmax_,
                     unsigned lmax_,
                     double charge_,
                     string atom_,
                     double ecs_percent_,
                     double ecs_alpha_ )
        : rmax( rmax_ ), points( points_ ), nmax( nmax_ ), lmax( lmax_ ),
          charge( charge_ ), folder( folder_ ), atom( atom_ ),
          ecs_percent( ecs_percent_ ), ecs_alpha( ecs_alpha_ )
    {
        assert( nmax > lmax + 1 );
        if ( rmin_ )
            rmin = *rmin_;
        else
            rmin = rmax / points;
    }

    void write() const;
    string grid_filename() const { return folder + "/grid.dat"; }
    template <typename T>
    void write_grid( const vector<T>& grid ) const
    {
        io::export_vector_binary( grid_filename(), grid );
    }
    template <typename T>
    vector<T> read_prototype() const
    {
        return io::import_vector_binary<T>( grid_filename() );
    }

    string l_filename( size_t l ) const
    {
        return folder + "/l_" + to_string( l ) + ".dat";
    }
    string l_filename_left( size_t l ) const
    {
        return folder + "/l_" + to_string( l ) + "_l.dat";
    }
    string l_filename_right( size_t l ) const
    {
        return folder + "/l_" + to_string( l ) + "_r.dat";
    }
    string prototype_filename() const { return folder + "/prototype.dat"; }
    void write_prototype( vector<BasisID> prototype ) const
    {
        io::export_vector_binary( prototype_filename(), prototype );
    }
    vector<BasisID> read_prototype() const
    {
        return io::import_vector_binary<BasisID>( prototype_filename() );
    }
    string print() const;

    double rmax;
    double rmin;
    size_t points;
    unsigned nmax;
    unsigned lmax;
    double charge;
    string folder;
    string atom;
    // the rest of these we ignore unless we are doing ecs stuff
    experimental::optional<double> ecs_percent;
    experimental::optional<double> ecs_alpha;
};

const BasisParameters make_BasisParameters( int argc, const char** argv );
}
