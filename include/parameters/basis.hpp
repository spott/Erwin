#pragma once

#include <utilities/io.hpp>
#include <utilities/types.hpp>

namespace Erwin
{

using namespace std;

struct BasisParameters {

    BasisParameters( string folder,
                     double rmax,
                     double rmin,
                     size_t points,
                     size_t nmax,
                     size_t lmax,
                     double charge,
                     string atom )
        : rmax( rmax ), rmin( rmin ), points( points ), nmax( nmax ),
          lmax( lmax ), charge( charge ), folder( folder ), atom( atom ),
          ecs_percent( 0 ), ecs_alpha( 0 )
    {
        assert( nmax > lmax + 1 );
    }

    BasisParameters( string folder,
                     double rmax,
                     double rmin,
                     size_t points,
                     size_t nmax,
                     size_t lmax,
                     double charge,
                     string atom,
                     double ecs_percent,
                     double ecs_alpha )
        : rmax( rmax ), rmin( rmin ), points( points ), nmax( nmax ),
          lmax( lmax ), charge( charge ), folder( folder ), atom( atom ),
          ecs_percent( ecs_percent ), ecs_alpha( ecs_alpha )
    {
        assert( nmax > lmax + 1 );
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
    size_t nmax;
    size_t lmax;
    double charge;
    string folder;
    string atom;
    // the rest of these we ignore unless we are doing ecs stuff
    double ecs_percent;
    double ecs_alpha;
};

const BasisParameters make_BasisParameters( int argc, const char** argv );
}
