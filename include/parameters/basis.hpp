#pragma once

#include <boost/program_options.hpp>
#include <fstream>
#include <sstream>
#include <string>

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
    }

    void write() const;
    string grid_filename() const { return folder + "/grid.dat"; }
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

void BasisParameters::write() const
{
    ofstream configfile{folder + "/BasisParameters.config"};
    configfile << print();
    configfile.close();
}

string BasisParameters::print() const
{
    stringstream ss;
    ss << "basis_rmax=" << rmax << endl;
    ss << "basis_rmin=" << rmin << endl;
    ss << "basis_nmax=" << nmax << endl;
    ss << "basis_lmax=" << lmax << endl;
    ss << "basis_points=" << points << endl;
    ss << "basis_charge=" << charge << endl;
    ss << "basis_folder=" << folder << endl;
    ss << "basis_atom=" << atom << endl;
    ss << "basis_ecs_percent=" << ecs_percent << endl;
    ss << "basis_ecs_alpha=" << ecs_alpha << endl;
    return ss.str();
}

const BasisParameters make_BasisParameters( int argc, const char** argv )
{
    namespace po = boost::program_options;

    string config_filename;
    po::options_description only_command_line;
    only_command_line.add_options()(
        "basis_config", po::value<string>( &config_filename ), "config file" );

    po::options_description basis( "Basis Options" );
    basis.add_options()( "basis_rmax",
                         po::value<double>()->default_value( 1000 ),
                         "the maximum r value" )(
        "basis_rmin", po::value<double>()->default_value( 0 ),
        "the minimum r value" )( "basis_points",
                                 po::value<size_t>()->default_value( 10000 ),
                                 "the number of points on the grid" )(
        "basis_nmax", po::value<size_t>()->default_value( 100 ),
        "the maximum principle atomic number" )(
        "basis_lmax", po::value<size_t>()->default_value( 10 ),
        "the maximum angular atomic number" )(
        "basis_charge", po::value<double>()->default_value( 0 ),
        "the charge on the nucleus" )(
        "basis_folder", po::value<string>()->default_value( "./" ),
        "the folder the basis should be saved" )(
        "basis_atom", po::value<string>()->default_value( "hydrogen" ),
        "the atom" )( "basis_ecs_percent",
                      po::value<double>()->default_value( 0 ),
                      "absorber size as percent of rmax size" )(
        "basis_ecs_alpha",
        po::value<double>()->default_value( math::PI<double> / 6. ),
        "alpha for external complex scaling" );


    po::variables_map vm;

    // we just want to get the config filename:
    po::store( po::command_line_parser( argc, argv )
                   .options( only_command_line )
                   .allow_unregistered()
                   .run(),
               vm );
    po::notify( vm );

    ifstream fs( config_filename );
    vm.clear();

    po::store( po::parse_config_file( fs, basis, true ), vm );
    po::store( po::command_line_parser( argc, argv )
                   .options( basis )
                   .allow_unregistered()
                   .run(),
               vm );
    po::notify( vm );

    if ( vm["basis_ecs_percent"].as<double>() == 0 or
         vm["basis_ecs_alpha"].as<double>() == 0 )
        return BasisParameters(
            vm["basis_folder"].as<string>(), vm["basis_rmax"].as<double>(),
            vm["basis_rmin"].as<double>(), vm["basis_points"].as<size_t>(),
            vm["basis_nmax"].as<size_t>(), vm["basis_lmax"].as<size_t>(),
            vm["basis_charge"].as<double>(), vm["basis_atom"].as<string>() );
    else
        return BasisParameters(
            vm["basis_folder"].as<string>(), vm["basis_rmax"].as<double>(),
            vm["basis_rmin"].as<double>(), vm["basis_points"].as<size_t>(),
            vm["basis_nmax"].as<size_t>(), vm["basis_lmax"].as<size_t>(),
            vm["basis_charge"].as<double>(), vm["basis_atom"].as<string>(),
            vm["basis_ecs_percent"].as<double>(),
            vm["basis_ecs_alpha"].as<double>() );
}
}
