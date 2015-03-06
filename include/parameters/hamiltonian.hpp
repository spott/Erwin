#pragma once

#include <boost/program_options.hpp>
#include <fstream>
#include <sstream>
#include <string>
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
    void write_prototype( vector<BasisID> prototype ) const
    {
        io::export_vector_binary( prototype_filename(), prototype );
    }
    vector<BasisID> read_prototype() const
    {
        return io::import_vector_binary<BasisID>( prototype_filename() );
    }
    string dipole_filename() const { return folder + "/Dipole.dat"; }
    string field_free_filename() const { return folder + "/Energy.dat"; }
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

void HamiltonianParameters::write() const
{
    ofstream configfile{folder + "/HamiltonianParameters.config"};
    configfile << print();
    configfile.close();
}

string HamiltonianParameters::print() const
{
    stringstream ss;
    ss << "hamiltonian_nmax=" << nmax << endl;
    ss << "hamiltonian_lmax=" << lmax << endl;
    ss << "hamiltonian_folder=" << folder << endl;
    ss << "hamiltonian_basis_config=" << basis.folder
       << "/BasisParameters.config" << endl;
    return ss.str();
}

vector<BasisID> shrink_prototype( vector<BasisID> rhs,
                                  BasisID largest_inclusive )
{
    vector<BasisID> out;
    copy_if( rhs.begin(), rhs.end(), out.begin(),
             [&largest_inclusive]( auto& a ) {
        return a.n <= largest_inclusive.n && a.l <= largest_inclusive.l &&
               abs( a.m ) <= largest_inclusive.m &&
               a.e.real() <= largest_inclusive.e.real();
    } );
    return out;
}

const HamiltonianParameters make_HamiltonianParameters( int argc,
                                                        const char** argv )
{
    namespace po = boost::program_options;

    string config_filename;
    po::options_description only_command_line;
    only_command_line.add_options()( "hamiltonian_config",
                                     po::value<string>( &config_filename ),
                                     "Hamiltonian config file" );

    po::options_description hamiltonian( "Hamiltonian Options" );
    hamiltonian.add_options()( "hamiltonian_nmax",
                               po::value<size_t>()->default_value( 100 ),
                               "the maximum principle atomic number" )(
        "hamiltonian_lmax", po::value<size_t>()->default_value( 10 ),
        "the maximum angular atomic number" )(
        "hamiltonian_mmax", po::value<size_t>()->default_value( 0 ),
        "the maximum magnetic atomic number" )(
        "hamiltonian_emax", po::value<size_t>()->default_value( 1000 ),
        "the maximum energy" )( "hamiltonian_folder",
                                po::value<string>()->default_value( "./" ),
                                "the folder the hamiltonian should be saved" )(
        "hamiltonian_basis_config", po::value<string>(),
        "the config file the basis is in" );


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

    po::store( po::parse_config_file( fs, hamiltonian, true ), vm );
    po::store( po::command_line_parser( argc, argv )
                   .options( hamiltonian )
                   .allow_unregistered()
                   .run(),
               vm );
    po::notify( vm );

    // now get the BasisParameters structure:
    const char* fake_commands[3] = {
        "placeholder",
        "--basis_config",
        vm["hamiltonian_basis_config"].as<string>().c_str()};
    auto basis = make_BasisParameters( 3, fake_commands );

    return HamiltonianParameters(
        vm["hamiltonian_folder"].as<string>(), basis,
        vm["basis_nmax"].as<size_t>(), vm["basis_lmax"].as<size_t>(),
        vm["basis_mmax"].as<size_t>(), vm["basis_emax"].as<size_t>() );
}
}
