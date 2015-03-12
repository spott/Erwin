#include <parameters/hamiltonian.hpp>
#include <boost/program_options.hpp>
#include <fstream>
#include <sstream>
#include <string>
#include <iterator>


namespace Erwin
{
using namespace std;

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
    if ( basis )
        ss << "hamiltonian_basis_config=" << basis->folder
           << "/BasisParameters.config" << endl;
    return ss.str();
}

vector<BasisID> shrink_prototype( vector<BasisID> rhs,
                                  BasisID largest_inclusive )
{
    vector<BasisID> out;
    copy_if( rhs.begin(), rhs.end(), back_inserter( out ),
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
    hamiltonian.add_options()

        ( "hamiltonian_nmax", po::value<unsigned>()->default_value( 100 ),
          "the maximum principle atomic number" )

            ( "hamiltonian_lmax", po::value<unsigned>()->default_value( 10 ),
              "the maximum angular atomic number" )

                ( "hamiltonian_mmax", po::value<unsigned>()->default_value( 0 ),
                  "the maximum magnetic atomic number" )

                    ( "hamiltonian_emax",
                      po::value<double>()->default_value( 1000. ),
                      "the maximum energy" )

                        ( "hamiltonian_folder",
                          po::value<string>()->default_value( "./" ),
                          "the folder the hamiltonian should be saved" )

                            ( "hamiltonian_basis_config", po::value<string>(),
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
    if ( !vm["hamiltonian_basis_config"].empty() ) {
        const char* fake_commands[3] = {
            "placeholder",
            "--basis_config",
            vm["hamiltonian_basis_config"].as<string>().c_str()};
        auto basis = make_BasisParameters( 3, fake_commands );

        return HamiltonianParameters(
            io::absolute_path( vm["hamiltonian_folder"].as<string>() ), basis,
            vm["hamiltonian_nmax"].as<unsigned>(),
            vm["hamiltonian_lmax"].as<unsigned>(),
            vm["hamiltonian_mmax"].as<unsigned>(),
            vm["hamiltonian_emax"].as<double>() );
    } else
        return HamiltonianParameters(
            io::absolute_path( vm["hamiltonian_folder"].as<string>() ),
            vm["hamiltonian_nmax"].as<unsigned>(),
            vm["hamiltonian_lmax"].as<unsigned>(),
            vm["hamiltonian_mmax"].as<unsigned>(),
            vm["hamiltonian_emax"].as<double>() );
}
}
