#include <boost/program_options.hpp>
#include <utilities/math.hpp>
#include <fstream>
#include <sstream>
#include <string>
#include <parameters/basis.hpp>

namespace Erwin
{

using namespace std;

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
    if ( ecs_percent ) ss << "basis_ecs_percent=" << *ecs_percent << endl;
    if ( ecs_alpha ) ss << "basis_ecs_alpha=" << *ecs_alpha << endl;
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
        "basis_rmin", po::value<double>()->default_value( 1e-6 ),
        "the minimum r value" )( "basis_points",
                                 po::value<size_t>()->default_value( 10000 ),
                                 "the number of points on the grid" )(
        "basis_nmax", po::value<unsigned>()->default_value( 100 ),
        "the maximum principle atomic number" )(
        "basis_lmax", po::value<unsigned>()->default_value( 10 ),
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

    if ( vm["basis_ecs_percent"].empty() or vm["basis_ecs_alpha"].empty() )
        return BasisParameters(
            io::absolute_path( vm["basis_folder"].as<string>() ),
            vm["basis_rmax"].as<double>(), vm["basis_rmin"].as<double>(),
            vm["basis_points"].as<size_t>(), vm["basis_nmax"].as<unsigned>(),
            vm["basis_lmax"].as<unsigned>(), vm["basis_charge"].as<double>(),
            vm["basis_atom"].as<string>() );
    else
        return BasisParameters(
            io::absolute_path( vm["basis_folder"].as<string>() ),
            vm["basis_rmax"].as<double>(), vm["basis_rmin"].as<double>(),
            vm["basis_points"].as<size_t>(), vm["basis_nmax"].as<unsigned>(),
            vm["basis_lmax"].as<unsigned>(), vm["basis_charge"].as<double>(),
            vm["basis_atom"].as<string>(), vm["basis_ecs_percent"].as<double>(),
            vm["basis_ecs_alpha"].as<double>() );
}
}
