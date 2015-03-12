#include <parameters/propagate.hpp>
#include <boost/program_options.hpp>

namespace Erwin
{

std::string PropagationParameters::print() const
{
    using namespace std;
    stringstream ss;
    ss << "propagation_ti=" << ti << endl;
    ss << "propagation_tf=" << tf << endl;
    ss << "propagation_dt=" << dt << endl;
    if ( initial_wavefunction_filename )
        ss << "propagation_wavefunction="
           << initial_wavefunction_filename.value() << endl;
    ss << "propagation_zeros=" << save_zeros << endl;
    return ss.str();
}

void PropagationParameters::write() const
{
    std::ofstream configfile{folder + "/PropagationParameters.config"};
    configfile << print();
    configfile.close();
}

std::istream& operator>>( std::istream& in, PropagationParameters::zeros& z )
{
    std::string token;
    in >> token;
    if ( token == "exact" )
        z = PropagationParameters::zeros::exact;
    else if ( token == "approximate" )
        z = PropagationParameters::zeros::approximate;
    else if ( token == "none" )
        z = PropagationParameters::zeros::none;
    else
        throw boost::program_options::validation_error(
            boost::program_options::validation_error::invalid_option_value );
    return in;
}
std::ostream& operator<<( std::ostream& out, const PropagationParameters::zeros& z )
{
    if ( z == PropagationParameters::zeros::exact )
        out << "exact";
    else if ( z == PropagationParameters::zeros::approximate )
        out << "approximate";
    else if ( z == PropagationParameters::zeros::none )
        out << "none";
    return out;
}

const PropagationParameters make_PropagationParameters( int argc,
                                                        const char** argv )
{
    namespace po = boost::program_options;
    using namespace std;

    string config_filename;
    po::options_description only_command_line;
    only_command_line.add_options()( "propagate_config",
                                     po::value<string>( &config_filename ),
                                     "Propagate config file" );

    po::options_description propagate( "Propagate Options" );
    propagate.add_options()( "propagate_ti",
                             po::value<double>()->default_value( 0 ),
                             "Initial time" )

        ( "propagate_tf", po::value<double>()->default_value( -1 ),
          "Final time" )

            ( "propagate_dt", po::value<double>()->default_value( .05 ),
              "time step" )

                ( "propagate_wavefunction", po::value<string>(),
                  "an initial wavefunction to propagate" )

                    ( "propagate_zeros",
                      po::value<PropagationParameters::zeros>()->default_value(
                          PropagationParameters::zeros::approximate ),
                      "When zeros are written out" );


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

    po::store( po::parse_config_file( fs, propagate, true ), vm );
    po::store( po::command_line_parser( argc, argv )
                   .options( propagate )
                   .allow_unregistered()
                   .run(),
               vm );

    po::notify( vm );

    if ( !vm["propagate_wavefunction"].empty() )
        return PropagationParameters(
            vm["propagate_ti"].as<double>(), vm["propagate_tf"].as<double>(),
            vm["propagate_dt"].as<double>(),
            vm["propagate_wavefunction"].as<string>(),
            vm["propagate_zeros"].as<PropagationParameters::zeros>() );
    else
        return PropagationParameters(
            vm["propagate_ti"].as<double>(), vm["propagate_tf"].as<double>(),
            vm["propagate_dt"].as<double>(),
            vm["propagate_zeros"].as<PropagationParameters::zeros>() );
}
}
