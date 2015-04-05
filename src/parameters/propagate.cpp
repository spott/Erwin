#include <parameters/propagate.hpp>
#include <boost/program_options.hpp>

namespace Erwin
{

std::string PropagationParameters::print() const
{
    using namespace std;
    stringstream ss;
    ss << "propagate_ti=" << ti << endl;
    ss << "propagate_tf=" << tf << endl;
    ss << "propagate_dt=" << dt << endl;
    ss << "propagate_folder=" << dt << endl;
    if ( initial_wavefunction_filename )
        ss << "propagate_wavefunction=" << *initial_wavefunction_filename
           << endl;
    return ss.str();
}

void PropagationParameters::write() const
{
    std::ofstream configfile{folder + "/PropagationParameters.config"};
    configfile << print();
    configfile.close();
}

petsc::TimeStepper::times PropagationParameters::time() const
{
    return {ti, tf, dt};
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

                    ( "propagate_folder",
                      po::value<string>()->default_value( "./" ),
                      "Folder for all propagation information" );

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

    po::store( po::command_line_parser( argc, argv )
                   .options( propagate )
                   .allow_unregistered()
                   .run(),
               vm );
    po::store( po::parse_config_file( fs, propagate, true ), vm );

    po::notify( vm );

    if ( !vm["propagate_wavefunction"].empty() )
        return PropagationParameters( vm["propagate_ti"].as<double>(),
                                      vm["propagate_tf"].as<double>(),
                                      vm["propagate_dt"].as<double>(),
                                      vm["propagate_wavefunction"].as<string>(),
                                      vm["propagate_folder"].as<string>() );
    else
        return PropagationParameters( vm["propagate_ti"].as<double>(),
                                      vm["propagate_tf"].as<double>(),
                                      vm["propagate_dt"].as<double>(),
                                      vm["propagate_folder"].as<string>() );
}
}
