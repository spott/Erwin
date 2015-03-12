#include <parameters/laser.hpp>


namespace Erwin
{

LaserEnvelope::~LaserEnvelope() {}

double sin_squared::operator()( double t ) const
{
    using namespace std;
    return t < period
               ? intensity * pow( sin( math::PI<double> * t / period ), 2. )
               : 0.0;
}

std::string sin_squared::print() const
{
    using namespace std;
    stringstream ss;
    ss << "laser_sin_squared_period=" << period << endl;
    ss << "laser_sin_squared_intensity=" << intensity << endl;
    return ss.str();
}

boost::program_options::options_description sin_squared::options()
{
    using namespace boost::program_options;
    options_description v( "Sin Squared Laser Pulse Options" );
    v.add_options()( "laser_sin_squared_period",
                     value<double>( &period )->default_value( 10 ),
                     "the period for a sin squared laser pulse" )(
        "laser_sin_squared_intensity",
        value<double>( &intensity )->default_value( 1e14 ),
        "the intensity for a sin squared laser pulse" );
    return v;
}

std::string LaserParameters::print() const
{
    using namespace std;
    stringstream ss;
    ss << "laser_frequency=" << frequency << endl;
    ss << "laser_cep=" << cep << endl;
    ss << envelope->print();
    return ss.str();
}

void LaserParameters::write() const
{
    std::ofstream configfile{folder + "/LaserParameters.config"};
    configfile << print();
    configfile.close();
}

const LaserParameters make_LaserParameters( int argc, const char** argv )
{
    namespace po = boost::program_options;
    using namespace std;

    string config_filename;
    po::options_description only_command_line;
    only_command_line.add_options()( "laser_config",
                                     po::value<string>( &config_filename ),
                                     "Laser config file" );

    po::options_description laser( "Laser Options" );
    laser.add_options()( "laser_frequency",
                         po::value<double>()->default_value( 0.057 ),
                         "the laser frequency (in atomic units)" )

        ( "laser_cep", po::value<double>()->default_value( 0 ),
          "the carrier envelope phase" )

            ( "laser_folder", po::value<string>()->default_value( "./" ),
              "the folder the various laser files should be written "
              "in" )

                ( "laser_filename", po::value<string>(),
                  "the file the efield will be read from" )

                    ( "laser_envelope_shape",
                      po::value<string>()->default_value( "sin_squared" ),
                      "the laser envelope shape" );


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

    po::store( po::parse_config_file( fs, laser, true ), vm );
    po::store( po::command_line_parser( argc, argv )
                   .options( laser )
                   .allow_unregistered()
                   .run(),
               vm );

    shared_ptr<LaserEnvelope> ptr;

    if ( vm["laser_envelope_shape"].as<string>() == "sin_squared" ) {
        ptr.reset( new sin_squared() );
        auto sins = ptr->options();
        po::store( po::parse_config_file( fs, sins, true ), vm );
        po::store( po::command_line_parser( argc, argv )
                       .options( sins )
                       .allow_unregistered()
                       .run(),
                   vm );
    }
    po::notify( vm );


    return LaserParameters( vm["laser_frequency"].as<double>(),
                            vm["laser_cep"].as<double>(),
                            vm["laser_folder"].as<string>(), ptr );
}
}
