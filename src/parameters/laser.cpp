#include <parameters/laser.hpp>


namespace Erwin
{

// to kill the weak vtable error... need to figure out why that happens
LaserEnvelope::~LaserEnvelope() {}

double sin_squared::operator()( double t, double frequency ) const
{
    using namespace std;
    return ( t < ( cycles * 2. * math::PI / frequency ) and t > 0 )
               ? sqrt( intensity ) *
                     pow( sin( frequency * t / ( cycles * 2. ) ), 2. )
               : 0.0;
}

std::string sin_squared::print() const
{
    using namespace std;
    stringstream ss;
    ss << "laser_sin_squared_cycles=" << cycles << endl;
    ss << "laser_sin_squared_intensity=" << intensity << endl;
    return ss.str();
}

boost::program_options::options_description sin_squared::options()
{
    using namespace boost::program_options;
    options_description v( "Sin Squared Laser Pulse Options" );
    v.add_options()( "laser_sin_squared_cycles",
                     value<double>( &cycles )->default_value( 10 ),
                     "the total period for a sin squared laser pulse" )(
        "laser_sin_squared_intensity",
        value<double>( &intensity )->default_value( .001 ),
        "the intensity for a sin squared laser pulse" );
    return v;
}

std::string LaserParameters::print() const
{
    using namespace std;
    stringstream ss;
    ss << "laser_frequency=" << frequency << endl;
    ss << "laser_cep=" << cep << endl;
    ss << "propagate_folder=" << folder << endl;
    ss << envelope->print();
    return ss.str();
}

void EfieldObserver::
operator()( const petsc::Matrix&, const petsc::Vector&, petsc::TimeStepper& ts )
{
    if ( interpolate_next ) {
        // find the time where ef(t) ~= 0:
        double tnext = ts.time();
        double tcurrent = ts.time() - ts.dt();
        double t = ( 2 * ts.time() - ts.dt() ) / 2.;
        while ( std::abs( efield( t ) ) > 1e-12 and tnext - t > 1e-16 and
                t - tcurrent > 1e-16 ) {
            // if sign(ef(t)) != sign(current_value) then we haven't gone
            // far enough.  advance half the distance to the next timestep
            if ( math::signum( efield( t ) ) ==
                 math::signum( current_value ) ) {
                auto tmp = t;
                t = ( tnext + t ) / 2.;
                tcurrent = tmp;
            }
            // else, we have gone to far, set t to miday between current
            // timestep and t;
            else {
                auto tmp = t;
                t = ( tcurrent + t ) / 2.;
                tnext = tmp;
            }
        }

        ts.interpolate( t )
            .to_file( folder + "/wf_" + std::to_string( zero ) + ".dat" );
        zero++;
        interpolate_next = false;
    }
    current_value = next_value;
    next_value = efield( ts.time() + ts.dt() );

    if ( math::signum( current_value ) != math::signum( next_value ) &&
         math::signum( current_value ) != 0 ) {
        interpolate_next = true;
    }

    outfile << current_value;
    double t = ts.time();
    draw.draw_point( &t, &current_value );
}

std::string EfieldObserver::last() const
{
    using namespace std;
    char buffer[12];
    std::snprintf( &buffer[0], 12, "%9.3e", current_value );
    return std::string( buffer );
}

std::string EfieldObserver::name() const { return "efield"; }

std::unique_ptr<Observable> LaserParameters::get_observer() const
{
    return std::unique_ptr<Observable>( new EfieldObserver( folder, ef ) );
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

            ( "propagate_folder", po::value<string>()->default_value( "./" ),
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

    po::store( po::command_line_parser( argc, argv )
                   .options( laser )
                   .allow_unregistered()
                   .run(),
               vm );
    po::store( po::parse_config_file( fs, laser, true ), vm );

    shared_ptr<LaserEnvelope> ptr;

    if ( vm["laser_envelope_shape"].as<string>() == "sin_squared" ) {
        ptr.reset( new sin_squared() );
        auto sins = ptr->options();
        po::store( po::command_line_parser( argc, argv )
                       .options( sins )
                       .allow_unregistered()
                       .run(),
                   vm );
        po::store( po::parse_config_file( fs, sins, true ), vm );
    }
    po::notify( vm );


    return LaserParameters( vm["laser_frequency"].as<double>(),
                            vm["laser_cep"].as<double>(),
                            vm["propagate_folder"].as<string>(), ptr );
}
}
