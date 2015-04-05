
#include <parameters/absorber.hpp>
#include <boost/program_options.hpp>

namespace Erwin
{

void MaskAbsorber::
modify( petsc::Matrix& A, petsc::Vector& U, petsc::Vector&, petsc::TimeStepper& )
{
#if defined(DEBUG)
    auto UU = U;
    current_value = UU.norm();
    UU *= mask;
    current_value -= UU.norm();
#endif

    A.diagonal_scale(mask);
}

std::string MaskAbsorber::last() const
{
    using namespace std;
    char buffer[10];
    std::snprintf( &buffer[0], 10, "%8.3e", current_value );
    return std::string( buffer );
}
std::string MaskAbsorber::name() const { return "absorbed"; }

std::string AbsorberParameters::print() const
{
    using namespace std;
    stringstream ss;
    ss << "absorber_n_size=" << n_size << endl;
    ss << "absorber_l_size=" << l_size << endl;
    if ( atype != type::custom ) ss << "absorber_type=" << atype << endl;
    return ss.str();
}

void AbsorberParameters::write() const
{
    std::ofstream configfile{folder + "/Absorber.config"};
    configfile << print();
    configfile.close();
}

petsc::Vector AbsorberParameters::mask( const std::vector<BasisID>& prototype,
                                        petsc::Vector v ) const
{
    // find the max l and n:
    unsigned n_max = 0, l_max = 0;
    for ( const auto& b : prototype ) {
        if ( b.n > n_max ) n_max = b.n;
        if ( b.l > l_max ) l_max = b.l;
    }
    assert( n_max > n_size );
    assert( l_max > l_size );
    unsigned n_start = n_max - n_size;
    unsigned l_start = l_max - l_size;

    petsc::map( v, [&]( PetscScalar a, unsigned i ) {
        unsigned n = prototype[i].n;
        unsigned l = prototype[i].l;
        a = 1;
        if ( n > n_start )
            a *= abs_function_n( static_cast<double>( n - n_start ) /
                                 this->n_size );
        if ( l > l_start )
            a *= abs_function_l( static_cast<double>( l - l_start ) /
                                 this->l_size );
        return a;
    } );

    return v;
}

std::unique_ptr<Observable>
AbsorberParameters::get_observer( const std::vector<BasisID>& prototype,
                                  petsc::Vector v ) const
{
    return std::unique_ptr<Observable>(
        new MaskAbsorber( this->mask( prototype, v ) ) );
}

std::istream& operator>>( std::istream& in, AbsorberParameters::type& z )
{
    std::string token;
    in >> token;
    if ( token == "cos_eigth" )
        z = AbsorberParameters::type::cos_eigth;
    else if ( token == "linear" )
        z = AbsorberParameters::type::linear;
    else
        throw boost::program_options::validation_error(
            boost::program_options::validation_error::invalid_option_value );
    return in;
}
std::ostream& operator<<( std::ostream& out, const AbsorberParameters::type& z )
{
    if ( z == AbsorberParameters::type::cos_eigth )
        out << "cos_eigth";
    else if ( z == AbsorberParameters::type::linear )
        out << "linear";
    return out;
}
const AbsorberParameters make_AbsorberParameters( int argc, const char** argv )
{
    namespace po = boost::program_options;
    using namespace std;

    string config_filename;
    po::options_description only_command_line;
    only_command_line.add_options()( "absorber_config",
                                     po::value<string>( &config_filename ),
                                     "Absorber config file" );

    po::options_description absorber( "Absorber Options" );
    absorber.add_options()( "absorber_n_size",
                            po::value<unsigned>()->default_value( 100 ),
                            "absorber size in n" )

        ( "absorber_l_size", po::value<unsigned>()->default_value( 50 ),
          "absorber size in l" )

            ( "propagate_folder", po::value<string>()->default_value( "./" ),
              "the folder the various laser files should be written "
              "in" )

                ( "absorber_type",
                  po::value<AbsorberParameters::type>()->default_value(
                      AbsorberParameters::type::cos_eigth ),
                  "the type of absorber: \"cos_eigth\" and \"linear\" are the "
                  "only "
                  "valid options" );

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
                   .options( absorber )
                   .allow_unregistered()
                   .run(),
               vm );
    po::store( po::parse_config_file( fs, absorber, true ), vm );

    po::notify( vm );

    return AbsorberParameters(
        vm["propagate_folder"].as<std::string>(),
        vm["absorber_n_size"].as<unsigned>(),
        vm["absorber_l_size"].as<unsigned>(),
        vm["absorber_type"].as<AbsorberParameters::type>() );
}
}
