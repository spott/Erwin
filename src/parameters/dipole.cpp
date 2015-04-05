
#include <parameters/dipole.hpp>
#include <numeric>

namespace Erwin
{

void DipoleObserver::operator()( const petsc::Matrix&,
                                 const petsc::Vector& U,
                                 petsc::TimeStepper& ts )
{
    // two copies of U:
    auto U1 = U.duplicate();
    auto U2 = U.duplicate();
    for ( auto i = sections.cbegin(); i != sections.cend(); ++i ) {
        petsc::map( U, i->second, U1 );
        for ( auto j = i; j != sections.cend(); ++j ) {
            petsc::map( U, j->second, U2 );
            PetscScalar ip = inner_product( U1, Dipole, U2 );
            ofstreams[i->first][j->first] << ip;
            if ( i == sections.cbegin() && j == i ) current_value = ip;
        }
    }
    double t[2] = {ts.time(), ts.time()};
    double y[2] = {current_value.real(), current_value.imag()};
    draw.draw_point( t, y );
}

std::string DipoleObserver::last() const
{
    using namespace std;
    char buffer[30];
    std::snprintf( &buffer[0], 30, "(%8.3e,%8.3e)", current_value.real(),
                   current_value.imag() );
    return std::string( buffer );
}
std::string DipoleObserver::name() const { return "dipole moment:"; }


std::string DipoleParameters::print() const
{
    using namespace std;
    stringstream ss;
    if ( !n_sections.empty() ) {
        ss << "dipole_n_sections=";
        for ( auto i = 0u; i < n_sections.size() - 1; ++i )
            ss << n_sections[i] << ",";
        ss << n_sections.back() << endl;
    }
    if ( !e_sections.empty() ) {
        ss << "dipole_e_sections=";
        for ( auto i = 0u; i < e_sections.size() - 1; ++i )
            ss << e_sections[i] << ",";
        ss << e_sections.back() << endl;
    }
    ss << "propagate_folder=" << folder << endl;
    return ss.str();
}


void DipoleParameters::write() const
{
    std::ofstream configfile{folder + "/DipoleParameters.config"};
    configfile << print();
    configfile.close();
}

std::unique_ptr<Observable>
DipoleParameters::get_observer( const petsc::Matrix& Dipole,
                                const std::vector<BasisID>& prototype ) const
{
    using namespace std;
    map<string, function<PetscScalar(PetscScalar, int)>> sections;
    sections.emplace( "all", []( PetscScalar a, int ) { return a; } );
    std::string alphabet = "abcdefghijklmnopqrstuvwxyz";
    std::string ALPHABET = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    vector<string> nlabels( n_sections.size() );
    vector<string> elabels( e_sections.size() );

    generate( nlabels.begin(), nlabels.end(), [&]() {
        static unsigned i = 0;
        return alphabet.substr( i++, 1 );
    } );
    generate( elabels.begin(), elabels.end(), [&]() {
        static unsigned i = 0;
        return ALPHABET.substr( i++, 1 );
    } );


    if ( !n_sections.empty() ) {
        auto n_sections_local = n_sections;
        if ( n_sections_local.front() != 1 )
            n_sections_local.insert( n_sections_local.begin()++, 1 );

        for ( auto n = 1u; n < n_sections_local.size(); n++ )
            sections.emplace( string( nlabels[n - 1] ),
                              [
                                b = n_sections_local[n],
                                a = n_sections_local[n - 1],
                                &prototype
                              ]( PetscScalar value, unsigned j ) {
                return prototype[j].n >= a && prototype[j].n < b ? value : 0.;
            } );
    }

    if ( !e_sections.empty() ) {

        auto e_sections_local = e_sections;
        if ( e_sections_local.front() > prototype[0].e.real() )
            e_sections_local.insert( e_sections_local.begin()++,
                                     prototype[0].e.real() - 1 );

        for ( auto e = 1u; e < e_sections_local.size(); e++ )
            sections.emplace( string( elabels[e - 1] ),
                              [
                                b = e_sections_local[e],
                                a = e_sections_local[e - 1],
                                &prototype
                              ]( PetscScalar value, unsigned j ) {
                return prototype[j].e.real() >= a && prototype[j].e.real() < b
                           ? value
                           : 0.;
            } );
    }

    return unique_ptr<Observable>(
        new DipoleObserver( folder, sections, Dipole ) );
}

const DipoleParameters make_DipoleParameters( int argc, const char** argv )
{
    namespace po = boost::program_options;
    using namespace std;

    string config_filename;
    po::options_description only_command_line;
    only_command_line.add_options()( "dipole_config",
                                     po::value<string>( &config_filename ),
                                     "Dipole config file" );

    po::options_description dipole( "Dipole Options" );
    dipole.add_options()( "dipole_decomposition_n",
                          po::value<vector<unsigned>>()->multitoken(),
                          "dipole decomposition in principle atomic number" )(
        "dipole_decomposition_e", po::value<vector<double>>()->multitoken(),
        "dipole decomposition in energy" )(
        "propagate_folder", po::value<string>()->default_value( "./" ),
        "the folder the various files will be written "
        "in" );

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
                   .options( dipole )
                   .allow_unregistered()
                   .run(),
               vm );
    po::store( po::parse_config_file( fs, dipole, true ), vm );
    po::notify( vm );

    if ( vm["dipole_decomposition_n"].empty() &&
         vm["dipole_decomposition_e"].empty() )
        return DipoleParameters( vm["propagate_folder"].as<string>() );
    else if ( !vm["dipole_decomposition_n"].empty() )
        return DipoleParameters(
            vm["propagate_folder"].as<string>(),
            vm["dipole_decomposition_n"].as<vector<unsigned>>() );
    else if ( !vm["dipole_decomposition_e"].empty() )
        return DipoleParameters(
            vm["propagate_folder"].as<string>(),
            vm["dipole_decomposition_e"].as<vector<double>>() );
    else
        return DipoleParameters(
            vm["propagate_folder"].as<string>(),
            vm["dipole_decomposition_e"].as<vector<double>>(),
            vm["dipole_decomposition_n"].as<vector<unsigned>>() );
}
}
