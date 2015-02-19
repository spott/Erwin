#include <time_independent/time_independent.hpp>
#include <utilities/math.hpp>

int main( int argc, const char** argv )
{
    using namespace petsc;
    using namespace std;
    using namespace Erwin;

    PetscContext pc( argc, argv );

    auto grid = math::make_equally_spaced_grid( 1000, 10 );
    auto potential = []( double r ) { return 1. / r; };

    SphericalHamiltonian<double> H( grid, potential, 0 );
    // Basis<SphericalHamiltonian<double>> B(
    //     static_cast<Hamiltonian<SphericalHamiltonian<double>>>( H ), 10 );
    auto B = make_basis( H, 10 );
    vector<BasisID> prototype;

    for ( int l = 0; l < 10; ++l ) {
        std::cout << "l: " << l << std::endl;
        H.l( l );
        B.find();
        B.e.print();
        B.save_basis( "stupid_string" );
        B.add_evalues( prototype );
    }
    B.save_grid( "stupid_string" );
}
