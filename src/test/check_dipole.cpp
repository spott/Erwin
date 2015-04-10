#include <petsc_cpp/Petsc.hpp>
#include <parameters/basis.hpp>
#include <parameters/hamiltonian.hpp>
#include <time_independent/BasisLoader.hpp>
#include <utilities/io.hpp>

int main( int argc, const char** argv )
{
    using namespace petsc;
    using namespace std;
    using namespace Erwin;

    PetscContext pc( argc, argv );

    auto parameters = make_BasisParameters( argc, argv );
    auto hparameters = make_HamiltonianParameters( argc, argv );
    if ( !pc.rank() ) cout << parameters.print();
    if ( !pc.rank() ) cout << hparameters.print();

    auto prototype = parameters.read_prototype();
    if ( !pc.rank() )
        for ( auto& a : prototype ) cout << a << endl;

    auto dipole = hparameters.read_dipole( pc.comm() );
    auto psi = dipole.get_right_vector();

    bool valid = false;
    BasisLoader<complex<double>> bl( parameters );

    Draw d1_l( 2, pc.comm() );
    Draw d1_g( 2, pc.comm() );
    Draw d2_l( 2, pc.comm() );
    Draw d2_g( 2, pc.comm() );

    d1_l.set_limits(0, 100, -5, 2); 
    d1_g.set_limits(0, 100, -5, 2); 
    d2_l.set_limits(0, 100, -5, 2); 
    d2_g.set_limits(0, 100, -5, 2); 
    while ( !valid ) {
        cout << "what l: ";
        unsigned l;
        unsigned n;
        if ( cin >> l && l <= parameters.lmax ) {
            cout << endl << "what n: ";
            if ( cin >> n && n >= l + 1 && n <= parameters.nmax ) {
                auto b = find_if(
                    prototype.begin(), prototype.end(),
                    [n, l]( BasisID a ) { return n == a.n && l == a.l; } );
                printf( "\nn = %d, l = %d (%18.16e,%18.16e)\n", n, l,
                        b->e.real(), b->e.imag() );
                psi.set_value( static_cast<int>( b - prototype.begin() ), 1. );
                auto v = dipole * psi;
                auto vv = v.to_vector();
                vector<complex<double>> l1_l;
                vector<complex<double>> l1_g;
                vector<double> g1_l;
                vector<double> g1_g;
                vector<complex<double>> l2_l;
                vector<complex<double>> l2_g;
                vector<double> g2_l;
                vector<double> g2_g;
                for ( auto a = 0u; a < vv.size(); ++a ) {
                    if ( prototype[a].l == l - 1 ) {
                        if ( vv[a].real() > 0 ) {
                            l1_g.push_back( vv[a] );
                            g1_g.push_back( prototype[a].n );
                        } else {
                            l1_l.push_back( vv[a] );
                            g1_l.push_back( prototype[a].n );
                        }
                    }
                    if ( prototype[a].l == l + 1 ) {
                        if ( vv[a].real() > 0 ) {
                            l2_g.push_back( vv[a] );
                            g2_g.push_back( prototype[a].n );
                        } else {
                            l2_l.push_back( vv[a] );
                            g2_l.push_back( prototype[a].n );
                        }
                    }
                }
                d1_l.set_title( "l - 1 < 0" );
                d1_l.set_grid( g1_l );
                d1_l.set_function( []( PetscScalar s, unsigned ) {
                    auto x = abs( s );
                    return x == 0.0 ? -5 : std::log10( x );
                } );
                d1_l.draw_vector( l1_l );
                d1_g.set_title( "l - 1 > 0" );
                d1_g.set_grid( g1_g );
                d1_g.set_function( []( PetscScalar s, unsigned ) {
                    auto x = abs( s );
                    return x == 0.0 ? -5 : std::log10( x );
                } );
                d1_g.draw_vector( l1_g );
                d2_l.set_title( "l + 1 < 0" );
                d2_l.set_grid( g2_l );
                d2_l.set_function( []( PetscScalar s, unsigned ) {
                    auto x = abs( s );
                    return x == 0.0 ? -5 : std::log10( x );
                } );
                d2_l.draw_vector( l2_l );
                d2_g.set_title( "l + 1 > 0" );
                d2_g.set_grid( g2_g );
                d2_g.set_function( []( PetscScalar s, unsigned ) {
                    auto x = abs( s );
                    return x == 0.0 ? -5 : std::log10( x );
                } );
                d2_g.draw_vector( l2_g );

                util::wait_for_key();
            }
        }
        cin.clear();
        cin.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
    }
}
