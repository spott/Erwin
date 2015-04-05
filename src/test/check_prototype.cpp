#include <petsc_cpp/Petsc.hpp>
#include <parameters/basis.hpp>
#include <time_independent/BasisLoader.hpp>
#include <utilities/io.hpp>

int main( int argc, const char** argv )
{
    using namespace petsc;
    using namespace std;
    using namespace Erwin;

    PetscContext pc( argc, argv );

    auto parameters = make_BasisParameters( argc, argv );
    if ( !pc.rank() ) cout << parameters.print();

    auto prototype = parameters.read_prototype();
    if ( !pc.rank() )
        for ( auto& a : prototype ) cout << a << endl;

    bool valid = false;
    BasisLoader<complex<double>> bl( parameters );
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
                auto left = bl.left( n, l );
                auto right = bl.right( n, l );
                Vector::draw( {left, right} );
                util::wait_for_key();
            }
        }
        cin.clear();
        cin.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
    }
}
