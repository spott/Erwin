#define SLEPC
#include <petsc_cpp/Petsc.hpp>
#include <parameters/eigenstates.hpp>

namespace Erwin
{
void EigenstateObserver::operator()( const petsc::Matrix& A,
                                     const petsc::Vector& U,
                                     petsc::TimeStepper& ts )
{
    using namespace std;
    if ( ts.step() % steps == 0 && ts.step() != 0 ) {
        // static auto A = AA;
        // A = AA;
        // A /= std::complex<double>(0,-1);
        auto uu = U;
        es.op( A );
        vector<petsc::Vector> deflation_space;

        bool found_gs = false;
        auto gs_vector = U;
        vector<complex<double>> ips;
        complex<double> max_ip;
        while ( !found_gs ) {
            space.insert( space.cbegin(), uu );
            es.set_initial_vectors( space );
            if ( deflation_space.size() > 0 )
                es.set_deflation_space( deflation_space );
            if ( !A.rank() ) printf( "space size: %lu\n", space.size() );
            es.shift_invert( complex<double>(0,-1) * state.e );
            es.print();
            es.solve();

            if ( es.num_converged() == 0 ) {
                printf("none converged, resetting.... \n");
                EPSReset( es.e_ );
                auto a = es.dimensions();
                es.dimensions( static_cast<int>( a[0]++ ),
                               static_cast<int>( a[2] * 2 ),
                               static_cast<int>( a[1] * 2 ) );
                continue;
            }
            space.clear();
            for ( auto i = es.begin(); i < es.end(); ++i ) {
                auto ip = petsc::inner_product( i->evector, U );
                // we want the initial space to be made of mostly the large
                // components.
                if ( ( ip * conj( ip ) ).real() > 1e-10 ) {
                    ips.push_back( ip );
                    space.push_back( i->evector );
                    if ( !A.rank() )
                        printf( "[%d: (%e,%e) -- (%e,%e)] ", i->nev,
                                i->evalue.real(), i->evalue.imag(), ip.real(),
                                ip.imag() );
                }
                // otherwise, we want to remove this particular eigenvector from
                // the solution
                else {
                    // remove the deflated vectors from the vector we add into
                    // the initial space
                    uu -= petsc::inner_product( i->evector, uu ) * uu;
                    deflation_space.push_back( i->evector );
                }

                // this eigenvector has a larger component than the last ones
                if ( ( max_ip * conj( max_ip ) ).real() <=
                     ( ip * conj( ip ) ).real() ) {
                    max_ip = ip;
                    current_value = ip;
                    // this isn't good... too much copying.
                    gs_vector = i->evector;
                }
            }
            // if the max component is big enough:
            if ( ( max_ip * conj( max_ip ) ).real() >= 1e-2 )
                // we don't need to restart;
                found_gs = true;
        }

        double t[2] = {ts.time(), ts.time()};
        double y[2] = {abs( current_value ), arg( current_value )};

        draw_gs_pop.draw_point( t, y );
        draw_gs.draw_vector( abs_square( gs_vector ) );
        gs_vector.to_file("gs_" + to_string(ts.step()) + "_" + to_string(ts.time()) + ".dat");
        if ( !A.rank() ) printf( "\n" );
    }
}

std::string EigenstateObserver::last() const
{
    using namespace std;
    char buffer[30];
    std::snprintf( &buffer[0], 30, "(%8.3e,%8.3e)", current_value.real(),
                   current_value.imag() );
    return std::string( buffer );
}
std::string EigenstateObserver::name() const { return "gs contribution: "; }
}
