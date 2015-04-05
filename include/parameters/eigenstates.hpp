#pragma once

#include <boost/program_options.hpp>
#include <utilities/types.hpp>
#include <time_dependent/observables.hpp>

namespace Erwin
{

struct EigenstateObserver final : Observable {
    EigenstateObserver( petsc::Matrix& A,
                        unsigned dim,
                        const std::string& folder_,
                        unsigned n_steps,
                        BasisID state_ )
        : es( A,
              dim,
              petsc::EigenvalueSolver::Which::target_imag,
              petsc::EigenvalueSolver::Type::nonhermitian ),
          gs_pop_file( folder_ + "/gs_pop.dat",
                       std::ios::binary | std::ios::ate | std::ios::out ),
          draw_gs_pop( 2, A.comm() ), draw_gs( 2, A.comm() ), steps( n_steps ),
          state( state_ ), folder( folder_ )
    {
        assert( steps > 0 );
        draw_gs_pop.set_title( "instantaneous ground state population" );
        draw_gs.set_title( "instantaneous ground state population" );
        draw_gs.set_function( []( PetscScalar d, unsigned ) {
            auto x = std::pow( std::abs( d ), 2 );
            return x == 0.0 ? -30 : std::log10( x );
        } );
        // TODO: create write for eigenstates.
        // TODO: add eigenstate.
    }

    void
    operator()( const petsc::Matrix& A, const petsc::Vector& U, petsc::TimeStepper& ts );

    void
        modify( petsc::Matrix&, petsc::Vector&,petsc::Vector&, petsc::TimeStepper& ) {}

    std::string last() const;
    std::string name() const;

  private:
    petsc::EigenvalueSolver es;
    std::vector<petsc::Vector> space;
    std::ofstream gs_pop_file;
    petsc::Draw draw_gs_pop;
    petsc::Draw draw_gs;
    const unsigned steps;
    const BasisID state;
    const std::string folder;
    PetscScalar current_value;
};
}
