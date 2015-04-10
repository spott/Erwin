#pragma once

#include <boost/program_options.hpp>
#include <time_dependent/observables.hpp>
#include <unordered_map>
#include <utilities/types.hpp>

namespace Erwin
{

struct DipoleObserver final : Observable {
    DipoleObserver(
        const std::string& folder_,
        std::map<std::string, std::function<PetscScalar(PetscScalar, int)>>
            sections_,
        const petsc::Matrix& d )
        : draw( 2, d.comm() ), sections( sections_ ), folder( folder_ ),
          current_value( 0 ), Dipole( d )
    {
        draw.set_title( "total dipole moment" );
        // create an ofstream for the combination of each function:
        for ( auto i = sections.cbegin(); i != sections.cend(); ++i ) {
            for ( auto j = i; j != sections.cend(); ++j ) {
                ofstreams[i->first][j->first].open(
                     folder + "/dipole_" + i->first + j->first +
                                       ".dat",
                                   std::ios_base::binary | std::ios_base::app |
                                       std::ios_base::out );
            }
        }
    }

    void
    operator()( const petsc::Matrix& A,const petsc::Vector& U, petsc::TimeStepper& ts );

    void
        modify( petsc::Matrix&, petsc::Vector&,petsc::Vector&, petsc::TimeStepper& ) {}

    std::string last() const;
    std::string name() const;

  private:
    petsc::Draw draw;
    std::unordered_map<std::string,
                       std::unordered_map<std::string, std::ofstream>>
        ofstreams;
    std::map<std::string, std::function<PetscScalar(PetscScalar, int)>>
        sections;
    const std::string folder;
    PetscScalar current_value;
    const petsc::Matrix& Dipole;
};

struct DipoleParameters {

    DipoleParameters( std::string folder_ ) : folder( folder_ ) {}

    DipoleParameters( std::string folder_, std::vector<unsigned> n )
        : n_sections( n ), folder( folder_ )
    {
    }
    DipoleParameters( std::string folder_, std::vector<double> e )
        : e_sections( e ), folder( folder_ )
    {
    }
    DipoleParameters( std::string folder_,
                      std::vector<double> e,
                      std::vector<unsigned> n )
        : n_sections( n ), e_sections( e ), folder( folder_ )
    {
    }

    std::string print() const;
    void write() const;
    std::unique_ptr<Observable>
    get_observer( const petsc::Matrix& Dipole,
                  const std::vector<BasisID>& prototype ) const;

  private:
    std::vector<unsigned> n_sections;
    std::vector<double> e_sections;
    std::string folder;
};

const DipoleParameters make_DipoleParameters( int argc, const char** argv );
}
