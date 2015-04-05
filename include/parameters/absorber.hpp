#pragma once

#include <time_dependent/observables.hpp>
#include <utilities/types.hpp>
#include <utilities/math.hpp>

namespace Erwin
{

struct MaskAbsorber final : Observable {
    MaskAbsorber( const petsc::Vector& mask_ )
        : current_value( 0 ), mask( mask_ )
    {
    }
    void
        operator()( const petsc::Matrix& , const petsc::Vector& , petsc::TimeStepper& ) {}

    void
        modify( petsc::Matrix& A, petsc::Vector& U, petsc::Vector& F, petsc::TimeStepper& ts );

    std::string last() const;
    std::string name() const;

  private:
    double current_value;
    const petsc::Vector mask;
};

struct AbsorberParameters {

    enum class type { cos_eigth, linear, custom };

    AbsorberParameters( std::string folder_,
                        unsigned n,
                        unsigned l,
                        std::function<PetscScalar(double)> abs_func_n,
                        std::function<PetscScalar(double)> abs_func_l )
        : folder( folder_ ), n_size( n ), l_size( l ),
          abs_function_n( abs_func_n ), abs_function_l( abs_func_l )
    {
    }
    AbsorberParameters( std::string folder_, unsigned n, unsigned l, type t )
        : folder( folder_ ), n_size( n ), l_size( l ), atype( t )
    {
        switch ( t ) {
            case ( type::cos_eigth ):

                abs_function_n = []( double n_ ) {
                    return std::pow( std::cos( n_ * math::PI / 2. ), 1. / 8. );
                };
                abs_function_l = []( double l_ ) {
                    return std::pow( std::cos( l_ * math::PI / 2. ), 1. / 8. );
                };
                break;
            case ( type::linear ):
                abs_function_n = []( double n_ ) { return n_; };
                abs_function_l = []( double l_ ) { return l_; };
                break;
            case ( type::custom ):
                throw std::invalid_argument( "AbsorberParameters::type::custom "
                                             "not supported in this "
                                             "AbsorberParameters constructor" );
        }
    }

    std::string print() const;
    void write() const;


    petsc::Vector mask( const std::vector<BasisID>& prototype,
                        petsc::Vector v ) const;

    std::unique_ptr<Observable>
    get_observer( const std::vector<BasisID>& prototype,
                  petsc::Vector v ) const;

  private:
    std::string folder;

  public:
    unsigned n_size;
    unsigned l_size;
    // abs_function_n/l take a value from 0 to 1 (1 == end of boundary, 0 =
    // beginning of boundary)
    std::function<PetscScalar(double)> abs_function_n;
    std::function<PetscScalar(double)> abs_function_l;

  private:
    type atype{type::custom};
};

std::istream& operator>>( std::istream& in, AbsorberParameters::type& z );
std::ostream& operator<<( std::ostream& out,
                          const AbsorberParameters::type& z );


const AbsorberParameters make_AbsorberParameters( int argc, const char** argv );
}
