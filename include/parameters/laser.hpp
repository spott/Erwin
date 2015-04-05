#pragma once

#include <boost/program_options.hpp>
#include <utilities/math.hpp>
#include <time_dependent/observables.hpp>
#include <functional>
#include <string>
#include <fstream>
#include <sstream>

namespace Erwin
{

struct LaserEnvelope {
    virtual double operator()( double t, double frequency ) const = 0;
    virtual std::string print() const = 0;
    virtual boost::program_options::options_description options() = 0;
    virtual ~LaserEnvelope();
};


struct sin_squared final : LaserEnvelope {
    double operator()( double t, double frequency ) const;
    std::string print() const;
    boost::program_options::options_description options();

    double intensity;
    double cycles;
};

struct EfieldObserver final : Observable {
    EfieldObserver( std::string folder_,
                    const std::function<double(double)>& ef )
        : draw( 1, PETSC_COMM_WORLD ), folder( folder_ ), current_value( 0 ),
          efield( ef )
    {
        draw.set_title("Efield");
        outfile.open( folder + "/efield.dat", std::ios_base::binary |
                                                  std::ios_base::app |
                                                  std::ios_base::out );
    }
    void
    operator()(const petsc::Matrix& A,const petsc::Vector& U, petsc::TimeStepper& ts );

    void modify(petsc::Matrix&, petsc::Vector&,petsc::Vector&, petsc::TimeStepper&) {}

    std::string last() const;
    std::string name() const;

  private:
    petsc::Draw draw;
    std::string folder;
    std::ofstream outfile;
    double current_value;
    const std::function<double(double)>& efield;
    double next_value{0};
    bool interpolate_next{false};
    unsigned zero{0};
};

struct LaserParameters {

    LaserParameters( double frequency_,
                     double cep_,
                     std::string folder_,
                     std::shared_ptr<LaserEnvelope> envelope_ )
        : frequency( frequency_ ), cep( cep_ ), envelope( envelope_ ),
          folder( folder_ ),
          ef( [this]( double t ) {
              return ( *this->envelope.get() )( t, this->frequency ) *
                     sin( this->frequency * t + this->cep );
          } )
    {
    }

    const std::function<double(double)>& efield() const
    {
        using namespace std;
        return ef;
    }

    std::unique_ptr<Observable> get_observer() const;

    std::string print() const;
    void write() const;

    double frequency;
    double cep;
    std::shared_ptr<LaserEnvelope> envelope;

  private:
    std::string folder;
    const std::function<double(double)> ef;
};


const LaserParameters make_LaserParameters( int argc, const char** argv );
}
