#pragma once

#include <string>
#include <petsc_cpp/Petsc.hpp>
#include <experimental/optional>

namespace Erwin
{
using namespace std;

struct PropagationParameters {
    enum class zeros : short { exact, approximate, none };

    PropagationParameters( double ti, double tf, double dt, zeros z )
        : ti( ti ), tf( tf ), dt( dt ), save_zeros( z )
    {
    }
    PropagationParameters(
        double ti, double tf, double dt, string wf_filename, zeros z )
        : ti( ti ), tf( tf ), dt( dt ),
          initial_wavefunction_filename( wf_filename ), save_zeros( z )
    {
    }

    string print() const;
    void write() const;
    petsc::Vector read_initial_wavefunction() const;
    petsc::TimeStepper::time time() const;


    double ti;
    double tf;
    double dt;
    experimental::optional<string> initial_wavefunction_filename;
    zeros save_zeros{zeros::approximate};
};

std::istream& operator>>( std::istream& in, PropagationParameters::zeros& z );
std::ostream& operator<<( std::ostream& out,
                          const PropagationParameters::zeros& z );

const PropagationParameters make_PropagationParameters( int argc,
                                                        const char** argv );
}
