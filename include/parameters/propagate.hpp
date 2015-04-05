#pragma once

#include <string>
#include <petsc_cpp/Petsc.hpp>
#include <experimental/optional>

namespace Erwin
{

struct PropagationParameters {

    PropagationParameters( double ti_,
                           double tf_,
                           double dt_,
                           std::string folder_ )
        : ti( ti_ ), tf( tf_ ), dt( dt_ ), folder( folder_ )
    {
    }
    PropagationParameters( double ti_,
                           double tf_,
                           double dt_,
                           std::string folder_,
                           std::string wf_filename_ )
        : ti( ti_ ), tf( tf_ ), dt( dt_ ), folder( folder_ ),
          initial_wavefunction_filename( wf_filename_ )
    {
    }

    std::string print() const;
    void write() const;
    petsc::Vector read_initial_wavefunction() const;
    petsc::TimeStepper::times time() const;


    double ti;
    double tf;
    double dt;
    std::string folder;
    std::experimental::optional<std::string> initial_wavefunction_filename;
};

const PropagationParameters make_PropagationParameters( int argc,
                                                        const char** argv );
}
