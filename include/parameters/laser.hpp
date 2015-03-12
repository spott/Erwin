#pragma once

#include <boost/program_options.hpp>
#include <utilities/math.hpp>
#include <functional>
#include <string>
#include <fstream>
#include <sstream>

namespace Erwin
{

struct LaserEnvelope {
    virtual double operator()( double t ) const = 0;
    virtual std::string print() const = 0;
    virtual boost::program_options::options_description options() = 0;
    virtual ~LaserEnvelope();
};


struct sin_squared final : LaserEnvelope {
    double operator()( double t ) const;
    std::string print() const;
    boost::program_options::options_description options();
    // virtual ~sin_squared() {}

    double intensity;
    double period;
};

struct LaserParameters {

    LaserParameters( double frequency_,
                     double cep_,
                     std::string folder_,
                     std::shared_ptr<LaserEnvelope> envelope_ )
        : frequency( frequency_ ), cep( cep_ ), folder( folder_ ),
          envelope( envelope_ )
    {
    }

    std::function<double(double)> efield() const
    {
        using namespace std;
        return [&]( double t ) {
            return ( *envelope )(t)*sin( frequency * t + cep );
        };
    }

    std::string print() const;
    void write() const;

    double frequency;
    double cep;
    std::string folder;
    std::shared_ptr<LaserEnvelope> envelope;
};


const LaserParameters make_LaserParameters( int argc, const char** argv );
}
