#pragma once

// stl:
#include <iostream>
#include <vector>
#include <array>

// boost:
#include <boost/serialization/access.hpp>
#include <boost/serialization/complex.hpp>

#include <petsc_cpp/Petsc.hpp>

namespace Erwin
{
template <typename iterator>
struct Range {
    iterator begin;
    iterator end;

    int size() const { return end - begin; }
};


/************************
 * BasisID:
 * Basis for simple atomic systems
 ************************/

struct BasisID {
    unsigned n, l;
    int m;
    PetscScalar e;

  private:
    friend class boost::serialization::access;
    template <class Archive>
    void serialize( Archive& ar, const unsigned int /*version*/ )
    {
        ar& n;
        ar& l;
        ar& m;
        ar& e;
    }
};


bool operator<( const BasisID& a, const BasisID& b );

bool operator==( const BasisID& a, const BasisID& b );

bool operator<=( const BasisID& a, const BasisID& b );

bool operator>( const BasisID& a, const BasisID& b );

bool operator!=( const BasisID& a, const BasisID& b );

std::istream& operator>>( std::istream& in, BasisID& b );        // input
std::ostream& operator<<( std::ostream& out, const BasisID& b ); // output
}

namespace std
{
template <>
struct hash<Erwin::BasisID> {
  public:
    size_t operator()( const Erwin::BasisID& a ) const
    {
        return static_cast<size_t>( a.n ) + 10000 * static_cast<size_t>( a.l ) +
               1000000u * static_cast<size_t>( a.m );
    }
};
}

// End BasisID definitions

/************************
 * Angular.  Part of BasisID:
 ************************/
namespace Erwin
{

struct Angular {
    unsigned l;
    int m;

    Angular( const BasisID& a ) : l( a.l ), m( a.m ) {}
};
}


/************************
 * For std::array
 ************************/
template <typename T, size_t N>
std::ostream& operator<<( std::ostream& out,
                          const std::array<T, N>& b ) // output an array
{
    for ( size_t i = 0; i < b.size(); i++ ) {
        if ( i != 0 )
            out << ", " << b[i];
        else
            out << b[i];
    }
    return out;
}

template <typename T, size_t N>
bool operator<( const std::array<T, N>& a, const std::array<T, N>& b )
{
    for ( size_t i = 0; i < N; ++i ) {
        if ( a[i] < b[i] )
            return true;
        else if ( a[i] > b[i] )
            return false;
    }
    return false;
}
// End std::array
