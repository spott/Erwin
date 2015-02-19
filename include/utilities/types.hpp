#pragma once

//stl:
#include<iostream>
#include<vector>
#include<array>

// boost:
#include<boost/serialization/access.hpp>
#include<boost/serialization/complex.hpp>

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
    PetscInt n, l, m;
    PetscScalar e;
    bool operator<( const BasisID b ) const
    {
        if ( this->l < b.l )
            return true;
        else if ( this->l == b.l && this->n < b.n )
            return true;
        else if ( this->l == b.l && this->n < b.n && this->m == b.m )
            return true;
        else
            return false;
    }

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
bool operator==( const BasisID& a, const BasisID& b )
{
    if ( a.l != b.l || a.n != b.n || a.e != b.e || a.m != b.m )
        return false;
    else
        return true;
}
bool operator!=( const BasisID& a, const BasisID& b )
{
    if ( a.l != b.l || a.n != b.n || a.e != b.e || a.m != b.m )
        return true;
    else
        return false;
}
std::istream& operator>>( std::istream& in, BasisID& b ) // input
{
    PetscReal er, ei;
    in >> b.n >> b.l >> b.m >> er >> ei;
    b.e = std::complex<double>( er, ei );
    return in;
}
std::ostream& operator<<( std::ostream& out, const BasisID& b ) // output
{
    out << b.n << ", " << b.l << ", " << b.m << ", " << b.e.real();
    return out;
}
}
  
namespace std {
template<>
struct hash< Erwin::BasisID > {
public:
  size_t operator()(const Erwin::BasisID &a) const 
    {
        return a.n + 10000 * a.l  + 1000000 * a.m;
    }
};
}

// End BasisID definitions


/************************
 * For std::array
 ************************/
template<typename T, size_t N>
std::ostream& operator<<(std::ostream &out, const std::array<T, N> &b) //output an array
{
    for (size_t i = 0; i < b.size(); i++)
    {
        if (i != 0)
            out << ", " << b[i];
        else
            out << b[i];
    }
    return out;
}

template< typename T, size_t N>
bool operator<( const std::array<T, N> &a, const std::array<T, N> &b )
{
    for (size_t i = 0; i < N; ++i)
    {
        if (a[i] < b[i])
            return true;
        else if (a[i] > b[i])
            return false;
    }
    return false;
}
// End std::array
