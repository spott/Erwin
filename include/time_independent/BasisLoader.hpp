#pragma once

#include <utilities/data_structures.hpp>
#include <utilities/io.hpp>
#include <parameters/basis.hpp>
#include <tuple>
#include <vector>
#include <boost/iostreams/device/mapped_file.hpp>
#include <petsc_cpp/Petsc.hpp>

namespace Erwin
{

using namespace std;
using namespace petsc;

template <typename Scalar>
struct BasisLoader;
template <>
struct BasisLoader<complex<double>> {
    BasisLoader( BasisParameters basis_ ) : basis( basis_ ) {}

    const Vector left( size_t n, size_t l )
    {
        auto& ll = left_l;
        auto& file = left_file;
        if ( ll == l && file.is_open() ) {
            const complex<double>* ptr =
                reinterpret_cast<const complex<double>*>( file.data() ) +
                basis.points * ( n - l - 1 );
            return Vector( ptr, basis.points, Vector::type::seq );
        } else {
            ll = l;
            if ( file.is_open() ) file.close();
            file.open( basis.l_filename_left( ll ) );
            assert( file.size() / sizeof( complex<double> ) % basis.points ==
                    0 );
            assert( file.size() / sizeof( complex<double> ) / basis.points >=
                    ( basis.nmax - l - 1 ) );
            const complex<double>* ptr =
                reinterpret_cast<const complex<double>*>( file.data() ) +
                basis.points * ( n - l - 1 );
            return Vector( ptr, basis.points, Vector::type::seq );
        }
    }
    const Vector right( size_t n, size_t l )
    {
        auto& ll = right_l;
        auto& file = right_file;
        if ( ll == l && file.is_open() ) {
            const complex<double>* ptr =
                reinterpret_cast<const complex<double>*>( file.data() ) +
                basis.points * ( n - l - 1 );
            return Vector( ptr, basis.points, Vector::type::seq );
        } else {
            ll = l;
            if ( file.is_open() ) file.close();
            file.open( basis.l_filename_right( ll ) );
            assert( file.size() / sizeof( complex<double> ) % basis.points ==
                    0 );
            assert( file.size() / sizeof( complex<double> ) / basis.points >=
                    ( basis.nmax - l - 1 ) );
            const complex<double>* ptr =
                reinterpret_cast<const complex<double>*>( file.data() ) +
                basis.points * ( n - l - 1 );
            return Vector( ptr, basis.points, Vector::type::seq );
        }
    }

  private:
    BasisParameters basis;
    boost::iostreams::mapped_file_source left_file;
    size_t left_l;
    boost::iostreams::mapped_file_source right_file;
    size_t right_l;
};


template <>
struct BasisLoader<double> {
    BasisLoader( BasisParameters basis_ ) : basis( basis_ ) {}

    Vector left( size_t n, size_t l )
    {
        auto& ll = left_l;
        auto& file = left_file;
        if ( ll == l && file.is_open() ) {
            const double* ptr = reinterpret_cast<const double*>( file.data() ) +
                                basis.points * ( n - l - 1 );
            auto a = Vector( basis.points, Vector::type::seq );
            for ( auto i = 0u; i < basis.points; ++i )
                a.set_value( static_cast<int>( i ), ptr[i] );
            a.assemble();
            return a;
        } else {
            ll = l;
            if ( file.is_open() ) file.close();
            file.open( basis.l_filename( ll ) );
            assert( file.size() / sizeof( double ) % basis.points == 0 );
            assert( file.size() / sizeof( double ) / basis.points >=
                    ( basis.nmax - l - 1 ) );
            const double* ptr = reinterpret_cast<const double*>( file.data() ) +
                                basis.points * ( n - l - 1 );
            auto a = Vector( basis.points, Vector::type::seq );
            for ( auto i = 0u; i < basis.points; ++i )
                a.set_value( static_cast<int>( i ), ptr[i] );
            a.assemble();
            return a;
        }
    }
    Vector right( size_t n, size_t l )
    {
        auto& ll = right_l;
        auto& file = right_file;
        if ( ll == l && file.is_open() ) {
            const double* ptr = reinterpret_cast<const double*>( file.data() ) +
                                basis.points * ( n - l - 1 );
            auto a = Vector( basis.points, Vector::type::seq );
            for ( auto i = 0u; i < basis.points; ++i )
                a.set_value( static_cast<int>( i ), ptr[i] );
            a.assemble();
            return a;
        } else {
            ll = l;
            if ( file.is_open() ) file.close();
            file.open( basis.l_filename( ll ) );
            assert( file.size() / sizeof( double ) % basis.points == 0 );
            assert( file.size() / sizeof( double ) / basis.points >=
                    ( basis.nmax - l - 1 ) );
            const double* ptr = reinterpret_cast<const double*>( file.data() ) +
                                basis.points * ( n - l - 1 );
            auto a = Vector( basis.points, Vector::type::seq );
            for ( auto i = 0u; i < basis.points; ++i )
                a.set_value( static_cast<int>( i ), ptr[i] );
            a.assemble();
            return a;
        }
    }

  private:
    BasisParameters basis;
    boost::iostreams::mapped_file_source left_file;
    size_t left_l;
    boost::iostreams::mapped_file_source right_file;
    size_t right_l;
};
}
