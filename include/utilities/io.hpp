#pragma once

// c stdlib
#include <unistd.h>

// stl
#include <iostream>
#include <fstream>
#include <cstdarg>
#include <vector>

#include <petsc_cpp/Petsc.hpp>

namespace Erwin
{

namespace io
{

    inline void wait_for_key()
    {
        std::cout << std::endl << "Press ENTER to continue..." << std::endl;
        std::cin.clear();
        std::cin.ignore( std::cin.rdbuf()->in_avail() );
        std::cin.get();
    }

    inline void printProgBar( double percent, std::ostream& os = std::cout )
    {
        percent *= 100;
        std::string bar;

        for ( size_t i = 0u; i < 50; i++ ) {
            if ( i < ( percent / 2 ) ) {
                bar.replace( i, 1, "=" );
            } else if ( i == ( size_t( percent ) / 2 ) ) {
                bar.replace( i, 1, ">" );
            } else {
                bar.replace( i, 1, " " );
            }
        }

        os << "\r"
              "[" << bar << "] ";
        os.width( 3 );
        os << percent << "%     " << std::flush;
    }


    inline std::string absolute_path( const std::string& rel_path )
    {
        if ( rel_path[0] == '.' ) {
            char* a = new char[1025];
            getcwd( a, 1025 );
            std::string cwd = std::string( a );
            delete[] a;
            return cwd.append( "/" ).append( rel_path );
        } else
            return rel_path;
    }


    inline bool file_exists( const std::string& fname )
    {
        bool ret;
        std::ifstream f( fname );
        if ( f.good() )
            ret = true;
        else
            ret = false;
        f.close();

        return ret;
    }

    inline void empty_file( const std::string& filename )
    {
        std::ofstream file;
        auto openmode = std::ios::trunc;
        file.open( filename.c_str(), openmode );
        if ( !file.is_open() || file.fail() ) {
            file.close();
            throw std::runtime_error( "Failed to erase file " + filename );
        }
        file.close();
    }

    template <typename T, typename U>
    inline void export_vector_binary( const std::string& filename,
                                      const std::vector<T>& out,
                                      const std::vector<U>& prefix )
    {
        static_assert( std::is_trivially_copyable<T>(),
                       "NO NO NO - T MUST BE TRIVIALLY COPYABLE!" );
        std::ofstream file;
        file.open( filename.c_str(), std::ios::binary | std::ios::out );
        if ( file.is_open() ) {
            if ( prefix.size() > 0 )
                file.write(
                    reinterpret_cast<const char*>( prefix.data() ),
                    static_cast<size_t>( sizeof( U ) * prefix.size() ) );
            file.write( reinterpret_cast<const char*>( &out[0] ),
                        static_cast<size_t>( sizeof( T ) * out.size() ) );
            file.close();
        } else {
            throw std::runtime_error( "error opening file " + filename +
                                      " does the folder exist?" );
        }
    }

    template <typename T, typename T2 = T, size_t block_size = 100>
    inline void export_vector_binary( const std::string& filename,
                                      const std::vector<T>& out,
                                      bool append = false )
    {
        static_assert( std::is_trivially_copyable<T>(),
                       "NO NO NO - T MUST BE TRIVIALLY COPYABLE!" );
        std::ofstream file;
        auto openmode = std::ios::binary | std::ios::out;
        if ( append ) openmode = openmode | std::ios::app;
        file.open( filename.c_str(), openmode );
        if ( file.is_open() ) {
            for ( auto i = out.begin(); i < out.end(); i += block_size ) {
                std::array<T2, block_size> ni;
                for ( auto j = i; j < ( ( out.end() - i < int( block_size ) )
                                            ? out.end()
                                            : i + block_size );
                      j++ )
                    ni[static_cast<size_t>( j - i )] = static_cast<T2>( *j );
                file.write( reinterpret_cast<const char*>( &ni ),
                            static_cast<std::streamsize>(
                                sizeof( T2 ) *
                                ( ( out.end() - i < int( block_size ) )
                                      ? static_cast<size_t>( out.end() - i )
                                      : block_size ) ) );
            }
            file.close();
        } else {
            throw std::runtime_error( "error opening file " + filename +
                                      " does the folder exist?" );
        }
    }


    template <typename T>
    std::vector<T> import_vector_binary( const std::string& filename )
    {
        static_assert( std::is_trivially_copyable<T>(),
                       "NO NO NO - T MUST BE TRIVIALLY COPYABLE!" );
        std::ifstream file;
        file.open( filename.c_str(), std::ios::binary | std::ios::ate );
        std::vector<T> vec;
        if ( file.is_open() ) {
            auto size = file.tellg();
            if ( size != 0 ) {
                vec.resize( static_cast<size_t>( size ) / sizeof( T ) );

                // make sure we haven't rounded unintentionally
                assert( static_cast<size_t>( size ) % sizeof( T ) == 0 );

                file.seekg( 0, std::ios::beg );

                file.read( reinterpret_cast<char*>( vec.data() ), size );
                file.close();
            } else {
                throw std::runtime_error( "file is empty " + filename );
            }
        } else {
            throw std::runtime_error( "file didn't open: " + filename );
        }

        return vec;
    };

    template <typename T>
    inline std::function<std::vector<T>()>
    import_vector_by_parts_fn( const std::string& filename,
                               const std::streamoff stride,
                               const std::streamoff start_ = 0 )
    {
        static_assert( std::is_trivially_copyable<T>(),
                       "NO NO NO - T MUST BE TRIVIALLY COPYABLE!" );
        using namespace std;
        ifstream file( filename.c_str(), ios::binary | ios::in );
        std::streamoff start{start_ * stride};
        std::streamoff end = [&]() {
            file.seekg( 0, ios_base::seekdir::end );
            return file.tellg();
        }();
        return [
            file( std::move( file ) ),
            stride,
            start,
            end,
            filename
        ]() mutable->vector<T>
        {

            if ( file.is_open() ) {
                if ( start >= end )
                    throw std::out_of_range(
                        "tried to read a stride out of range " + filename +
                        " @ " + to_string( start ) );
                vector<T> v( stride );
                file.seekg( start );
                file.read( reinterpret_cast<char*>( v.data() ),
                           stride * sizeof( T ) );

                start += stride * sizeof( T );
                return v;
            } else {

                throw std::runtime_error( "error opening file " + filename +
                                          " does the folder exist?" );
            }
        };
    }

    template <typename T>
    inline std::function<petsc::Vector()>
    import_SeqVector_by_parts_fn( const std::string& filename,
                                  const std::streamoff stride,
                                  const std::streamoff start_ = 0 )
    {
        using namespace std;
        ifstream file( filename.c_str(), ios::binary | ios::in );
        std::streamoff start{start_ * stride};
        std::streamoff end = [&]() {
            file.seekg( 0, ios_base::seekdir::end );
            return file.tellg();
        }();

        return [
            file( std::move( file ) ),
            stride,
            start,
            end,
            filename
        ]() mutable->petsc::Vector
        {
            if ( file.is_open() ) {
                if ( start >= end )
                    throw std::out_of_range(
                        "tried to read a stride out of range " + filename +
                        " @ " + to_string( start ) );
                auto v = make_unique<vector<complex<double>>>( stride );
                file.seekg( start );
                for ( auto i = 0; i < stride; ++i )
                    file.read( reinterpret_cast<char*>( v->data() + i ),
                               sizeof( T ) );

                start += stride * sizeof( T );
                return petsc::Vector( move( v ), petsc::Vector::type::seq );
            } else {

                throw std::runtime_error( "error opening file " + filename +
                                          " does the folder exist?" );
            }
        };
    }
}
}
