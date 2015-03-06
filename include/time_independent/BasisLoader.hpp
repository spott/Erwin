#pragma once

#include <utilities/data_structures.hpp>
#include <utilities/io.hpp>

namespace Erwin
{

using namespace std;
using namespace petsc;

template <typename Scalar>
struct BasisLoader {

    BasisLoader( string left_filename,
                 string right_filename,
                 size_t npoints,
                 size_t que_size = 20 )
        : left( make_async_loader( que_size,
                                   io::import_vector_by_parts_fn<double>(
                                       left_filename, npoints ) ) ),
          right( make_async_loader( que_size,
                                    io::import_vector_by_parts_fn<double>(
                                        right_filename, npoints ) ) )
    {
    }

    BasisLoader( string left_filename,
                 size_t left_start,
                 string right_filename,
                 size_t right_start,
                 size_t npoints,
                 size_t que_size = 20 )
        : left( make_async_loader( que_size,
                                   io::import_vector_by_parts_fn<double>(
                                       left_filename, npoints, left_start ) ) ),
          right(
              make_async_loader( que_size,
                                 io::import_vector_by_parts_fn<double>(
                                     right_filename, npoints, right_start ) ) )
    {
    }

    vector<double> left_vector() { return left.pop_front(); }
    vector<double> right_vector() { return right.pop_front(); }

  private:
    AsyncLoader<vector<Scalar>> left;
    AsyncLoader<vector<Scalar>> right;
};
}
