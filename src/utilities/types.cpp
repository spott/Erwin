
#include <utilities/types.hpp>


namespace Erwin
{

bool operator<( const BasisID& a, const BasisID& b )
{
    if ( a.l < b.l )
        return true;
    else if ( a.l == b.l && a.n < b.n )
        return true;
    else if ( a.l == b.l && a.n == b.n && a.m < b.m )
        return true;
    else
        return false;
}

bool operator==( const BasisID& a, const BasisID& b )
{
    if ( a.l != b.l || a.n != b.n || a.e != b.e || a.m != b.m )
        return false;
    else
        return true;
}

bool operator<=( const BasisID& a, const BasisID& b )
{
    return ( a < b ) || ( a == b );
}

bool operator>( const BasisID& a, const BasisID& b ) { return !( a <= b ); }

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
    out << b.n << ", " << b.l << ", " << b.m << ", " << b.e;
    return out;
}
}
