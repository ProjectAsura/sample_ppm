//---------------------------------------------------------------------------------------
// File : math_util.h
// Desc : Math Utility.
//---------------------------------------------------------------------------------------

#ifndef __MATH_UTIL_H__
#define __MATH_UTIL_H__

//---------------------------------------------------------------------------------------
// Includes
//---------------------------------------------------------------------------------------
#include <cmath>


//---------------------------------------------------------------------------------------
// Constant Values
//---------------------------------------------------------------------------------------
static const double     D_PI  = 3.1415926535897932384626433832795;
static const double     D_INF = 1e20;
static const double     D_EPS = 1e-4;


double halton( const int b, int j );


/////////////////////////////////////////////////////////////////////////////////////////
// Vector3 structure
/////////////////////////////////////////////////////////////////////////////////////////
struct Vector3
{
    double  x;
    double  y;
    double  z;

    Vector3()
    : x( 0.0 ), y( 0.0 ), z( 0.0 )
    { /* DO_NOTHING */ }

    Vector3( double _x, double _y, double _z )
    : x( _x ), y( _y ), z( _z )
    { /* DO_NOTHING */ }

    inline Vector3 operator + ( const Vector3& v ) const
    { return Vector3( x + v.x, y + v.y, z + v.z ); }

    inline Vector3 operator - ( const Vector3& v ) const
    { return Vector3( x - v.x, y - v.y, z - v.z ); }

    inline Vector3 operator + ( double v ) const
    { return Vector3( x + v, y + v, z + v ); }

    inline Vector3 operator - ( double v ) const
    { return Vector3( x - v, y - v, z - v ); }

    inline Vector3 operator * ( double v ) const
    { return Vector3( x * v, y * v, z * v ); }

    inline Vector3 operator / ( double v ) const
    { return Vector3( x / v, y / v, z / v ); }

    inline Vector3 operator + () const
    { return (*this); }

    inline Vector3 operator - () const
    { return Vector3( -x, -y, -z ); }

    inline Vector3& operator += ( const Vector3& v )
    { x += v.x; y += v.y; z += v.z; return (*this); }

    inline Vector3& operator -= ( const Vector3& v )
    { x -= v.x; y -= v.y; z -= v.z; return (*this); }

    inline Vector3& operator *= ( double v )
    { x *= v; y *= v; z *= v; return (*this); }
};


inline Vector3 mul( const Vector3& a, const Vector3& b )
{ return Vector3( a.x * b.x, a.y * b.y, a.z * b.z ); }

inline Vector3 normalize( const Vector3& v )
{
    double mag = sqrt( v.x * v.x + v.y * v.y + v.z * v.z );
    return Vector3( v.x / mag, v.y / mag, v.z / mag );
}

inline double dot( const Vector3& a, const Vector3& b )
{ return ( a.x * b.x + a.y * b.y + a.z * b.z ); }

inline Vector3 cross( const Vector3& a, const Vector3& b )
{
    return Vector3(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x );
}

inline Vector3 reflect( const Vector3& a, const Vector3& b )
{
    auto dot = a.x * b.x + a.y * b.y + a.z * b.z;
    return Vector3(
        a.x - 2.0 * dot * b.x,
        a.y - 2.0 * dot * b.y,
        a.z - 2.0 * dot * b.z );
}

/////////////////////////////////////////////////////////////////////////////////////////
// Ray structure
/////////////////////////////////////////////////////////////////////////////////////////
struct Ray
{
    Vector3 pos;
    Vector3 dir;

    Ray()
    : pos(), dir()
    { /* DO_NOTHING */ }

    Ray( const Vector3& _position, const Vector3& _direction )
    : pos( _position ), dir( _direction )
    { /* DO_NOTHING */ }
};


/////////////////////////////////////////////////////////////////////////////////////////
// Sphre structure
/////////////////////////////////////////////////////////////////////////////////////////
struct Sphere
{
    Vector3 pos;
    double  r;

    Sphere( const Vector3& _pos, double _radius )
    : pos( _pos ), r( _radius )
    { /* DO_NOTHING */ }

    inline double intersect( const Ray& ray ) const
    {
        auto diff = pos - ray.pos;
        auto b = dot( diff, ray.dir );
        auto det = ( b * b ) - dot( diff, diff ) + r * r;

        if (det < 0)
            return D_INF;

        det = sqrt(det);
        auto t1 = b - det;
        if ( t1 >  D_EPS )
            return t1;

        auto t2 = b + det;
        if ( t2 > D_EPS )
            return t2;

        return D_INF;
    }
};


/////////////////////////////////////////////////////////////////////////////////////////
// BoundingBox structure
/////////////////////////////////////////////////////////////////////////////////////////
struct BoundingBox
{
    Vector3 mini;
    Vector3 maxi;

    inline void merge( const Vector3& value )
    {
        mini.x = ( value.x < mini.x ) ? value.x : mini.x;
        mini.y = ( value.y < mini.y ) ? value.y : mini.y;
        mini.z = ( value.z < mini.z ) ? value.z : mini.z;

        maxi.x = ( value.x > maxi.x ) ? value.x : maxi.x;
        maxi.y = ( value.y > maxi.y ) ? value.y : maxi.y;
        maxi.z = ( value.z > maxi.z ) ? value.z : maxi.z;
    }

    inline void reset()
    {
        mini = Vector3(  D_INF,  D_INF,  D_INF );
        maxi = Vector3( -D_INF, -D_INF, -D_INF );
    }
};


/////////////////////////////////////////////////////////////////////////////////////////
// Random class
/////////////////////////////////////////////////////////////////////////////////////////
class Random
{
public:
    Random()
    { seed( 123456789 ); }

    Random( const unsigned int _seed )
    { seed( _seed ); }

    Random( const Random& v )
    : m_x( v.m_x ), m_y( v.m_y ), m_z( v.m_z ), m_w( v.m_w )
    { /* DO_NOTHING */ }

    inline void seed( const unsigned int _seed )
    {
        m_x = 123456789;
        m_y = 362436069;
        m_z = 521288629;
        m_w = ( _seed <= 0 ) ? 88675123 : _seed;
    }

    inline unsigned int get_uint()
    {
        unsigned int t = m_x ^ ( m_x << 11 );
        m_x = m_y;
        m_y = m_z;
        m_z = m_w;
        m_w = ( m_w ^ ( m_w >> 19 ) ) ^ ( t ^ ( t >> 8 ) );
        return m_w;
    }

    inline double get_double()
    { return static_cast<double>( get_uint() ) / 0xffffffffui32; }

    inline Random& operator = ( const Random& v )
    {
        m_x = v.m_x;
        m_y = v.m_y;
        m_z = v.m_z;
        m_w = v.m_w;
        return (*this);
    }

private:
    unsigned int m_x;
    unsigned int m_y;
    unsigned int m_z;
    unsigned int m_w;
};


#endif//__MATH_UTIL_H__
