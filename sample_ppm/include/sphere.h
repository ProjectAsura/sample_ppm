//------------------------------------------------------------------------------
// File : sphere.h
// Desc : Sphere Object.
//------------------------------------------------------------------------------

#ifndef __SPHERE_H__
#define __SPHERE_H__


//------------------------------------------------------------------------------
// Includes
//------------------------------------------------------------------------------
#include <math_util.h>


////////////////////////////////////////////////////////////////////////////////
// MaterialType
////////////////////////////////////////////////////////////////////////////////
enum MaterialType
{
    Matte = 0,  // Diffuse
    Mirror,     // Specular
    Glass,      // Refract
};


////////////////////////////////////////////////////////////////////////////////
// RenderSphere structure
////////////////////////////////////////////////////////////////////////////////
struct RenderSphere : public Sphere
{
    Vector3         color;
    MaterialType    type;

    RenderSphere( double r, Vector3 pos, Vector3 col, MaterialType _type )
    : Sphere( pos, r )
    , color ( col )
    , type  ( _type )
    { /* DO_NOTHING */ }
};


#endif//__SPHERE_H__