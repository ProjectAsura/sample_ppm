//---------------------------------------------------------------------------------------
// File : bitmap.h
// Desc : Bitmap Module.
// Copyright(c) Project Asura. All right reserved.
//---------------------------------------------------------------------------------------

#ifndef __BITMAP_H__
#define __BITMAP_H__


bool save_to_bmp(
    const char*     filename,
    const int       width,
    const int       height,
    const double*   colors,
    const double    gamma );


#endif//__BITMAP_H__
