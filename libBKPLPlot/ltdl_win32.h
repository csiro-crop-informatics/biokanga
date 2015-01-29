// $Id: ltdl_win32.h 11970 2011-10-14 08:30:19Z airwin $
//
//      Contains all prototypes for driver functions.
//
//  Copyright (C) 2008  Werner Smekal
//
//  This file is part of PLplot.
//
//  PLplot is free software; you can redistribute it and/or modify
//  it under the terms of the GNU Library General Public License as published
//  by the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  PLplot is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Library General Public License for more details.
//
//  You should have received a copy of the GNU Library General Public License
//  along with PLplot; if not, write to the Free Software
//  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
//
//

#ifndef __LTDL_WIN32_H__
#define __LTDL_WIN32_H__

#ifdef _WIN32
#include <windows.h>
#include "pldll.h"

struct __dlhandle
{
    HINSTANCE        hinstLib;
    struct __dlhandle* previousHandle;
};
typedef struct __dlhandle*  lt_dlhandle;
typedef void                lt_ptr;

void lt_dlinit( void );

void lt_dlexit( void );

lt_dlhandle lt_dlopenext( char* dllname );

const char* lt_dlerror();

void* lt_dlsym( lt_dlhandle dlhandle, const char* symbol );

int lt_dlmakeresident( lt_dlhandle handle );

#endif // _WIN32
#endif // __LTDL_WIN32_H__
