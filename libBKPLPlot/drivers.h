// $Id: drivers.h 12275 2012-11-17 22:49:34Z hbabcock $
//
//      Contains all prototypes for driver functions.
//
//  Copyright (C) 2004  Andrew Roach
//  Copyright (C) 2005  Thomas J. Duck
//  Copyright (C) 2006  Andrew Ross
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

#ifndef __DRIVERS_H__
#define __DRIVERS_H__

#include "plstrm.h"

#ifdef __cplusplus
extern "C" {
#endif


void plD_dispatch_init_jpg( PLDispatchTable *pdt );
void plD_dispatch_init_jpeg( PLDispatchTable *pdt );
void plD_dispatch_init_pbm( PLDispatchTable *pdt );
void plD_dispatch_init_png( PLDispatchTable *pdt );
void plD_dispatch_init_gif( PLDispatchTable *pdt );
void plD_dispatch_init_null( PLDispatchTable *pdt );
void plD_dispatch_init_mem( PLDispatchTable *pdt );
void plD_dispatch_init_svg( PLDispatchTable *pdt );

// Prototypes for plot buffer calls.

void plbuf_init( PLStream * );
void plbuf_line( PLStream *, short, short, short, short );
void plbuf_polyline( PLStream *, short *, short *, PLINT );
void plbuf_eop( PLStream * );
void plbuf_bop( PLStream * );
void plbuf_tidy( PLStream * );
void plbuf_state( PLStream *, PLINT );
void plbuf_esc( PLStream *, PLINT, void * );
void * plbuf_save( PLStream *, void * );
void * plbuf_switch( PLStream *, void * );
void plbuf_restore( PLStream *, void * );

void plRemakePlot( PLStream * );

#ifdef __cplusplus
}
#endif

#endif  // __DRIVERS_H__
