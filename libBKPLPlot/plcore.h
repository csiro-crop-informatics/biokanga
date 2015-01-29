// $Id: plcore.h 12275 2012-11-17 22:49:34Z hbabcock $
//
//      Contains declarations for core plplot data structures.  This file
//      should be included only by plcore.cpp.
//
//  Copyright (C) 2004  Andrew Roach
//  Copyright (C) 2005  Thomas J. Duck
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

#ifndef __PLCORE_H__
#define __PLCORE_H__

#include "plplotP.h"
#include "drivers.h"
#include "plDevs.h"
#include "disptab.h"

typedef void ( *PLDispatchInit )( PLDispatchTable *pdt );

#ifdef HAVE_LIBUNICODE
#include <unicode.h>
#endif


// Static function prototypes

static const char     *utf8_to_ucs4( const char *ptr, PLUNICODE *unichar );
static void     grline( short *, short *, PLINT );
static void     grpolyline( short *, short *, PLINT );
static void     grfill( short *, short *, PLINT );
static void     grgradient( short *, short *, PLINT );
static void     plSelectDev( void );
static void     pldi_ini( void );
static void     calc_diplt( void );
static void     calc_didev( void );
static void     calc_diori( void );
static void     calc_dimap( void );
static void     plgdevlst( const char **, const char **, int *, int );

static void     plInitDispatchTable( void );

static void     plLoadDriver( void );

// Static variables

static PLINT xscl[PL_MAXPOLY], yscl[PL_MAXPOLY];

static PLINT initfont = 1;      // initial font: extended by default

static PLINT lib_initialized = 0;

//--------------------------------------------------------------------------
// Allocate a PLStream data structure (defined in plstrm.h).
//
// This struct contains a copy of every variable that is stream dependent.
// Only the first [index=0] stream is statically allocated; the rest
// are dynamically allocated when you switch streams (yes, it is legal
// to only initialize the first element of the array of pointers).
//--------------------------------------------------------------------------

static PLStream pls0;                             // preallocated stream
static PLINT    ipls;                             // current stream number

static PLStream *pls[PL_NSTREAMS] = { &pls0 };    // Array of stream pointers



#ifdef __cplusplus
extern "C" {
#endif

// Current stream pointer.  Global, for easier access to state info
PLStream *m_plsc  = &pls0;

#ifdef __cplusplus
}
#endif

// Only now can we include this
#include "pldebug.h"

//--------------------------------------------------------------------------
// Initialize dispatch table.
//
// Each device is selected by the appropriate define, passed in from the
// makefile.  When installing plplot you may wish to exclude devices not
// present on your system in order to reduce screen clutter.
//
// If you hit a <CR> in response to the plinit() prompt, you get the FIRST
// one active below, so arrange them accordingly for your system (i.e. all
// the system-specific ones should go first, since they won't appear on
// most systems.)
//--------------------------------------------------------------------------

static PLDispatchTable **dispatch_table = 0;
static int             npldrivers       = 0;

static PLDispatchInit  static_device_initializers[] = {

#if defined ( PLD_jpg )
    plD_dispatch_init_jpg,
#endif

#if defined ( PLD_png )
    plD_dispatch_init_png,
#endif
#if defined ( PLD_jpeg )
    plD_dispatch_init_jpeg,
#endif

#if defined ( PLD_null )
    plD_dispatch_init_null,
#endif

#if defined ( PLD_svg )
    plD_dispatch_init_svg,
#endif

    NULL
};

static int             nplstaticdevices = ( sizeof ( static_device_initializers ) /
                                            sizeof ( PLDispatchInit ) ) - 1;
static int             npldynamicdevices = 0;

#endif  // __PLCORE_H__
