// $Id: plmap.c 12827 2013-12-09 13:20:01Z andrewross $
//
//      Continental Outline and Political Boundary Backgrounds
//
//      Some plots need a geographical background such as the global
//      surface temperatures or the population density.  The routine
//      plmap() will draw one of the following backgrounds: continental
//      outlines, political boundaries, the United States, and the United
//      States with the continental outlines.  The routine plmeridians()
//      will add the latitudes and longitudes to the background.  After
//      the background has been drawn, one can use a contour routine or a
//      symbol plotter to finish off the plot.
//
//      Copyright (C) 1991, 1993, 1994  Wesley Ebisuzaki
//      Copyright (C) 1994, 2000, 2001  Maurice LeBrun
//      Copyright (C) 1999  Geoffrey Furnish
//      Copyright (C) 2000, 2001, 2002  Alan W. Irwin
//      Copyright (C) 2001  Andrew Roach
//      Copyright (C) 2001, 2004  Rafael Laboissiere
//      Copyright (C) 2002  Vincent Darley
//      Copyright (C) 2003  Joao Cardoso
//
//      This file is part of PLplot.
//
//      PLplot is free software; you can redistribute it and/or modify
//      it under the terms of the GNU Library General Public License
//      as published by the Free Software Foundation; version 2 of the
//      License.
//
//      PLplot is distributed in the hope that it will be useful, but
//      WITHOUT ANY WARRANTY; without even the implied warranty of
//      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//      GNU Library General Public License for more details.
//
//      You should have received a copy of the GNU Library General Public
//      License along with this library; if not, write to the Free Software
//      Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
//      USA
//
//
#include "stdafx.h"
#define DEBUG
#define NEED_PLDEBUG

#include "plplotP.h"

//--------------------------------------------------------------------------
// void plmap(void (*mapform)(PLINT, PLFLT *, PLFLT *), const char *type,
//            PLFLT minlong, PLFLT maxlong, PLFLT minlat, PLFLT maxlat);
//
// plot continental outline in world coordinates
//
// v1.4: machine independant version
// v1.3: replaced plcontinent by plmap, added plmeridians
// v1.2: 2 arguments:  mapform, type of plot
//
// mapform(PLINT n, PLFLT *x, PLFLT *y) is a routine to transform the
// coordinate longitudes and latitudes to a plot coordinate system.  By
// using this transform, we can change from a longitude, latitude
// coordinate to a polar stereographic project, for example.  Initially,
// x[0]..[n-1] are the longitudes and y[0]..y[n-1] are the corresponding
// latitudes.  After the call to mapform(), x[] and y[] should be replaced
// by the corresponding plot coordinates.  If no transform is desired,
// mapform can be replaced by NULL.
//
// type is a character string. The value of this parameter determines the
// type of background. The possible values are,
//
//      "globe"		continental outlines
//      "usa"		USA and state boundaries
//      "cglobe"	continental outlines and countries
//      "usaglobe"	USA, state boundaries and continental outlines
// alternatively the filename of a shapefile can be passed if PLplot has
// been compiled with shapelib. In this case either the base name of the
// file can be passed or the filename including the .shp or .shx suffix.
// Only the .shp and .shx files are used.
//
// minlong, maxlong are the values of the longitude on the left and right
// side of the plot, respectively. The value of minlong must be less than
// the values of maxlong, and the values of maxlong-minlong must be less
// or equal to 360.
//
// minlat, maxlat are the minimum and maximum latitudes to be plotted on
// the background.  One can always use -90.0 and 90.0 as the boundary
// outside the plot window will be automatically eliminated.  However, the
// program will be faster if one can reduce the size of the background
// plotted.
//--------------------------------------------------------------------------

#define MAP_FILE    ".map"
#define OpenMap     plLibOpenPdfstrm
#define CloseMap    pdf_close
#define OFFSET      ( 180 * 100 )
#define SCALE       100.0
#define W_BUFSIZ    ( 32 * 1024 )

void
plmap( void ( *mapform )( PLINT, PLFLT *, PLFLT * ), const char *type,
       PLFLT minlong, PLFLT maxlong, PLFLT minlat, PLFLT maxlat )
{
plwarn( "Use of the old plplot map file format is deprecated.\nIt is recommended that the shapelib library be used to provide map support.\n" );
}

//--------------------------------------------------------------------------
// void plmeridians(void (*mapform)(PLINT, PLFLT *, PLFLT *),
//		    PLFLT dlong, PLFLT dlat, PLFLT minlong, PLFLT maxlong,
//		    PLFLT minlat, PLFLT maxlat);
//
// Plot the latitudes and longitudes on the background.  The lines
// are plotted in the current color and line style.
//
// mapform(PLINT n, PLFLT *x, PLFLT *y) is a routine to transform the
// coordinate longitudes and latitudes to a plot coordinate system.  By
// using this transform, we can change from a longitude, latitude
// coordinate to a polar stereographic project, for example.  Initially,
// x[0]..x[n-1] are the longitudes and y[0]..y[n-1] are the corresponding
// latitudes.  After the call to mapform(), x[] and y[] should be replaced
// by the corresponding plot coordinates.  If no transform is desired,
// mapform can be replaced by NULL.
//
// dlat, dlong are the interval in degrees that the latitude and longitude
// lines are to be plotted.
//
// minlong, maxlong are the values of the longitude on the left and right
// side of the plot, respectively. The value of minlong must be less than
// the values of maxlong, and the values of maxlong-minlong must be less
// or equal to 360.
//
// minlat, maxlat are the minimum and maximum latitudes to be plotted on
// the background.  One can always use -90.0 and 90.0 as the boundary
// outside the plot window will be automatically eliminated.  However, the
// program will be faster if one can reduce the size of the background
// plotted.
//--------------------------------------------------------------------------

#define NSEG    100

void
plmeridians( void ( *mapform )( PLINT, PLFLT *, PLFLT * ),
             PLFLT dlong, PLFLT dlat,
             PLFLT minlong, PLFLT maxlong, PLFLT minlat, PLFLT maxlat )
{
    PLFLT yy, xx, temp, x[2], y[2], dx, dy;

    if ( minlong > maxlong )
    {
        temp    = minlong;
        minlong = maxlong;
        maxlong = temp;
    }
    if ( minlat > maxlat )
    {
        temp   = minlat;
        minlat = maxlat;
        maxlat = temp;
    }
    dx = ( maxlong - minlong ) / NSEG;
    dy = ( maxlat - minlat ) / NSEG;

    // latitudes

    for ( yy = dlat * ceil( minlat / dlat ); yy <= maxlat; yy += dlat )
    {
        if ( mapform == NULL )
        {
            plpath( NSEG, minlong, yy, maxlong, yy );
        }
        else
        {
            for ( xx = minlong; xx < maxlong; xx += dx )
            {
                y[0] = y[1] = yy;
                x[0] = xx;
                x[1] = xx + dx;
                ( *mapform )( 2, x, y );
                plline( 2, x, y );
            }
        }
    }

    // longitudes

    for ( xx = dlong * ceil( minlong / dlong ); xx <= maxlong; xx += dlong )
    {
        if ( mapform == NULL )
        {
            plpath( NSEG, xx, minlat, xx, maxlat );
        }
        else
        {
            for ( yy = minlat; yy < maxlat; yy += dy )
            {
                x[0] = x[1] = xx;
                y[0] = yy;
                y[1] = yy + dy;
                ( *mapform )( 2, x, y );
                plline( 2, x, y );
            }
        }
    }
}

