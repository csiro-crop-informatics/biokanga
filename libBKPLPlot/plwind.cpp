// $Id: plwind.c 12968 2014-01-29 01:42:57Z airwin $
//
//      Routines for setting up world coordinates of the current viewport.
//
// Copyright (C) 2004  Alan W. Irwin
//
// This file is part of PLplot.
//
// PLplot is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published
// by the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// PLplot is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with PLplot; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
//
#include "stdafx.h"
#include "plplotP.h"

#define  dtr    0.01745329252

//--------------------------------------------------------------------------
// void plwind()
//
// Set up world coordinates of the viewport boundaries (2d plots).
//--------------------------------------------------------------------------

void
c_plwind( PLFLT xmin, PLFLT xmax, PLFLT ymin, PLFLT ymax )
{
    PLFLT    dx, dy, mmxmi, mmxma, mmymi, mmyma;
    PLFLT    xvpwxmin, xvpwxmax, xvpwymin, xvpwymax;
    PLWindow w;

    if ( m_plsc->level < 2 )
    {
        plabort( "plwind: Please set up viewport first" );
        return;
    }

// Best to just warn and recover on bounds errors

    if ( xmin == xmax )
    {
        plwarn( "plwind: Invalid window limits in x." );
        xmin--; xmax++;
    }
    if ( ymin == ymax )
    {
        plwarn( "plwind: Invalid window limits in y." );
        ymin--; ymax++;
    }

    m_plsc->vpwxmi = xmin;
    m_plsc->vpwxma = xmax;
    m_plsc->vpwymi = ymin;
    m_plsc->vpwyma = ymax;

// The true plot window is made slightly larger than requested so that
// the end limits will be on the graph
// Get the (slightly extended) window limits.
    plP_xgvpw( &xvpwxmin, &xvpwxmax, &xvpwymin, &xvpwymax );

// Compute the scaling between coordinate systems

    dx = xvpwxmax - xvpwxmin;
    dy = xvpwymax - xvpwymin;

    m_plsc->wpxscl = ( m_plsc->vppxma - m_plsc->vppxmi ) / dx;
    m_plsc->wpxoff = ( xmax * m_plsc->vppxmi - xmin * m_plsc->vppxma ) / dx;
    m_plsc->wpyscl = ( m_plsc->vppyma - m_plsc->vppymi ) / dy;
    m_plsc->wpyoff = ( ymax * m_plsc->vppymi - ymin * m_plsc->vppyma ) / dy;

    mmxmi = plP_dcmmx( m_plsc->vpdxmi );
    mmxma = plP_dcmmx( m_plsc->vpdxma );
    mmymi = plP_dcmmy( m_plsc->vpdymi );
    mmyma = plP_dcmmy( m_plsc->vpdyma );

// Set transformation variables for world coordinates to mm

    m_plsc->wmxscl = ( mmxma - mmxmi ) / dx;
    m_plsc->wmxoff = ( xmax * mmxmi - xmin * mmxma ) / dx;
    m_plsc->wmyscl = ( mmyma - mmymi ) / dy;
    m_plsc->wmyoff = ( ymax * mmymi - ymin * mmyma ) / dy;

// Set transformation variables for world coordinates to device coords

    m_plsc->wdxscl = m_plsc->wmxscl * m_plsc->xpmm / ( m_plsc->phyxma - m_plsc->phyxmi );
    m_plsc->wdxoff = m_plsc->wmxoff * m_plsc->xpmm / ( m_plsc->phyxma - m_plsc->phyxmi );
    m_plsc->wdyscl = m_plsc->wmyscl * m_plsc->ypmm / ( m_plsc->phyyma - m_plsc->phyymi );
    m_plsc->wdyoff = m_plsc->wmyoff * m_plsc->ypmm / ( m_plsc->phyyma - m_plsc->phyymi );

// Register plot window attributes

    w.dxmi = m_plsc->vpdxmi;
    w.dxma = m_plsc->vpdxma;
    w.dymi = m_plsc->vpdymi;
    w.dyma = m_plsc->vpdyma;

    w.wxmi = xvpwxmin;
    w.wxma = xvpwxmax;
    w.wymi = xvpwymin;
    w.wyma = xvpwymax;

    plP_swin( &w );

// Go to level 3

    m_plsc->level = 3;
}

//--------------------------------------------------------------------------
// void plw3d()
//
// Set up a window for three-dimensional plotting. The data are mapped
// into a box with world coordinate size "basex" by "basey" by "height",
// with the base being symmetrically positioned about zero. Thus
// the mapping between data 3-d and world 3-d coordinates is given by:
//
//   x = xmin   =>   wx = -0.5*basex
//   x = xmax   =>   wx =  0.5*basex
//   y = ymin   =>   wy = -0.5*basey
//   y = ymax   =>   wy =  0.5*basey
//   z = zmin   =>   wz =  0.0
//   z = zmax   =>   wz =  height
//
// The world coordinate box is then viewed from position "alt"-"az",
// measured in degrees. For proper operation, 0 <= alt <= 90 degrees,
// but az can be any value.
//--------------------------------------------------------------------------

void
c_plw3d( PLFLT basex, PLFLT basey, PLFLT height, PLFLT xmin,
         PLFLT xmax, PLFLT ymin, PLFLT ymax, PLFLT zmin,
         PLFLT zmax, PLFLT alt, PLFLT az )
{
    PLFLT xmin_adjusted, xmax_adjusted, ymin_adjusted, ymax_adjusted, zmin_adjusted, zmax_adjusted, d;
    PLFLT cx, cy, saz, caz, salt, calt, zscale;

    if ( m_plsc->level < 3 )
    {
        plabort( "plw3d: Please set up 2-d window first" );
        return;
    }
    if ( basex <= 0.0 || basey <= 0.0 || height <= 0.0 )
    {
        plabort( "plw3d: Invalid world coordinate boxsize" );
        return;
    }
    if ( xmin == xmax || ymin == ymax || zmin == zmax )
    {
        plabort( "plw3d: Invalid axis range" );
        return;
    }
    if ( alt < 0.0 || alt > 90.0 )
    {
        plabort( "plw3d: Altitude must be between 0 and 90 degrees" );
        return;
    }

    d             = 1.0e-5 * ( xmax - xmin );
    xmax_adjusted = xmax + d;
    xmin_adjusted = xmin - d;
    d             = 1.0e-5 * ( ymax - ymin );
    ymax_adjusted = ymax + d;
    ymin_adjusted = ymin - d;
    d             = 1.0e-5 * ( zmax - zmin );
    zmax_adjusted = zmax + d;
    zmin_adjusted = zmin - d;
    cx            = basex / ( xmax_adjusted - xmin_adjusted );
    cy            = basey / ( ymax_adjusted - ymin_adjusted );
    zscale        = height / ( zmax_adjusted - zmin_adjusted );
    saz           = sin( dtr * az );
    caz           = cos( dtr * az );
    salt          = sin( dtr * alt );
    calt          = cos( dtr * alt );

    m_plsc->domxmi = xmin_adjusted;
    m_plsc->domxma = xmax_adjusted;
    m_plsc->domymi = ymin_adjusted;
    m_plsc->domyma = ymax_adjusted;
    m_plsc->zzscl  = zscale;
    m_plsc->ranmi  = zmin_adjusted;
    m_plsc->ranma  = zmax_adjusted;

    m_plsc->base3x = basex;
    m_plsc->base3y = basey;
    m_plsc->basecx = 0.5 * ( xmin_adjusted + xmax_adjusted );
    m_plsc->basecy = 0.5 * ( ymin_adjusted + ymax_adjusted );
// Mathematical explanation of the 3 transformations of coordinates:
// (I) Scaling:
//     x' = cx*(x-x_mid) = cx*(x-m_plsc->basecx)
//     y' = cy*(y-y_mid) = cy*(y-m_plsc->basecy)
//     z' = zscale*(z-zmin_adjusted) = zscale*(z-m_plsc->ranmi)
// (II) Rotation about z' axis clockwise by the angle of the azimut when
//      looking from the top in a right-handed coordinate system.
//     x''          x'
//     y'' =  M_1 * y'
//     z''          z'
//    where the rotation matrix M_1 (see any mathematical physics book such
//    as Mathematical Methods in the Physical Sciences by Boas) is
//    caz          -saz       0
//    saz           caz       0
//     0             0        1
// (III) Rotation about x'' axis by 90 deg - alt to bring z''' axis
//      coincident with line of sight and x''' and y''' corresponding to
//      x and y coordinates in the 2D plane of the plot.
//     x'''          x''
//     y''' =  M_2 * y''
//     z'''          z''
//    where the rotation matrix M_2 is
//     1            0         0
//     0           salt      calt
//     0          -calt      salt
// Note
//     x'''          x'
//     y''' =  M *   y'
//     z'''          z'
// where M = M_2*M_1 is given by
//          caz      -saz     0
//     salt*saz  salt*caz    calt
//    -calt*saz -calt*caz    salt
// plP_w3wcx and plP_w3wcy take the combination of the m_plsc->basecx,
// m_plsc->basecy, m_plsc->ranmi, m_plsc->cxx, m_plsc->cxy, m_plsc->cyx, m_plsc->cyy, and
// m_plsc->cyz data stored here to implement the combination of the 3
// transformations to determine x''' and y''' from x, y, and z.
//
    m_plsc->cxx = cx * caz;
    m_plsc->cxy = -cy * saz;
    m_plsc->cyx = cx * saz * salt;
    m_plsc->cyy = cy * caz * salt;
    m_plsc->cyz = zscale * calt;
    m_plsc->czx = -cx * calt * saz;
    m_plsc->czy = -cy * calt * caz;
    m_plsc->czz = zscale * salt;
}
