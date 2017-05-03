/*
 * CSIRO Open Source Software License Agreement (GPLv3)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENSE for the complete license information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */

// Please note that the CBKPLPlot class is essentially just a object wrapper around the plotting functionality provided
// by the PLPlot library and full credit is given to the various contributing authors
//

#pragma once

#include "../libBKPLPlot/plplot.h"

const int cMaxPlotTitleLen = 80;	// plot titles are of this maximum length
const int cMaxPlotXTitleLen = 80;	// plot X axis titles are of this maximum length
const int cMaxPlotYTitleLen = 80;	// plot Y axis titles are of this maximum length

const int cMaxPointLabelLen = 30;	// if point labels are used then these can be at most this long excluding terminating '\0'

const int cMaxPlotSeries = 1000;		// allow at most this many series per plot
const int cMaxPlotNumPoints = 1000000;  // and this many points per series
const int cMaxLineLegends = 40;			// show legends for first this many series lines only

typedef enum TAG_eBKPlotColor {
	BKPLblack = 0, // black (default background) 
	BKPLred,	// (default foreground) 
	BKPLyellow, 
	BKPLgreen, 
	BKPLaquamarine, 
	BKPLpink, 
	BKPLwheat, 
	BKPLgrey, 
	BKPLbrown, 
	BKPLblue, 
	BKPLBlueViolet, 
	BKPLcyan, 
	BKPLturquoise, 
	BKPLmagenta, 
	BKPLsalmon, 
	BKPLwhite 
} etBKPlotColor;

// Pointer to character string specifying options for this axis. The string can include any combination of the following letters (upper or lower case) in any order: 
//layout options are one or more of the following chars -
// a: Draws axis, X-axis is horizontal line (y=0), and Y-axis is vertical line (x=0). 
// b: Draws bottom (X) or left (Y) edge of frame. 
// c: Draws top (X) or right (Y) edge of frame. 
// d: Plot labels as date / time. Values are assumed to be seconds since the epoch (as used by gmtime). 
// f: Always use fixed point numeric labels. 
// g: Draws a grid at the major tick interval. 
// h: Draws a grid at the minor tick interval. 
// i: Inverts tick marks, so they are drawn outwards, rather than inwards. 
// l: Labels axis logarithmically. This only affects the labels, not the data, and so it is necessary to compute the logarithms of data points before passing them to any of the drawing routines. 
// m: Writes numeric labels at major tick intervals in the unconventional location (above box for X, right of box for Y). 
// n: Writes numeric labels at major tick intervals in the conventional location (below box for X, left of box for Y). 
// o: Use custom labelling function to generate axis label text. The custom labelling function can be defined with the plslabelfunc command. 
// s: Enables subticks between major ticks, only valid if t is also specified. 
// t: Draws major ticks. 
// u: Exactly like "b" except don't draw edge line. 
// w: Exactly like "c" except don't draw edge line. 
// x: Exactly like "t" (including the side effect of the numerical labels for the major ticks) except exclude drawing the major and minor tick marks. 
// v: Only applies to Y-axis == Write numeric labels for vertical axis parallel to the base of the graph, rather than parallel to the axis



typedef enum {
// The "just" parameter control how the axes will be scaled:
//
//       just=-1 : The scales will not be set, the user must set up the scale
//                   before calling plenv() using plsvpa(), plvasp() or other;
//       just= 0 : The scales will be set up to optimize plot area;
//       just= 1 : The scales will be the same;
//       just= 2 : The axes will be equal, the plot box will be square.
//
} etPlotAxisRelScale;

typedef enum {
//	axis=-2 : draw no box, no tick marks, no numeric tick labels, no axes.
//	axis=-1 : draw box only.
//	axis= 0 : Draw box, ticks, and numeric tick labels.
//	axis= 1 : Also draw coordinate axes at X=0, and Y=0.
//	axis= 2 : Also draw a grid at major tick positions in both coordinates.
//	axis= 3 : Same as 2, but the grid will be also at the minor ticks.
//	axis=10 : Same as 0 except Logarithmic X tick marks. (The X data have
//      to be converted to logarithms separately.)
//	axis=11 : Same as 1 except Logarithmic X tick marks. (The X data have
//      to be converted to logarithms separately.)
//	axis=12 : Same as 2 except Logarithmic X tick marks. (The X data have
//      to be converted to logarithms separately.)
//      axis=13 : Same as 12, but the grid will be also at the minor ticks.
//	axis=20 : Same as 0 except Logarithmic Y tick marks. (The Y data have
//      to be converted to logarithms separately.)
//	axis=21 : Same as 1 except Logarithmic Y tick marks. (The Y data have
//      to be converted to logarithms separately.)
//	axis=22 : Same as 2 except Logarithmic Y tick marks. (The Y data have
//      to be converted to logarithms separately.)
//      axis=23 : Same as 22, but the grid will be also at the minor ticks.
//	axis=30 : Same as 0 except Logarithmic X,Y tick marks. (The X,Y data have
//      to be converted to logarithms separately.)
//	axis=31 : Same as 1 except Logarithmic X,Y tick marks. (The X,Y data have
//      to be converted to logarithms separately.)
//	axis=32 : Same as 2 except Logarithmic X,Y tick marks. (The X,Y data have
//      to be converted to logarithms separately.)
//      axis=33 : Same as 32, but the grid will be also at the minor ticks.
//	axis=40 : Same as 0 except date / time X tick marks.
//	axis=41 : Same as 1 except date / time X tick marks.
//	axis=42 : Same as 2 except date / time X tick marks.
//      axis=43 : Same as 42, but the grid will be also at the minor ticks.
//	axis=50 : Same as 0 except date / time Y tick marks.
//	axis=51 : Same as 1 except date / time Y tick marks.
//	axis=52 : Same as 2 except date / time Y tick marks.
//      axis=53 : Same as 52, but the grid will be also at the minor ticks.
//	axis=60 : Same as 0 except date / time X,Y tick marks.
//	axis=61 : Same as 1 except date / time X,Y tick marks.
//	axis=62 : Same as 2 except date / time X,Y tick marks.
//      axis=63 : Same as 62, but the grid will be also at the minor ticks.
//      axis=70 : Same as 0 except custom X,Y labels.
//      axis=71 : Same as 1 except custom X,Y labels.
//      axis=72 : Same as 2 except custom X,Y labels.
//      axis=73 : Same as 72, but the grid will be also at the minor ticks.
} ePlotAxisAs;

typedef struct TAG_sPlotSeries {
	struct TAG_sPlotSeries *pNext;  // next series
	int SeriesID;					// uniquely identifies this series
	char szTitle[cMaxPlotTitleLen+1];	// title for this series
	etBKPlotColor Colour;	  	    // colour
	char Symbol;					// line includes this symbol, '\0' if no symbol
	int Style;						// line style, 1 is solid
    int Thickness;					// line thickness
	bool bLabels;					// true if each point will have associated labels
	int SeriesNumPoints;			// actual number of points added to this series
	PLFLT *pValues;					// to hold values for this series
    char *pszLabels;				// point labels
} tsPlotSeries;

class CBKPLPlot
{
	char m_szTitle[cMaxPlotTitleLen+1];	// plot title
	char m_szXTitle[cMaxPlotXTitleLen+1]; // plot X axis title
	char m_szYTitle[cMaxPlotYTitleLen+1]; // plot Y axis title
	PLFLT m_xmin;						// min world co-ordinates X value
	PLFLT m_xmax;						// max world co-ordinates X value
	PLFLT m_ymin;						// min world co-ordinates Y value
	PLFLT m_ymax;						// max world co-ordinates Y value
    etPlotAxisRelScale m_just;
    ePlotAxisAs m_axis;

	char m_szxopt[19];					// character string specifying options for horizontal x axis
	PLFLT m_xtick;						// World coordinate interval between major ticks on the x axis. If it is set to zero, PLplot automatically generates a suitable tick interval
	int   m_nxsub;						// Number of subintervals between major x axis ticks for minor ticks. If it is set to zero, PLplot automatically generates a suitable minor tick interval
	char m_szyopt[19];					// character string specifying options for horizontal y axis
	PLFLT m_ytick;						// World coordinate interval between major ticks on the y axis. If it is set to zero, PLplot automatically generates a suitable tick interval
	int   m_nysub;						// Number of subintervals between major y axis ticks for minor ticks. If it is set to zero, PLplot automatically generates a suitable minor tick interval

	int m_ValuesPerPoint;				// number of values per point
	int m_MaxNumSeries;					// hold at most this number of series
    int m_NumSeries;                    // actual number of series added 
	int m_MaxPerSeriesPoints;			// allocate to hold at most this number of points per series
	int m_MaxSeriesNumPoints;			// max number of points added to any series
	tsPlotSeries *m_pSeries;			// linked added series
	tsPlotSeries *m_pLastAccessedSeries;	// cached ptr to last accessed series

	static int BKPLPlotErrHandler(const char *);
	
	void plfbox( PLFLT x0, PLFLT y0 );

	PLINT zdefined( PLFLT x, PLFLT y );

	void f2mnmx( PLFLT **f, PLINT nnx, PLINT nny, PLFLT *fnmin, PLFLT *fnmax ); // Returns min & max of input 2d array.


	static int m_gPLPSerialiseCreated;		// initialised to 0, and set to be non-zero after first call to CreateSerialise()

	void AcquireSerialise(void);		// serialise access to underlying PLPlot functionality - I'm not convinced that these are thread safe!
	void ReleaseSerialise(void);

public:
	CBKPLPlot();
	~CBKPLPlot();

	static int CreateSerialise(void);		// create and initialise serialisation for access into PLPlot functionality - I'm not convinced that these are thread safe!

	void Reset(void);					// reset state back to that immediately following class instantiation

	int									// < 1 if errors, 1 if initialised  
	Init(int ValuesPerPoint,			// number of values per point 
				int NumSeries,			// initialise for this max number of series
				int NumPoints,			// initialise for this number of points per series
				char *pszTitle=NULL,	// top title
				char *pszXTitle=NULL,	// X axis title
				char *pszYTitle=NULL,  // Y axis title
				ePlotAxisAs Axis = (ePlotAxisAs)3,  // axis markings and plot grid
				etPlotAxisRelScale Just = (etPlotAxisRelScale)1); // axis scale

	int									// < 1 if errors, else returned series identifer (1..cMaxPlotSeries) to use when subsequently referencing that series
	AddSeries(char *pszTitle,			// add series with this title
					bool bLabels = false,	// true if each point will have associated label
			         etBKPlotColor Colour = BKPLbrown,	// line colour
					 int Style = 1,			// line style, 1 is solid
                     int Thickness = 2);    // line thickness

	int									// returns total number of points added in this series or if negative then error 
		AddPoint(int SeriesID,			// add to this series
		   PLFLT Value,					// the value
		   char *pszLabel = NULL);		// optional label 

	int									// returns total number of points added in this series or if negative then error
		AddXYPoint(int Series,			// add to this series
			   PLFLT XValue,			// the X value
			   PLFLT YValue,            // the Y value
			   char *pszLabel = NULL);	// optional label

	int									// returns total number of points added in this series or if negative then error
		AddXYZPoint(int Series,			// add to this series
		   PLFLT XValue,				// the X value
           PLFLT YValue,                // the Y value
           PLFLT ZValue,               // the Z value
			char *pszLabel = NULL);		// optional label

	void SetWorldCoords(PLFLT xmin,	// min world co-ordinates X value
					    PLFLT xmax,	// max world co-ordinates X value
						PLFLT ymin,	// min world co-ordinates Y value
						PLFLT ymax);	// max world co-ordinates Y value

	void SetAxis(int Axis,				// 0: X axis, 1: Y axis
				PLFLT tick,             // 	World coordinate interval between major ticks on the axis. If it is set to zero, PLplot automatically generates a suitable tick interval
				int nsub,				// Number of subintervals between major axis ticks for minor ticks. If it is set to zero, PLplot automatically generates a suitable minor tick interval
				char *pszOptions      // Pointer to character string specifying options for this axis. The string can include any combination of the following letters (upper or lower case) in any order: 
										//layout options are one or more of the following chars -
										// a: Draws axis, X-axis is horizontal line (y=0), and Y-axis is vertical line (x=0). 
										// b: Draws bottom (X) or left (Y) edge of frame. 
										// c: Draws top (X) or right (Y) edge of frame. 
										// d: Plot labels as date / time. Values are assumed to be seconds since the epoch (as used by gmtime). 
										// f: Always use fixed point numeric labels. 
										// g: Draws a grid at the major tick interval. 
										// h: Draws a grid at the minor tick interval. 
										// i: Inverts tick marks, so they are drawn outwards, rather than inwards. 
										// l: Labels axis logarithmically. This only affects the labels, not the data, and so it is necessary to compute the logarithms of data points before passing them to any of the drawing routines. 
										// m: Writes numeric labels at major tick intervals in the unconventional location (above box for X, right of box for Y). 
										// n: Writes numeric labels at major tick intervals in the conventional location (below box for X, left of box for Y). 
										// o: Use custom labelling function to generate axis label text. The custom labelling function can be defined with the plslabelfunc command. 
										// s: Enables subticks between major ticks, only valid if t is also specified. 
										// t: Draws major ticks. 
										// u: Exactly like "b" except don't draw edge line. 
										// w: Exactly like "c" except don't draw edge line. 
										// x: Exactly like "t" (including the side effect of the numerical labels for the major ticks) except exclude drawing the major and minor tick marks. 
				);						// v: Only applies to Y-axis == Write numeric labels for vertical axis parallel to the base of the graph, rather than parallel to the axis


	int PlotLineGraph(char *pszFile);	    // write SVG line graph plot to this file

	int PlotBarChartGraph(char *pszFile,	// write SVG histogram plot to this file
						int SeriesID = 1);	// plot this series

	int PlotPhredScores(char *pszFile);	    // write SVG Phred score distributions graph plot to this file

};
 