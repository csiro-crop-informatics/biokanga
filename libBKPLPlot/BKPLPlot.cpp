/*
 * CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENCE for the complete licence information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */

// Please note that the CBKPLPlot class is essentially just a object wrapper around the plotting functionality provided
// by the PLPlot library and full credit is given to the various contributing authors
// 
#include "stdafx.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#if _WIN32
#include <process.h>
#include "../libbiokanga/commhdrs.h"
#else
#include <sys/mman.h>
#include <pthread.h>
#include "../libbiokanga/commhdrs.h"
#endif

#include "BKPLPlot.h"


CBKPLPlot::CBKPLPlot()
{
m_pSeries = NULL;
CreateSerialise();
Reset();
}


CBKPLPlot::~CBKPLPlot()
{
Reset();
}

void 
CBKPLPlot::Reset(void)					// reset state back to that immediately following class instantiation
{
if(m_pSeries != NULL)
	{
	tsPlotSeries *pNext;
	do {
		if(m_pSeries->pValues != NULL)
			delete m_pSeries->pValues;
		if(m_pSeries->pszLabels != NULL) 
			delete m_pSeries->pszLabels;
		pNext = m_pSeries->pNext;
		delete m_pSeries;
		}
	while((m_pSeries = pNext) != NULL);
	}
m_szTitle[0] = '\0';
m_szXTitle[0] = '\0';
m_szYTitle[0] = '\0';
m_xmin = 0.0;
m_xmax = 0.0;
m_ymin = 0.0;
m_ymax = 0.0;
m_just = (etPlotAxisRelScale)0;
m_axis = (ePlotAxisAs)0;
m_ValuesPerPoint = 0;
m_MaxNumSeries = 0;
m_NumSeries = 0; 
m_MaxPerSeriesPoints = 0;
m_MaxSeriesNumPoints = 0;

m_szxopt[0] = '\0';	
m_xtick = 0.0;	
m_nxsub = 0;
m_szyopt[0] = '\0';	
m_ytick = 0.0;
m_nysub = 0;

m_ValuesPerPoint = 0;
m_MaxNumSeries = 0;
m_NumSeries = 0;
m_MaxPerSeriesPoints = 0;
m_MaxSeriesNumPoints = 0;
m_pLastAccessedSeries = NULL; 
}


#ifdef _WIN32
static CRITICAL_SECTION m_gPLPCritSect;
#else
static pthread_spinlock_t m_ghPLPSpinLock;
#endif

int CBKPLPlot::m_gPLPSerialiseCreated = 0;   // initialised to 0, and set to be non-zero after first call to CreateSerialise()

int
CBKPLPlot::CreateSerialise(void)	
{
if(m_gPLPSerialiseCreated != 0)
	{
#ifdef _WIN32
	Sleep(2000);		// short sleep in case 1st thread to call CreateSerialise() is still initialising the serialisation critical section or spinlock
#else
	sleep(2);
#endif
	return(eBSFSuccess);
	}
m_gPLPSerialiseCreated = 0x0fffffff;

#ifdef _WIN32
if (!InitializeCriticalSectionAndSpinCount(&m_gPLPCritSect, 1000))
	{
#else
if(pthread_spin_init(&m_ghPLPSpinLock,PTHREAD_PROCESS_PRIVATE)!=0)
	{
#endif
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Fatal: unable to create mutex");
	return(eBSFerrInternal);
	}
return(eBSFSuccess);
}

void
CBKPLPlot::AcquireSerialise(void)
{
int SpinCnt = 1000;
#ifdef _WIN32
while(!TryEnterCriticalSection(&m_gPLPCritSect))
	{
	if (SpinCnt -= 1)
		continue;
	SwitchToThread();
	SpinCnt = 100;
	}
#else
while (pthread_spin_trylock(&m_ghPLPSpinLock) == EBUSY)
	{
	if (SpinCnt -= 1)
		continue;
	pthread_yield();
	SpinCnt = 100;
	}
#endif
}

void
CBKPLPlot::ReleaseSerialise(void)
{
#ifdef _WIN32
LeaveCriticalSection(&m_gPLPCritSect);
#else
pthread_spin_unlock(&m_ghPLPSpinLock);
#endif
}


int 
CBKPLPlot::Init(int ValuesPerPoint,		// number of values per point 
				int NumSeries,			// initialise for this max number of series
				int NumPoints,			// initialise for this number of points per series
				char *pszTitle,	// top title
				char *pszXTitle,	// X axis title
				char *pszYTitle,  // Y axis title
				ePlotAxisAs Axis,		// axis markings and plot grid
				etPlotAxisRelScale Just) // axis scale
{
if(ValuesPerPoint < 1 || ValuesPerPoint > 3 ||
	NumSeries < 1 || NumSeries > cMaxPlotSeries ||
	NumPoints < 1 || NumPoints > cMaxPlotNumPoints)
	return(eBSFerrParams);

Reset();

m_ValuesPerPoint = ValuesPerPoint;
m_MaxNumSeries = NumSeries;
m_MaxPerSeriesPoints = NumPoints;

if(pszTitle != NULL && pszTitle[0] != '\0')
	{
	strncpy(m_szTitle,pszTitle,sizeof(m_szTitle));
	m_szTitle[sizeof(m_szTitle)-1] = '\0';
	}

if(pszXTitle != NULL && pszXTitle[0] != '\0')
	{
	strncpy(m_szXTitle,pszXTitle,sizeof(m_szXTitle));
	m_szXTitle[sizeof(m_szXTitle)-1] = '\0';
	}

if(pszYTitle != NULL && pszYTitle[0] != '\0')
	{
	strncpy(m_szYTitle,pszYTitle,sizeof(m_szYTitle));
	m_szYTitle[sizeof(m_szYTitle)-1] = '\0';
	}

m_xmin = 0.0;
m_xmax = 1.0;
m_ymin = 0.0;
m_ymax = 1.0;

m_axis = Axis;
m_just = Just;

return(eBSFSuccess);
}

void 
CBKPLPlot::SetWorldCoords(PLFLT xmin,	// min world co-ordinates X value
					    PLFLT xmax,		// max world co-ordinates X value
						PLFLT ymin,		// min world co-ordinates Y value
						PLFLT ymax)		// max world co-ordinates Y value
{
m_xmin = xmin;
m_xmax = xmax;
m_ymin = ymin;
m_ymax = ymax;
}


void 
CBKPLPlot::SetAxis(int Axis,			// 0: X axis, 1: Y axis
				PLFLT tick,             // 	World coordinate interval between major ticks on the axis. If it is set to zero, PLplot automatically generates a suitable tick interval
				int nsub,				// Number of subintervals between major axis ticks for minor ticks. If it is set to zero, PLplot automatically generates a suitable minor tick interval
				char *pszOptions		// Pointer to character string specifying options for this axis. The string can include any combination of the following letters (upper or lower case) in any order: 
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
				)						// v: Only applies to Y-axis == Write numeric labels for vertical axis parallel to the base of the graph, rather than parallel to the axis
{
if(Axis == 0)
	{
	m_xtick = tick;
	m_nxsub = nsub;
	strncpy(m_szxopt,pszOptions,sizeof(m_szxopt));
	m_szxopt[sizeof(m_szxopt)-1] = '\0';
	}
else
	{
	m_ytick = tick;
	m_nysub = nsub;
	strncpy(m_szyopt,pszOptions,sizeof(m_szyopt));
	m_szyopt[sizeof(m_szyopt)-1] = '\0';
	}
}

// AddSeries
// Series must be added before any point values can be associated to that series
int											// negative if errors otherwise returned series identifer to use when subsequently referencing that series
CBKPLPlot::AddSeries(char *pszTitle,		// add series with this title
					 bool bLabels,			// true if each point will have associated labels
                     etBKPlotColor Colour,	// line colour
					 int Style,				// line style
                     int Thickness)         // line thickness
{
tsPlotSeries *pSeries;
tsPlotSeries *pLastSeries;
PLFLT *pValues;
char *pszLables;

if(pszTitle == NULL || pszTitle[0] == '\0')
	return(eBSFerrParams);

// check if series previously added, if so then return existing identifier
if(m_pSeries != NULL)
	{
	pSeries = m_pSeries;
	do {
		if(!stricmp(pSeries->szTitle,pszTitle))
			return(pSeries->SeriesID);
		pLastSeries = pSeries; 
		}
	while((pSeries = pSeries->pNext) != NULL);
	if(m_NumSeries == m_MaxNumSeries)
		return(eBSFerrMaxDatasets);
	}
if((pSeries = new tsPlotSeries)==NULL)
	return(eBSFerrMem);
memset(pSeries,0,sizeof(tsPlotSeries));
if((pValues = new PLFLT[m_ValuesPerPoint * m_MaxPerSeriesPoints])==NULL)
	{
	delete pSeries;
	return(eBSFerrMem);
	}

if(bLabels)
	{
	if((pszLables = new char [(cMaxPointLabelLen + 1) * m_MaxPerSeriesPoints])==NULL)
		{
		delete pSeries;
		delete pValues;
		return(eBSFerrMem);
		}
	}
else
	pszLables = NULL;

if(m_pSeries == NULL)
	{
	m_pSeries = pSeries;
	m_NumSeries = 1;
	}
else
	{
	pLastSeries->pNext = pSeries;	
	m_NumSeries += 1;
	}
pSeries->SeriesID = m_NumSeries;
strncpy(pSeries->szTitle,pszTitle,cMaxPlotTitleLen);
pSeries->szTitle[cMaxPlotTitleLen] = '\0';
pSeries->Colour = Colour;
pSeries->Style = Style;
pSeries->Thickness = Thickness;
pSeries->bLabels = bLabels;
pSeries->pValues = pValues;
pSeries->pszLabels = pszLables;
m_pLastAccessedSeries = NULL;
return(m_NumSeries);
}


int										// returns total number of points added in this series or if negative then error 
CBKPLPlot::AddPoint(int SeriesID,		// add to this series
		   PLFLT XValue,				// the value
		   char *pszLabel)			    // optional label
{
tsPlotSeries *pSeries;
PLFLT *pValue;
char *pszLabels;

if(m_ValuesPerPoint != 1)
	return(eBSFerrDataPtType);

if(m_pSeries == NULL || SeriesID < 1 || SeriesID > m_MaxNumSeries)
	return(eBSFerrDataset);

if(m_pLastAccessedSeries == NULL || m_pLastAccessedSeries->SeriesID != SeriesID)
	{
	pSeries = m_pSeries;
	do {
		if(pSeries->SeriesID == SeriesID)
			break;
		}
	while((pSeries = pSeries->pNext)!=NULL);
	if(pSeries == NULL)
		return(eBSFerrDataset);
	m_pLastAccessedSeries = pSeries;
	}
if(m_pLastAccessedSeries->SeriesNumPoints == m_MaxPerSeriesPoints)
	return(eBSFerrMaxEntries);

pValue = &m_pLastAccessedSeries->pValues[m_pLastAccessedSeries->SeriesNumPoints];
*pValue = XValue;
if(m_pLastAccessedSeries->bLabels && m_pLastAccessedSeries->pszLabels)
	{
	pszLabels = &m_pLastAccessedSeries->pszLabels[m_pLastAccessedSeries->SeriesNumPoints * (cMaxPointLabelLen + 1)];
	if(pszLabel != NULL && pszLabel[0] != '\0')
		{
		strncpy(pszLabels,pszLabel,cMaxPointLabelLen);
		pszLabels[cMaxPointLabelLen] = '\0';
		}
	else
		sprintf(pszLabels,"P:%d",m_pLastAccessedSeries->SeriesNumPoints+1);
	}

m_pLastAccessedSeries->SeriesNumPoints += 1;
if(m_MaxSeriesNumPoints < m_pLastAccessedSeries->SeriesNumPoints)
	m_MaxSeriesNumPoints = m_pLastAccessedSeries->SeriesNumPoints;
return(m_pLastAccessedSeries->SeriesNumPoints);
}	

int										// returns total number of points added in this series or if negative then error 
CBKPLPlot::AddXYPoint(int SeriesID,		// add to this series
		   PLFLT XValue,				// the X value
           PLFLT YValue,                // the Y value
		   char *pszLabel)			    // optional label
{
tsPlotSeries *pSeries;
PLFLT *pValue;
char *pszLabels;

if(m_ValuesPerPoint != 2)
	return(eBSFerrDataPtType);

if(m_pSeries == NULL || SeriesID < 1 || SeriesID > m_MaxNumSeries)
	return(eBSFerrDataset);

if(m_pLastAccessedSeries == NULL || m_pLastAccessedSeries->SeriesID != SeriesID)
	{
	pSeries = m_pSeries;
	do {
		if(pSeries->SeriesID == SeriesID)
			break;
		}
	while((pSeries = pSeries->pNext)!=NULL);
	if(pSeries == NULL)
		return(eBSFerrDataset);
	m_pLastAccessedSeries = pSeries;
	}
if(m_pLastAccessedSeries->SeriesNumPoints == m_MaxPerSeriesPoints)
	return(eBSFerrMaxEntries);

pValue = &m_pLastAccessedSeries->pValues[m_pLastAccessedSeries->SeriesNumPoints];
*pValue = XValue;
pValue += m_MaxPerSeriesPoints;
*pValue = YValue;
if(m_pLastAccessedSeries->bLabels && m_pLastAccessedSeries->pszLabels)
	{
	pszLabels = &m_pLastAccessedSeries->pszLabels[m_pLastAccessedSeries->SeriesNumPoints * (cMaxPointLabelLen + 1)];
	if(pszLabel != NULL && pszLabel[0] != '\0')
		{
		strncpy(pszLabels,pszLabel,cMaxPointLabelLen);
		pszLabels[cMaxPointLabelLen] = '\0';
		}
	else
		sprintf(pszLabels,"P:%d",m_pLastAccessedSeries->SeriesNumPoints+1);
	}
m_pLastAccessedSeries->SeriesNumPoints += 1;
if(m_MaxSeriesNumPoints < m_pLastAccessedSeries->SeriesNumPoints)
	m_MaxSeriesNumPoints = m_pLastAccessedSeries->SeriesNumPoints;
return(m_pLastAccessedSeries->SeriesNumPoints);
}	

int										// returns total number of points added in this series or if negative then error
CBKPLPlot::AddXYZPoint(int SeriesID,	// add to this series
		   PLFLT XValue,				// the X value
           PLFLT YValue,                // the Y value
           PLFLT ZValue,                // the Z value
		   char *pszLabel)			    // optional label
{
tsPlotSeries *pSeries;
PLFLT *pValue;
char *pszLabels;

if(m_ValuesPerPoint != 3)
	return(eBSFerrDataPtType);

if(m_pSeries == NULL || SeriesID < 1 || SeriesID > m_MaxNumSeries)
	return(eBSFerrDataset);

if(m_pLastAccessedSeries == NULL || m_pLastAccessedSeries->SeriesID != SeriesID)
	{
	pSeries = m_pSeries;
	do {
		if(pSeries->SeriesID == SeriesID)
			break;
		}
	while((pSeries = pSeries->pNext)!=NULL);
	if(pSeries == NULL)
		return(eBSFerrDataset);
	m_pLastAccessedSeries = pSeries;
	}
if(m_pLastAccessedSeries->SeriesNumPoints == m_MaxPerSeriesPoints)
	return(eBSFerrMaxEntries);

pValue = &m_pLastAccessedSeries->pValues[m_pLastAccessedSeries->SeriesNumPoints];
*pValue = XValue;
pValue += m_MaxPerSeriesPoints;
*pValue = YValue;
pValue += m_MaxPerSeriesPoints;
*pValue = ZValue;
if(m_pLastAccessedSeries->bLabels && m_pLastAccessedSeries->pszLabels)
	{
	pszLabels = &m_pLastAccessedSeries->pszLabels[m_pLastAccessedSeries->SeriesNumPoints * (cMaxPointLabelLen + 1)];
	if(pszLabel != NULL && pszLabel[0] != '\0')
		{
		strncpy(pszLabels,pszLabel,cMaxPointLabelLen);
		pszLabels[cMaxPointLabelLen] = '\0';
		}
	else
		sprintf(pszLabels,"P:%d",m_pLastAccessedSeries->SeriesNumPoints+1);
	}
m_pLastAccessedSeries->SeriesNumPoints += 1;
if(m_MaxSeriesNumPoints < m_pLastAccessedSeries->SeriesNumPoints)
	m_MaxSeriesNumPoints = m_pLastAccessedSeries->SeriesNumPoints;
return(m_pLastAccessedSeries->SeriesNumPoints);
}	

const int cSigPlotErrHandler = 21;  // thrown when PLPlot is unable to write SVG to output
int
CBKPLPlot::BKPLPlotErrHandler(const char *)
{
throw cSigPlotErrHandler;
return(1);
}
			


int
CBKPLPlot::PlotLineGraph(char *pszFile)				// write SVG plot to this file
{
int SeriesIdx;
tsPlotSeries *pSeries;

PLINT        opt_array[cMaxLineLegends];
const char   *text[cMaxLineLegends];
const char   *symbols[cMaxLineLegends];
PLINT        text_colors[cMaxLineLegends];
PLINT        line_colors[cMaxLineLegends];
PLINT        line_styles[cMaxLineLegends];
PLFLT        line_widths[cMaxLineLegends];
PLINT        symbol_numbers[cMaxLineLegends], symbol_colors[cMaxLineLegends];
PLFLT        symbol_scales[cMaxLineLegends];
PLFLT        legend_width, legend_height;
PLINT		 nlegend;

if(pszFile == NULL || pszFile[0] == '\0')
	return(eBSFerrParams);

if(m_NumSeries == 0 || m_MaxSeriesNumPoints == 0)
	return(eBSFerrNoEntries);

AcquireSerialise();
    // set plplot error handler
plsexit(BKPLPlotErrHandler);

try { 
	   // Initialize plplot
	plsfnam(pszFile);		// write SVG to this file
	plscolbg(255,255,255);  // white background, R G B
	plinit();
	pladv(0);				// Advance to subpage "page" or to the next page if "page" = 0

	plvstaextd(5.0,							// bottom margin in character heights
			   7.0,							// left margin in character heights
			   15.0,						// top margin in character heights
			   5.0);						// right margin in character heights

	// Set up world coordinates of the viewport boundaries (2d plots)
	plwind( m_xmin,							// xmin: The world x coordinate of the left-hand edge of the viewport  
			m_xmax,							// xmax: The world x coordinate of the right-hand edge of the viewport
			m_ymin,							// ymin: The world y coordinate of the bottom edge of the viewport
			m_ymax);						// ymax: The world y coordinate of the top edge of the viewport
  
	// This draws a box around the current viewport, complete with axes, ticks, numeric labels, and grids
	plcol0(BKPLbrown);
	plwidth(0.10);
	plbox(m_szxopt,						// xopt: draw bottom axis (a) and top (c) with grid lines (g), draw ticks (t) outside (i), numeric labels
			m_xtick,					// xtick: World coordinate interval between major ticks on the x axis. If it is set to zero, PLplot automatically generates a suitable tick interval
			m_nxsub,					// nxsub: Number of subintervals between major x axis ticks for minor ticks
			m_szyopt,					// yopt: draw left axis (a) and right (c), labels left (n), horizontal labels (v), draw ticks (t), and grid lines (g) 
			m_ytick,					// ytick: World coordinate interval between major ticks on the y axis. If it is set to zero, PLplot automatically generates a suitable tick interval
			m_nysub);					// nysub: Number of subintervals between major y axis ticks for minor ticks. If it is set to zero, PLplot automatically generates a suitable minor tick interval
    plcol0(BKPLbrown);

	plwidth(0.20);							// set pen width (0 is the minimum)

		// Create a labelled box to hold the plot.
	pllab( m_szXTitle[0] != '\0' ? m_szXTitle : "X", 
				m_szYTitle[0] != '\0' ? m_szYTitle :"Y", 
				m_szTitle[0] != '\0' ? m_szTitle : "Line Plot" );

		// Plot data added with AddXYPoint()
	pSeries = m_pSeries;
	nlegend = 0;
	for(SeriesIdx = 0; SeriesIdx < m_NumSeries && pSeries != NULL; SeriesIdx++, pSeries = pSeries->pNext)
		{
		if(pSeries->SeriesNumPoints == 0)
			continue;
		pllsty(pSeries->Style);
		plcol0(pSeries->Colour);
		plwidth(pSeries->Thickness);

		plline(pSeries->SeriesNumPoints,pSeries->pValues, &pSeries->pValues[m_MaxPerSeriesPoints]);

		if(m_NumSeries > 1)
			{
			if(SeriesIdx < cMaxLineLegends)
				{
	
				if(pSeries->Symbol == '\0')
					{
					opt_array[nlegend]   = PL_LEGEND_LINE;
					symbol_colors[nlegend]  = 0;
					symbol_scales[nlegend]  = 0;
					symbol_numbers[nlegend] = 0;
					symbols[nlegend]        = "";
					}
				else
					{
					opt_array[nlegend]   = PL_LEGEND_LINE | PL_LEGEND_SYMBOL;
					symbol_colors[nlegend]  = 3;
					symbol_scales[nlegend]  = 1.;
					symbol_numbers[nlegend] = 4;
					symbols[nlegend]        = (const char *)&pSeries->Symbol;
					}
				text_colors[nlegend] = 2;
				text[nlegend]        = pSeries->szTitle;
				symbols[nlegend] = "";
				line_colors[nlegend] = pSeries->Colour;
				line_styles[nlegend] = pSeries->Style;
				line_widths[nlegend] = pSeries->Thickness;
				nlegend += 1;
				}
			}
		}
	if(m_NumSeries > 1)
		{
		plscol0a( 15, 32, 32, 32, 0.70 );
		pllegend( &legend_width, &legend_height,
			PL_LEGEND_BACKGROUND | PL_LEGEND_BOUNDING_BOX, 0,
			0.0, 0.0, 0.1, 
			15,                                 // background colour from col0
			1, 1, 0, 0,
			nlegend, opt_array,
			1.0, 0.75, 1.25,					// text_offset, text_scale, text_spacing
			1.,									// text_justification
			text_colors, (const char **) text,
			NULL, NULL, NULL, NULL,
			line_colors, line_styles, line_widths,
			symbol_colors, symbol_scales, symbol_numbers, (const char **) symbols );
		}
	}

catch(int e)
	{
	if(e != cSigPlotErrHandler)
		throw e;
	plend();
	ReleaseSerialise();
	return(eBSFerrOpnFile);
	}

    // Close PLplot library
plend();
ReleaseSerialise();
return(eBSFSuccess);
}


int
CBKPLPlot::PlotBarChartGraph(char *pszFile,				// write SVG plot to this file
							int SeriesID)				// plot this series
{
int          i;
char         string[100];
PLFLT        *pValue;
tsPlotSeries *pSeries;

 static PLFLT pos[]   = { 0.0, 0.25, 0.5, 0.75, 1.0 };
 static PLFLT red[]   = { 0.0, 0.25, 0.5, 1.0, 1.0 };
 static PLFLT green[] = { 1.0, 0.5, 0.5, 0.5, 1.0 };
 static PLFLT blue[]  = { 1.0, 1.0, 0.5, 0.25, 0.0 };

// Initialize plplot
if(pszFile == NULL || pszFile[0] == '\0' || SeriesID == 0)
	return(eBSFerrParams);

if(m_NumSeries == 0 || m_MaxSeriesNumPoints == 0 || SeriesID > m_NumSeries)
	return(eBSFerrNoEntries);

pSeries = m_pSeries;
while(pSeries != NULL && pSeries->SeriesID != SeriesID)
	pSeries = pSeries->pNext;
if(pSeries == NULL)
	return(eBSFerrNoEntries);
pValue = pSeries->pValues;
	

AcquireSerialise();

    // set plplot error handler
plsexit(BKPLPlotErrHandler);

try { 
	   // Initialize plplot
	plsfnam(pszFile);						// write SVG to this file
	plscolbg(255,255,255);					// white background, R G B
	plinit();

    pladv( 0 );								// Advance to subpage "page" or to the next page if "page" = 0
 
	// Defines a "standard" viewport with seven character heights for the left margin and four character heights everywhere else.
	plvstaextd(5.0,							// bottom margin in character heights
			   7.0,							// left margin in character heights
			   15.0,						// top margin in character heights
			   5.0);						// right margin in character heights
 
	// Set up world coordinates of the viewport boundaries (2d plots)
	plwind( m_xmin,							// xmin: The world x coordinate of the left-hand edge of the viewport  
			m_xmax,							// xmax: The world x coordinate of the right-hand edge of the viewport
			m_ymin,							// ymin: The world y coordinate of the bottom edge of the viewport
			m_ymax);			
  
	// This draws a box around the current viewport, complete with axes, ticks, numeric labels, and grids
	plcol0(BKPLbrown);
	plwidth(0.10);
	plbox(m_szxopt,						// xopt: draw bottom axis (a) and top (c) with grid lines (g), draw ticks (t) outside (i), numeric labels
			m_xtick,					// xtick: World coordinate interval between major ticks on the x axis. If it is set to zero, PLplot automatically generates a suitable tick interval
			m_nxsub,					// nxsub: Number of subintervals between major x axis ticks for minor ticks
			m_szyopt,					// yopt: draw left axis (a) and right (c), labels left (n), horizontal labels (v), draw ticks (t), and grid lines (g) 
			m_ytick,					// ytick: World coordinate interval between major ticks on the y axis. If it is set to zero, PLplot automatically generates a suitable tick interval
			m_nysub);					// nysub: Number of subintervals between major y axis ticks for minor ticks. If it is set to zero, PLplot automatically generates a suitable minor tick interval
    plcol0(BKPLbrown);

	plwidth(0.20);							// set pen width (0 is the minimum)

		// Create a labelled box to hold the plot.
	pllab( m_szXTitle[0] != '\0' ? m_szXTitle : "X", 
				m_szYTitle[0] != '\0' ? m_szYTitle :"Y", 
				m_szTitle[0] != '\0' ? m_szTitle : "Bar Plot" );

	// Set color map 1 colors using a piece-wise linear relationship between intensity [0,1] (cmap 1 index) and position in HLS or RGB color space
    plscmap1l( 1,							// itype: 1: RGB, 0: HLS
				5,							// npts: number of control points
				pos,						// intensity: intensity index for each control point (between 0.0 and 1.0, in ascending order) 
				red,						// coord1: first coordinate (H or R) for each control point
				green,						// coord2: second coordinate (L or G) for each control point
				blue,						// coord3: third coordinate (S or B) for each control point
				NULL );						// alt_hue_path: alternative interpolation method flag for each control point. (alt_hue_path[i] refers to the interpolation interval between the i and i + 1 control points).

    for ( i = 0; i < max(1,pSeries->SeriesNumPoints); i++, pValue++ )
		{
        plcol1( 0.5 );					// Set color, map 1.  Argument is a float between 0. and 1		
//        plpsty( i % 8 );					// Set fill pattern, using one of the predefined patterns, 0..8
		if(*pValue)
			{
			plfbox( ( i + m_xmin), *pValue );   // draw bar
//			sprintf( string, "%.1f", *pValue );

			// Prints out "text" at world cooordinate within the viewport. The text may be at any angle "angle" relative to the horizontal
//			plschr(0,0.7);
//			plptex( (  m_xmin + i),		// x: x coordinate of reference point of string
//					(*pValue + 0.025 ),			// y: y coordinate of reference point of string
//					1.0,						// dx: Together with dy, this specifies the inclination of the string. The baseline of the string is parallel to a line joining (x, y) to (x+dx, y+dy)
//					0.0,						// dy: Together with dx, this specifies the inclination of the string
//					0.5,						// just: Specifies the position of the string relative to its reference point. If just=0., the reference point is at the left and if just=1., it is at the right of the string. Other values of just give intermediate justifications
//					string );					// text: The string to be written out
//			plschr(0,1.0);
			}

		if(pSeries->SeriesNumPoints < 20 || !(i % 2))
			{
			if(pSeries->bLabels && pSeries->pszLabels)
				strcpy(string,&pSeries->pszLabels[i * (cMaxPointLabelLen + 1)]);
			else
				sprintf( string, "%d", (int)(i + m_xmin));

			// Write text relative to viewport boundaries
			plcol0(BKPLbrown);
			plmtex( "b",						//	side: Bottom of viewport (b), add 'v' for text written at right angles to edge
					 1.1,						//  disp: Position of the reference point of string, measured outwards from the specified viewport edge in units of the current character height. Use negative disp to write within the viewport
					(1.0/(2.0 * pSeries->SeriesNumPoints)) + (i/(double)pSeries->SeriesNumPoints),		//  pos: Position of the reference point of string along the specified edge, expressed as a fraction of the length of the edge
					0.5,						//  just: Specifies the position of the string relative to its reference point. If just=0., the reference point is at the left and if just=1., it is at the right of the string. Other values of just give intermediate justifications
					string );					//  text: The string to be written out

			}
		}
    }
catch(int e)
	{
	if(e != cSigPlotErrHandler)
		throw e;
	plend();
	ReleaseSerialise();
	return(eBSFerrOpnFile);
	}

    // Close PLplot library
plend();
ReleaseSerialise();
return(eBSFSuccess);
}

void
CBKPLPlot::plfbox( PLFLT x0, PLFLT y0 )
{
    PLFLT x[4], y[4];

    x[0] = x0;
    y[0] = 0.0;
    x[1] = x0;
    y[1] = y0;
    x[2] = x0 + 1.0;
    y[2] = y0;
    x[3] = x0 + 1.0;
    y[3] = 0.0;
    plfill( 4, x, y );			// Pattern fills the polygon bounded by the input points
    plcol0( 1 );
    pllsty( 1 );
    plline( 4, x, y );
}

PLINT
CBKPLPlot::zdefined( PLFLT x, PLFLT y )
{
    PLFLT z = sqrt( x * x + y * y );

    return z < 0.4 || z > 0.6;
}



//--------------------------------------------------------------------------
// Phred score distribution shade plot
//--------------------------------------------------------------------------

int
CBKPLPlot::PlotPhredScores(char *pszFile)	    // write SVG Phred score distributions graph plot to this file)
{
tsPlotSeries *pSeries;
PLFLT *pValue;
int ns;        // number of shade levels
int nx;        // number of data points in x
int ny;        // number of data points in y

    int        i, j;
    PLFLT      **z, zmin, zmax;
    PLFLT      *shedge;
    PLFLT      fill_width = 2., cont_width = 0.;
    PLFLT      colorbar_width, colorbar_height;
    PLINT      cont_color = 0;
#define NUM_AXES    1
    PLINT      n_axis_opts  = NUM_AXES;
    const char *axis_opts[] = { "bcvtm", };
    PLINT      num_values[NUM_AXES];
    PLFLT      *values[NUM_AXES];
    PLFLT      axis_ticks[NUM_AXES] = { 0.0, };
    PLINT      axis_subticks[NUM_AXES] = {  0,  };
#define NUM_LABELS    1
    PLINT      n_labels     = NUM_LABELS;
    PLINT      label_opts[] = { PL_COLORBAR_LABEL_BOTTOM, };
    const char *labels[] = { "Proportion", };

// Initialize plplot
if(pszFile == NULL || pszFile[0] == '\0')
	return(eBSFerrParams);

if(m_NumSeries != 42 || m_MaxSeriesNumPoints == 0)
	return(eBSFerrNoEntries);

pSeries = m_pSeries;
if(pSeries == NULL)
	return(eBSFerrNoEntries);
nx = pSeries->SeriesNumPoints;
ny = 42;
ns = ny/2;
	
AcquireSerialise();
    // set plplot error handler
plsexit(BKPLPlotErrHandler);

try {
	plsfnam(pszFile);						// write SVG to this file
	// Load colour palettes
    plspal0( "cmap0_black_on_white.pal" );
    plspal1( "cmap1_blue_yellow.pal", 1 );
	// Reduce colors in cmap 0 so that cmap 1 is useful on a 16-color display
    plscmap0n( 3 );
	
    plinit();				// Initialize plplot
	pladv(0);				// Advance to subpage "page" or to the next page if "page" = 0

	// Allocate data structures
    shedge = (PLFLT *) calloc( (size_t) ( ns + 1 ), sizeof ( PLFLT ) );
    plAlloc2dGrid( &z, nx, ny );

	// Set up data array
	pSeries = m_pSeries;
	for ( j = 0; j < ny; j++,pSeries = pSeries->pNext)
		{
		pValue = pSeries->pValues;
		for ( i = 0; i < nx; i++ )
			{
			z[i][j] = *pValue++;
			}
		}

    f2mnmx( z, nx, ny, &zmin, &zmax );

    for ( i = 0; i < ns + 1; i++ )
        shedge[i] = zmin + ( zmax - zmin ) * (PLFLT) i / (PLFLT) ns;

    plvpor( 0.1, 0.9, 0.1, 0.9 );
    plwind( 0.0, nx, 0.0, ny);

//	plcol0(BKPLbrown);
	plwidth(0.10);
	plbox(m_szxopt,						// xopt: draw bottom axis (a) and top (c) with grid lines (g), draw ticks (t) outside (i), numeric labels
			m_xtick,					// xtick: World coordinate interval between major ticks on the x axis. If it is set to zero, PLplot automatically generates a suitable tick interval
			m_nxsub,					// nxsub: Number of subintervals between major x axis ticks for minor ticks
			m_szyopt,					// yopt: draw left axis (a) and right (c), labels left (n), horizontal labels (v), draw ticks (t), and grid lines (g) 
			m_ytick,					// ytick: World coordinate interval between major ticks on the y axis. If it is set to zero, PLplot automatically generates a suitable tick interval
			m_nysub);					// nysub: Number of subintervals between major y axis ticks for minor ticks. If it is set to zero, PLplot automatically generates a suitable minor tick interval

	plwidth(0.20);							// set pen width (0 is the minimum)
		// Create a labelled box to hold the plot.
	pllab( m_szXTitle[0] != '\0' ? m_szXTitle : "X", 
				m_szYTitle[0] != '\0' ? m_szYTitle :"Y", 
				m_szTitle[0] != '\0' ? m_szTitle : "Colour Map Plot" );
    plpsty( 0 );

    plshades( (const PLFLT * const *) z,   // ** pointer to array to be plotted. The array must have been declared as PLFLT a[nx][ny] 
		nx, ny,								// dimension x, y	 
		NULL,                               // user function for specifying regions to be excluded 
		0.0,                                // grid coordinate x min 
		(double)nx,							// grid coordinate x max
		0.0,								// grid coordinate y min
		(double)ny,							// grid coordinate y max
        shedge,								// Pointer to array containing the data levels corresponding to the edges of each shaded region that will be plotted by this function. To work properly the levels should be monotonic.
		ns + 1,								// Number of shades plus 1 (i.e., the number of shade edge values in shedge.
		fill_width,							// Defines line width used by the fill pattern
        cont_color,                         // Defines pen color used for contours defining edges of shaded regions. The pen color is only temporary set for the contour drawing. Set this value to zero or less if no shade edge contours are wanted 
		cont_width,							// Defines line width used for contours defining edges of shaded regions. This value may not be honored by all drivers. The pen width is only temporary set for the contour drawing. Set this value to zero or less if no shade edge contours are wanted.
        plfill,								// Routine used to fill the region
		1,									// Set rectangular to true if rectangles map to rectangles after coordinate transformation with pltrl. Otherwise, set rectangular to false.
		NULL,								// Pointer to function that defines transformation between indices in array z and the world coordinates (C only).
		NULL );								// Extra parameter to help pass information to pltr0, pltr1, pltr2, or whatever routine that is externally supplied.

    // Smaller text
    plschr( 0.0, 0.75 );
    // Small ticks on the vertical axis
    plsmaj( 0.0, 0.5 );
    plsmin( 0.0, 0.5 );

    num_values[0] = ns + 1;
    values[0]     = shedge;
    plcolorbar( &colorbar_width, &colorbar_height,
        PL_COLORBAR_SHADE | PL_COLORBAR_SHADE_LABEL, 0,
        0.005, 0.0, 0.0375, 0.875, 0, 1, 1, 0.0, 0.0,
        cont_color, cont_width,
        n_labels, label_opts, labels,
        n_axis_opts, axis_opts,
        axis_ticks, axis_subticks,
        num_values, (const PLFLT * const *) values );

    // Reset text and tick sizes
    plschr( 0.0, 1.0 );
    plsmaj( 0.0, 1.0 );
    plsmin( 0.0, 1.0 );

	// Clean up
    free( (void *) shedge );
    plFree2dGrid( z, nx, ny );
	}

catch(int e)
	{
	if(e != cSigPlotErrHandler)
		throw e;
	plend();
	ReleaseSerialise();
	return(eBSFerrOpnFile);
	}

    // Close PLplot library
plend();
ReleaseSerialise();
return(eBSFSuccess);
}

//--------------------------------------------------------------------------
// f2mnmx
//
// Returns min & max of input 2d array.
//--------------------------------------------------------------------------

void
CBKPLPlot::f2mnmx( PLFLT **f, PLINT nnx, PLINT nny, PLFLT *fnmin, PLFLT *fnmax ) // Returns min & max of input 2d array.
{
    int i, j;

    *fnmax = f[0][0];
    *fnmin = *fnmax;

    for ( i = 0; i < nnx; i++ )
    {
        for ( j = 0; j < nny; j++ )
        {
            *fnmax = MAX( *fnmax, f[i][j] );
            *fnmin = MIN( *fnmin, f[i][j] );
        }
    }
}


