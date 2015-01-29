// EmbeddedPalettes.h
// Contains the colour information for the standard palettes distributed with PLPlot as initialised structures
// Saves having to install the original palette files and configuring the plotting application for location of these palette files
// Please contact stuart.stephen@csiro.au for details

#if !defined ( _EMBEDDEDPALETTES_H )
#define _EMBEDDEDPALETTES_H

const int cMaxEmbeddedPaletteColors = 16;    // max supported colours in embedded palettes
typedef struct {
			unsigned int Red;
			unsigned int Green;
			unsigned int Blue;
			PLFLT Alpha;
	} tsCol0ColorInfo;

typedef struct {
    PLFLT r;
	PLFLT g;
	PLFLT b;
	PLFLT a; 
	PLFLT pos; 
	PLBOOL alt_hue_path;
	} tsCol1ColorInfo;


typedef struct TAG_sEmbeddedCol0Palette {
	const char *pszPaletteName;	// same as the palette file from which this palette was extracted
	int NumColors;			// number of colours in this palette (max 16)
	tsCol0ColorInfo ColorInfo[cMaxEmbeddedPaletteColors];
	} tsEmbeddedCol0Palette;

typedef struct TAG_sEmbeddedCol1Palette {
	const char *pszPaletteName;	// same as the palette file from which this palette was extracted
	int NumColors;			// number of colours in this palette (max 16)
	int Format;             // 0 original, 1 new V2 hls , 2 new V2 rgb
	tsCol1ColorInfo ColorInfo[cMaxEmbeddedPaletteColors];
	} tsEmbeddedCol1Palette;


tsEmbeddedCol1Palette EmbeddedCol1Palettes[] = {
	{ "cmap1_default.pal",
	6, 0,
	{{(double)0x55/255.0,(double)0x00/255.0,(double)0xff/255.0, 1.0, 00.0 * 0.01, 0},
	 {(double)0x11/255.0,(double)0x00/255.0,(double)0x33/255.0, 1.0, 44.0 * 0.01, 0},
	 {(double)0x01/255.0,(double)0x00/255.0,(double)0x05/255.0, 1.0, 50.0 * 0.01, 0},
	 {(double)0x05/255.0,(double)0x00/255.0,(double)0x00/255.0, 1.0, 50.0 * 0.01, 0},
	 {(double)0x33/255.0,(double)0x00/255.0,(double)0x00/255.0, 1.0, 56.0 * 0.01, 0},
	 {(double)0xff/255.0,(double)0x00/255.0,(double)0x00/255.0, 1.0, 100.0 * 0.01, 0}}},


	{ "cmap1_blue_red.pal",  // v2 hsl
    2,1,
	{{240.0, 0.5, 1.0, 1.0, 0.0, 0},
	 {  0.0, 0.5, 1.0, 1.0, 1.0, 0}}},

	{ "cmap1_highfreq.pal",  // v2 hsl
    2,1,
	{{240.0, 0.3, 0.5, 1.0, 0.0, 0},
	 {240.0, 1.0, 0.5, 1.0, 1.0, 0}}},

	{ "cmap1_lowfreq.pal",  // v2 hsl
	4,1,
	{{ 240.0, 0.5, 1.0, 1.0, 0.000, 0},
	 { 240.0, 0.5, 0.0, 1.0, 0.499, 0},
	 {  60.0, 0.5, 0.0, 1.0, 0.501, 0},
	 {  60.0, 0.5, 1.0, 1.0, 1.000, 0}}},


	{ "cmap1_grey.pal",
	2,2,
	{{0.0,0.0,0.0,1.0,0.0,0},
	 {1.0,1.0,1.0,1.0,1.0,0}}},

	{ "cmap1_blue_yellow.pal",
	7,2,
	{{   1.000, 1.0,  1.00, 1.0,0.000, 0}, 
     {   0.001, 0.998,1.00, 1.0,0.001, 0},
     {   0.250, 0.500,1.00, 1.0,0.250, 0},
     {   0.500, 0.500,0.50, 1.0,0.500, 0},
     {   1.000, 0.500,0.25, 1.0,0.750, 0},
     {   1.000, 0.998,0.001,1.0,0.999, 0},
     {   0.000, 0.000,0.000,1.0,1.000, 0}}},

	{ "cmap1_radar.pal",					
	16,2,
	{{0.00,0.00,0.00,0.0,0.00,1},
	 {0.00,0.93,0.93,1.0,0.07,0},
	 {0.00,0.63,0.96,1.0,0.13,0},
	 {0.00,0.00,0.93,1.0,0.20,0},
	 {0.00,1.00,0.00,1.0,0.27,0},
	 {0.00,0.78,0.00,1.0,0.33,0},
	 {0.00,0.56,0.00,1.0,0.40,0},
	 {1.00,1.00,0.00,1.0,0.47,0},
	 {0.91,0.75,0.00,1.0,0.53,0},
	 {1.00,0.56,0.00,1.0,0.60,0},
	 {1.00,0.00,0.00,1.0,0.67,0},
	 {0.84,0.00,0.00,1.0,0.73,0},
	 {0.75,0.00,0.00,1.0,0.80,1},
	 {1.00,0.00,1.00,1.0,0.87,0},
	 {0.60,0.33,0.79,1.0,0.93,1},
	 {1.00,1.00,1.00,1.0,1.00,0}}}
	};
const int cNumEmbeddedCol1Pallets = sizeof(EmbeddedCol1Palettes)/sizeof(tsEmbeddedCol1Palette);


tsEmbeddedCol0Palette EmbeddedCol0Palettes[] = {
	{ "cmap0_default.pal",
	16,
	{{0x0000,0x0000,0x0000,1.0},
	{0x00ff,0x0000,0x0000,1.0},
	{0x00ff,0x00ff,0x0000,1.0},
	{0x0000,0x00ff,0x0000,1.0},
	{0x007f,0x00ff,0x00d4,1.0},
	{0x00ff,0x00c0,0x00cb,1.0},
	{0x00f5,0x00de,0x00b3,1.0},
	{0x00be,0x00be,0x00be,1.0},
	{0x00a5,0x002a,0x002a,1.0},
	{0x0000,0x0000,0x00ff,1.0},
	{0x008a,0x002b,0x00e2,1.0},
	{0x0000,0x00ff,0x00ff,1.0},
	{0x0040,0x00e0,0x00d0,1.0},
	{0x00ff,0x0000,0x00ff,1.0},
	{0x00fa,0x0080,0x0072,1.0},
	{0x00ff,0x00ff,0x00ff,1.0}}
	},

	{ "cmap0_alternate.pal",
	16,
	{{0x00ff,0x00ff,0x00ff,1.0},
	{0x0000,0x0000,0x0000,1.0},
	{0x0000,0x0000,0x00ff,1.0},
	{0x00ff,0x0000,0x0000,1.0},
	{0x00a5,0x002a,0x002a,1.0},
	{0x00fa,0x0080,0x0072,1.0},
	{0x00ff,0x00c0,0x00cb,1.0},
	{0x007f,0x00ff,0x00d4,1.0},
	{0x00f5,0x00de,0x00b3,1.0},
	{0x0040,0x00e0,0x00d0,1.0},
	{0x00be,0x00be,0x00be,1.0},
	{0x0000,0x00ff,0x00ff,1.0},
	{0x0000,0x00ff,0x0000,1.0},
	{0x00ff,0x00ff,0x0000,1.0},
	{0x00ff,0x0000,0x00ff,1.0},
	{0x008a,0x002b,0x00e2,1.0}}
	},

	{ "cmap0_black_on_white.pal",
	16,
	{{0x00ff,0x00ff,0x00ff,1.0},
	{0x0000,0x0000,0x0000,1.0},
	{0x0000,0x0000,0x0000,1.0},
	{0x0000,0x0000,0x0000,1.0},
	{0x0000,0x0000,0x0000,1.0},
	{0x0000,0x0000,0x0000,1.0},
	{0x0000,0x0000,0x0000,1.0},
	{0x0000,0x0000,0x0000,1.0},
	{0x0000,0x0000,0x0000,1.0},
	{0x0000,0x0000,0x0000,1.0},
	{0x0000,0x0000,0x0000,1.0},
	{0x0000,0x0000,0x0000,1.0},
	{0x0000,0x0000,0x0000,1.0},
	{0x0000,0x0000,0x0000,1.0},
	{0x0000,0x0000,0x0000,1.0},
	{0x0000,0x0000,0x0000,1.0}}},

	{ "cmap0_white_bg.pal",
	16,
	{{0x00ff,0x00ff,0x00ff,1.0},
	{0x00ff,0x0000,0x0000,1.0},
	{0x00ff,0x00ff,0x0000,1.0},
	{0x0000,0x00ff,0x0000,1.0},
	{0x007f,0x00ff,0x00d4,1.0},
	{0x00ff,0x00c0,0x00cb,1.0},
	{0x00f5,0x00de,0x00b3,1.0},
	{0x00be,0x00be,0x00be,1.0},
	{0x00a5,0x002a,0x002a,1.0},
	{0x0000,0x0000,0x00ff,1.0},
	{0x008a,0x002b,0x00e2,1.0},
	{0x0000,0x00ff,0x00ff,1.0},
	{0x0040,0x00e0,0x00d0,1.0},
	{0x00ff,0x0000,0x00ff,1.0},
	{0x00fa,0x0080,0x0072,1.0},
	{0x0000,0x0000,0x0000,1.0}}}
	};
const int cNumEmbeddedCol0Pallets = sizeof(EmbeddedCol0Palettes)/sizeof(tsEmbeddedCol0Palette);


#endif   // _EMBEDDEDPALETTES_H
