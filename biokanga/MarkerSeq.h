/*
 * CSIRO Open Source Software License Agreement (GPLv3)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENSE for the complete license information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */

#pragma once

const int cMarkerBuffAlloc = 0x0fffff;	// allocation size for m_pMarkerBuff
const int cMaxMarkerSeqLen = 1000;		// max allowed generated marker sequence length

class CMarkerSeq
{
	int m_PMode;				// currently default processing only is supported
	char m_szInLociFile[_MAX_PATH];		// Loci file specifying the SNP or region loci in assembly for which marker sequences are to be generated (CSV, BED or SAM)
	char m_szInFastaFile[_MAX_PATH];		// multifasta assembly file from which marker sequences are to be generated
	char m_szOutFile[_MAX_PATH];			// marker sequences written to this file
	int m_hOutMarkerFile;					// marker output file handle

	int m_Ftype;				// input element file format: 0 - auto, 1 - CSV, 2 - BED, 3 - SAM
	teCSVFormat m_CSVFormat;	// expected input CSV loci file format
	int m_MaxNs;				// filter out input sequences having higher than this number of indeterminate bases per 100bp (default is 0, range 0..5)
	int m_Extend5;				// extend markers 5' bases from centroid of loci (default is 0, range 0..250)
	int m_Extend3;				// extend markers 3' bases from centroid of loci (default is 0, range 5..250)
	int m_MinSeqLen;			// filter out marker sequences which are less than this length (default is 10bp, range 10..1000)
	int m_MaxSeqLen;			// filter out marker sequences which are longer than this length (default is 1000bp, minseqlen..1000)

	bool m_bNoOverlaps;			// filter out marker sequences which overlap with other marker sequences

	int m_NumEls;				// number of element loci parsed and accepted from m_szInLociFile 

	CHyperEls *m_pHypers;		// used to contain accepted loci
	CFasta *m_pFasta;			// multifasta file containing assembly sequences from which markers are to be processed

	UINT8 *m_pSeq;				// allocated to hold sequences read from input multifasta file
	size_t m_AllocdSeqMem;		// this size allocation to m_pSeq;

	UINT8 *m_pMarkerBuff;		// allocated to buffer marker sequences to be written out to file
	size_t m_AllocdMarkerBuff;	// this size allocation to m_pSeq;
	UINT32 MarkerBuffOfs;		// offset in m_pMarkerBuff to next write 

public:
	CMarkerSeq(void);
	~CMarkerSeq(void);
	int Reset(void);
	int
		ProcessMarkerSeqs(int PMode,			// currently default processing only is supported
					int Ftype,					// input element file format: 0 - auto, 1 - CSV, 2 - BED, 3 - SAM
					teCSVFormat CSVFormat,		// expected input CSV loci file format
					int MaxNs,					// filter out marker sequences having higher than this number of indeterminate bases (default is 0, range 0..5)
					int Extend5,				// extend markers 5' bases from centroid of loci (default is 0, range 0..250)
					int Extend3,				// extend markers 3' bases from centroid of loci (default is 0, range 0..250)
					int MinSeqLen,				// filter out marker sequences which are less than this length (default is 10bp, range 10..1000)
					int MaxSeqLen,				// filter out marker sequences which are longer than this length (default is 1000bp, range minseqlen..1000)
					bool bNoOverlaps,			// filter out marker sequences which overlap with other marker sequences
					char *pszInLociFile,		// Loci file specifying the SNP or region loci in assembly from which marker sequences are to be generated (CSV, BED or SAM)
					char *pszInFastaFile,		// multifasta assembly file from which marker sequences are to be generated
					char *pszOutFile);			// marker sequences written to this file
};

