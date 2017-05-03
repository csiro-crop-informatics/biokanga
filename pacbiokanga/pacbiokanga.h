/*
 * CSIRO Open Source Software License Agreement (GPLv3)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENSE for the complete license information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */

#pragma once

#define USEPBfilter 1

#include "SQLiteSummaries.h"

typedef struct TAG_sSubProcess {
	const char *pszName;	// subprocess name
	const char *pszBriefDescr; // subprocess brief description
	const char *pszFullDescr; // subprocess full description
	int (* SubFunct)(int argc, char* argv[]);
} tsSubProcess;

extern const char *cpszProgVer;				// will be incremented with each release
extern CStopWatch gStopWatch;				// time keeper
extern CDiagnostics gDiagnostics;			// for writing diagnostics messages to log file
extern CSQLiteSummaries gSQLiteSummaries;	// for writing processing result summaries to SQLite database
extern int	gExperimentID;					// SQLite experiment identifier
extern int gProcessID;						// SQLite processor identifier
extern int	gProcessingID;					// SQLite processing identifier

extern char gszProcName[_MAX_FNAME];		// process name
extern tsSubProcess *gpszSubProcess;		// selected subprocess