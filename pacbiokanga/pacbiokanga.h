// Copyright 2014 CSIRO  ( http://www.csiro.au/ ) 
//   Licensed under the Apache License, Version 2.0 (the "License");
//   you may not use this file except in compliance with the License.
//   You may obtain a copy of the License at
//       http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing, software
//   distributed under the License is distributed on an "AS IS" BASIS,
//   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//   See the License for the specific language governing permissions and
//   limitations under the License
//   Please contact stuart.stephen@csiro.au for support or 
//   to submit modifications to this source

#pragma once

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