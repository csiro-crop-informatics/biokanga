/*
 * CSIRO Open Source Software License Agreement (GPLv3)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENSE for the complete license information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */

#pragma once

const int cMaxTextBuff = 4000;			// buffer upto this length results values text
const int cMaxNameLen = 50;				// maximum experiment/process/groupas name length
const int cMaxDescrText = 1000;			// max length descriptive text
const int cMaxSQLStatement = 500;		// allow for this length SQL statements

typedef enum eSQLliteSummParamTypes {
	 ePTBool = 0,
	 ePTInt32,
	 ePTUint32,
	 ePTInt64,
	 ePTUint64,
	 ePTDouble,
	 ePTText
} teSQLliteSummParamTypes;

typedef struct TAG_sSummStmsSQL {
	char *pTblName;					// table name
	char *pszCreateTbl;				// SQL statement used to create the table
	char *pszInsert;				// SQL statement used to insert row into table
	sqlite3_stmt *pPrepInsert;		// prepared insert row statement; NULL if not prepared
	char *pszOpenCreateIndexes;		// SQL statement used to create indexes on this table immediately after opening this database, essentially these are used whilst populating
	char *pszCreateIndexes;			// SQL statement used to create indexes on this table before closing the database
	char *pszDropIndexes;			// SQL statement used to drop indexes on this table
} tsSummStmSQL;

class CSQLiteSummaries
{
	bool m_bInitialised;				// set true after successful initialisation during class construction
	sqlite3 *m_pDB;						// pts to instance of SQLite
	static tsSummStmSQL m_StmSQL[7];	// SQLite table and index statements

	
	char *RemoveQuotes(char *pszRawText);	// remove quotes which may be a risk for aql injection attacks...

	int					// length - strlen() - of timestamp string 
		GetTimeStamp(char *pszTimeStamp);	// copy timestamp into this

	int Init(void);							// initialisation performed during class construction

	sqlite3 *OpenDatabase(char *pszDatabase,		// database to open, if not already existing then will be created
							bool bReplace = false);			// if database already exists then replace
	int CloseDatabase(void);

	static int ExecCallbackID(void *pCallP1, // callback function processing identifier (4th arg to sqlite3_exec())
					int NumCols,			// number of result columns 
					char **ppColValues,		// array of ptrs to column values 
					char **ppColName);		// array of ptrs to column names

	int							// eBSFSuccess if pszValueText contains textified value
		ValueToText(teSQLliteSummParamTypes ParamType,	// result value type
			int ValueSize,			// result value is of this byte length
			void *pValue,			// result value
			int MaxTextLen,			// truncate returned value text at this many chars - 1 (allows for terminating '\0') NOTE: error if insufficent length to hold max value for type
			char *pszValueText);		// user allocated to hold returned value as text
#ifdef _WIN32
	CRITICAL_SECTION m_hSCritSect;	// used to serialise access to SQLite functionality
#else
	pthread_spinlock_t m_hSpinLock;	// used to serialise access to SQLite functionality
#endif

	inline void SQLiteSerialise(void);	// serialise access to SQLite functionality
	inline void SQLiteRelease(void);

public:
	CSQLiteSummaries(void);
	~CSQLiteSummaries(void);

	int													// returned experiment identifier
			StartExperiment(char *pszDatabase,			// summary results to this SQLite database
							bool bReplace,				// if false then append to existing, if true then replace any existing database
							bool bContinue,				// if true then treat as if continuing pre-existing experiment
							char *pszExprimentName,		// experiment name
							char *pszExperimentTitle,	// experiment title
							char *pszExperimentDescr);  // describes experiment

	int					// returned process identifier 
			AddProcess(char *pszProcessName,			// process name
							char *pszProcessTitle,		// process title
							char *pszProcessDescr);		// describes process


	int													// uiniqely identifies this starting experiment process instance
		    StartProcessing(int ExperimentID,			// identifier returned by StartExperiment()
						 int ProcessID,					// identifier as returned by AddProcess()
						char *pszProcessVersion);		// process version

	int													// uiniqely identifies this added log text
			AddLog(int ProcessingID,					// identifier returned by StartProcessing()
						 const char *pszFormat,...);	// printf style format
	int
			AddParameter(int ProcessingID,		// identifier returned by StartProcessing()
						 teSQLliteSummParamTypes ParamType,	// parameter type
						 int ValueSize,					// parameter value is of this byte length
						 const char *pszParamName,		// parameter name
						 void *pParmValue);				// parameter value

	int
			AddResult(int ProcessingID,		// identifier returned by StartProcessing()
						const char *GroupAs,			// result is part of this grouping
						 teSQLliteSummParamTypes ParamType,	// result value type
						 int ValueSize,					// result value is of this byte length
						 const char *pszResultName,		// result name
						 void *pResultValue);			// result value

	int
			AddResultXY(int ProcessingID,			// identifier returned by StartProcessing()
						 const char *GroupAs,			// result is part of this grouping
						 teSQLliteSummParamTypes ParamXType,	// result value X type
						 int ValueXSize,				// result value is of this byte length
						 const char *pszResultXName,	// result X name
						 void *pResultXValue,			// result X value
 						 teSQLliteSummParamTypes ParamYType,	// result value Y type
						 int ValueYSize,				// result value is of this byte length
						 const char *pszResultYName,	// result Y name
						 void *pResultYValue);			// result Y value

	int					// returned process identifier
			EndProcessing(int ProcessingID,					// identifier returned by StartProcessing()
						  int ResultCode);					// processing result code

	int
			EndExperiment(int ExperimentID);			// identifier returned by StartExperiment()
};

