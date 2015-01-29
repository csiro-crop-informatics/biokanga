#pragma once

const int cFiltRefIDsAllocChunk = 10000;	// allocate for RefIDs in this many instances

class CFilterRefIDs : CCSVFile
{
	int m_NumFilterRefIDs;			// number of RefIDs in m_pFilterRefIDs
	int m_AllocdFilterRefIDs;		// allocd RefIDs
	int *m_pFilterRefIDs;			// array of RefIDs
	static int SortRefIDs( const void *arg1, const void *arg2);
public:
	CFilterRefIDs(void);
	~CFilterRefIDs(void);
	int Open(char *pszFile);		// open and load RefIDs from file, expected to be in first field
	void Reset(void);				// resets state back to that immediately following class instantiation
	bool Locate(int RefID);			// returns true if RefID is in RefIDs loaded
};
