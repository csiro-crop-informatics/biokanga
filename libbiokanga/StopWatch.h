#pragma once
#include "./commdefs.h"
class CStopWatch  
{
	char m_szDuration[30];			// used to return stopwatch as a formatted string - "h:mm:ss.ms seconds\n"
	INT64 m_Freq;					// stopwatch cnt freq in ticks/sec
	INT64 m_Elapsed;				// total elapsed ticks since stopwatch intially started
	INT64 m_Started;				// when stopwatch was initially started
	int m_StartCnt;					// keeps count of starts w/o corresponding stop 
public:
	CStopWatch();
	virtual ~CStopWatch();
	void Reset(void);				// resets stopwatch
	void Start(void);				// starts stopwatch
	void Stop(void);				// stops stopwatch
	char *Read(void);				// returns elapsed time as string - "h:mm:ss.ms seconds\n"
	unsigned long ReadUSecs(unsigned long *pMicroSecs=NULL);  // read stopwatch, returns seconds with optional Microsecs 
};
