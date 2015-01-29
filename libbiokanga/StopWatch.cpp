// Copyright 2013 CSIRO  ( http://www.csiro.au/ ) 
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

#include "stdafx.h"
#ifdef _WIN32
#include "./commhdrs.h"
#else
#include "./commhdrs.h"
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CStopWatch::CStopWatch()
{
m_StartCnt = 0;
m_Elapsed = 0;
m_Started = 0;
#ifdef _WIN32
if(!QueryPerformanceFrequency((LARGE_INTEGER *)&m_Freq))
	m_Freq = 0;
#else
m_Freq = sysconf(_SC_CLK_TCK);
#endif
}

CStopWatch::~CStopWatch()
{

}

void 
CStopWatch::Start(void)
{
#ifdef _WIN32
if(!m_StartCnt++)	// if first time then get current time
	QueryPerformanceCounter((LARGE_INTEGER *)&m_Started);
#else
struct tms Times;
if(!m_StartCnt++)       // if first time then get current time
	m_Started = times(&Times);
#endif
}

void
CStopWatch::Stop(void)
{
#ifdef _WIN32
INT64 Now;
if(m_StartCnt && m_StartCnt-- == 1)	// if startcnt is decremented down to 0 then accumulate time difference
	{
	QueryPerformanceCounter((LARGE_INTEGER *)&Now);
	m_Elapsed += Now - m_Started;
	}
#else
struct tms Times;
clock_t Now;
if(m_StartCnt && m_StartCnt-- == 1)     // if startcnt is decremented down to 0 then accumulate time difference
        {
        Now = (INT64)times(&Times);
        m_Elapsed += Now - m_Started;
        }
#endif
}

void
CStopWatch::Reset(void)
{
m_StartCnt = 0;
m_Elapsed = 0;
m_Started = 0;
}

char *
CStopWatch::Read(void)
{
unsigned long Secs;
unsigned long USecs;
Secs = ReadUSecs(&USecs);
sprintf(m_szDuration,"%6.2d:%2.2d:%2.2d.%3.3d seconds\n",(int)(Secs/3600),(int)((Secs % 3600)/60),(int)(Secs%60),(int)(USecs/1000));
return(m_szDuration);
}

unsigned long		// returned number of seconds
CStopWatch::ReadUSecs(unsigned long *pMicroSecs) // optional microsecs (secs.microsecs == elapsed time)
{
#ifdef _WIN32
INT64 Now;
INT64 Elapsed;
unsigned long Secs;
Elapsed = m_Elapsed;
if(m_StartCnt)			// if still timing...
	{
	QueryPerformanceCounter((LARGE_INTEGER *)&Now);
	Elapsed += Now - m_Started;
	}
Secs = (unsigned long)(Elapsed / m_Freq);	// unlikely that any process would stay up for ulong secs!
if(pMicroSecs != NULL)
	*pMicroSecs = (unsigned long)(((Elapsed - ((INT64)Secs * m_Freq)) * (INT64)1000000)/m_Freq);
return(Secs);
#else
struct tms Times;
INT64 Now;
INT64 Elapsed;
unsigned long Secs;
Elapsed = m_Elapsed;
if(m_StartCnt)                  // if still timing...
        {
        Now = (INT64)times(&Times);
        Elapsed += Now - m_Started;
        }
Secs = (unsigned long)(Elapsed / m_Freq);       // unlikely that any process would stay up for ulong secs!
if(pMicroSecs != NULL)
        *pMicroSecs = (unsigned long)(((Elapsed - ((INT64)Secs * m_Freq)) * (INT64)1000000)/m_Freq);
return(Secs);
#endif
}
