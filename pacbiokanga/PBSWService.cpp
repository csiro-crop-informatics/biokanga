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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#if _WIN32
#include <process.h>
#include "../libbiokanga/commhdrs.h"
#include <winsock2.h>

#else
#include <sys/mman.h>
#include <pthread.h>
#include <sys/ioctl.h>
#include <sys/socket.h>
#include <netdb.h>
typedef struct sockaddr_storage SOCKADDR_STORAGE;
#include "../libbiokanga/commhdrs.h"
#endif

#include "pacbiokanga.h"
#include "SSW.h"
#include "BKScommon.h"
#include "BKSProvider.h"
#include "PBSWService.h"

#ifdef _WIN32
int ProcSWService(int argc, char* argv[])
{
	// determine my process name
_splitpath(argv[0], NULL, NULL, gszProcName, NULL);
#else
int
ProcSWService(int argc, char** argv)
{
	// determine my process name
	CUtility::splitpath((char *)argv[0], NULL, gszProcName);
#endif

	int iFileLogLevel;			// level of file diagnostics
	int iScreenLogLevel;		// level of file diagnostics
	char szLogFile[_MAX_PATH];	// write diagnostics to this file

	char szHostName[cMaxHostNameLen];			// connecting on this host name or IPv4/IPv5 address to service providers
	char szServiceName[cMaxServiceNameLen];		// connect on this service name or port to service providers
	int MaxServInsts;							// max number of service instances supported
	int MaxServMemGB;							// max allocatable memory available for all service instances 
	int MaxConnWait;							// wait for at most this many minutes for connection 
	int PMode;

	// command line args
	struct arg_lit  *help = arg_lit0("hH", "help", "print this help and exit");
	struct arg_lit  *version = arg_lit0("v", "version,ver", "print version information and exit");
	struct arg_int *FileLogLevel = arg_int0("f", "FileLogLevel", "<int>", "Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
	struct arg_file *LogFile = arg_file0("F", "log", "<file>", "diagnostics log file");
	struct arg_int *pmode = arg_int0("m", "mode", "<int>", "Marker processing mode: 0 - default");
	struct arg_int *maxservinsts = arg_int0("n", "instances", "<int>", "max number of service instances supported (0 - defaults to number of CPU cores, max of 511)");
	struct arg_int *maxservmemgb = arg_int0("m", "mem", "<int>", "max allocatable physical memory (GB) available for all service instances (0 - defaults to 75% host memory, max of 1000)");
	struct arg_int *maxconnwait = arg_int0("w", "wait", "<int>", "wait for at most this many minutes for service requester connection (0 - defaults to 15min, max of 240)");

	struct arg_str  *host = arg_str0("u", "rmihost", "<string>", "Connect to this service requester host name or IPv4/IPv5 address (default 127.0.0.1)");
	struct arg_str  *service = arg_str0("U", "rmiservice", "<string>", "Connect to service requester on this service name or port (default 7869)");

	struct arg_end *end = arg_end(100);
	void *argtable[] = { help,version,FileLogLevel,LogFile,
		pmode,maxconnwait,maxservinsts,maxservmemgb,host,service,
		end };
	char **pAllArgs;
	int argerrors;
	argerrors = CUtility::arg_parsefromfile(argc, (char **)argv, &pAllArgs);
	if (argerrors >= 0)
		argerrors = arg_parse(argerrors, pAllArgs, argtable);

	/* special case: '--help' takes precedence over error reporting */
	if (help->count > 0)
	{
		printf("\n%s BKSTestProvider, Version %s\nOptions ---\n", gszProcName, cpszProgVer);
		arg_print_syntax(stdout, argtable, "\n");
		arg_print_glossary(stdout, argtable, "  %-25s %s\n");
		printf("\nNote: Parameters can be entered into a parameter file, one parameter per line.");
		printf("\n      To invoke this parameter file then precede it's name with '@'");
		printf("\n      e.g. %s @myparams.txt\n", gszProcName);
		printf("\nPlease report any issues regarding usage of %s to stuart.stephen@csiro.au\n\n", gszProcName);
		exit(1);
	}

	/* special case: '--version' takes precedence error reporting */
	if (version->count > 0)
	{
		printf("\n%s Version %s\n", gszProcName, cpszProgVer);
		exit(1);
	}

	if (!argerrors)
	{
		if (FileLogLevel->count && !LogFile->count)
		{
			printf("\nError: FileLogLevel '-f%d' specified but no logfile '-F<logfile>\n'", FileLogLevel->ival[0]);
			exit(1);
		}

		iScreenLogLevel = iFileLogLevel = FileLogLevel->count ? FileLogLevel->ival[0] : eDLInfo;
		if (iFileLogLevel < eDLNone || iFileLogLevel > eDLDebug)
		{
			printf("\nError: FileLogLevel '-l%d' specified outside of range %d..%d\n", iFileLogLevel, eDLNone, eDLDebug);
			exit(1);
		}

		if (LogFile->count)
		{
			strncpy(szLogFile, LogFile->filename[0], _MAX_PATH);
			szLogFile[_MAX_PATH - 1] = '\0';
		}
		else
		{
			iFileLogLevel = eDLNone;
			szLogFile[0] = '\0';
		}

		// now that log parameters have been parsed then initialise diagnostics log system
		if (!gDiagnostics.Open(szLogFile, (etDiagLevel)iScreenLogLevel, (etDiagLevel)iFileLogLevel, true))
		{
			printf("\nError: Unable to start diagnostics subsystem\n");
			if (szLogFile[0] != '\0')
				printf(" Most likely cause is that logfile '%s' can't be opened/created\n", szLogFile);
			exit(1);
		}

		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Version: %s", cpszProgVer);

		PMode = (int)(pmode->count ? pmode->ival[0] : 0);
		if (PMode < 0 || PMode > 0)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Processing mode '-m%d' specified outside of range %d..%d\n", PMode, 0, 0);
			exit(1);
			}

		MaxConnWait = (int)(maxconnwait->count ? maxconnwait->ival[0] : 15);
		if (MaxConnWait < 0 || MaxConnWait > 240)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Maximum connection wait '-w%d' specified outside of range %d..%d\n", MaxConnWait, 1,240);
			exit(1);
			}
		if(MaxConnWait == 0)
			MaxConnWait = 15;
	
		MaxServInsts = (int)(maxservinsts->count ? maxservinsts->ival[0] : 0);
		if (MaxServInsts < 0 || MaxServInsts > cMaxServiceInsts)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Maximum service instances '-n%d' specified outside of range %d..%d\n", MaxServInsts, 0,cMaxServiceInsts);
			exit(1);
			}

		MaxServMemGB = (int)(maxservmemgb->count ? maxservmemgb->ival[0] : 0);
		if (MaxServMemGB < 0 || MaxServMemGB > 1000)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Maximum allocatable physical memory '-m%d' specified outside of range %d..%d\n", MaxServMemGB, 0, 1000);
			exit(1);
			}

// show user current resource limits
#ifndef _WIN32
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "Resources: %s",CUtility::ReportResourceLimits());
#endif

		if(MaxServInsts == 0 || MaxServMemGB == 0)
			{
#ifdef _WIN32
			if(MaxServInsts == 0)
				{
				SYSTEM_INFO SystemInfo;
				GetSystemInfo(&SystemInfo);
				MaxServInsts = SystemInfo.dwNumberOfProcessors;
				}
			if(MaxServMemGB == 0)
				{
				MEMORYSTATUSEX MemStatus;
				MemStatus.dwLength = sizeof (MemStatus);
				GlobalMemoryStatusEx(&MemStatus);
				MaxServMemGB = (int)(((MemStatus.ullAvailPhys + 0x01fffffff)*(UINT64)75) / ((UINT64)100 * 0x03fffffff)); 
				}
#else
			if(MaxServInsts == 0)
				MaxServInsts = sysconf(_SC_NPROCESSORS_CONF);
			if(MaxServMemGB == 0)
				MaxServMemGB = (int)(((((UINT64)sysconf(_SC_PHYS_PAGES) * (UINT64)sysconf(_SC_PAGESIZE)) + 0x01fffffff)*(UINT64)75) / ((UINT64)100 * 0x03fffffff));
#endif
			if(MaxServMemGB < 1)
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Insufficent physical memory on this host");
				exit(1);
				}
			}



		szHostName[0] = '\0';
		if (host->count)
			{
			strncpy(szHostName, host->sval[0], sizeof(szHostName) - 1);
			szHostName[sizeof(szHostName) - 1] = '\0';
			CUtility::TrimQuotedWhitespcExtd(szHostName);
			}
		if (szHostName[0] == '\0')
			strcpy(szHostName, "127.0.0.1");

		szServiceName[0] = '\0';
		if (service->count)
			{
			strncpy(szServiceName, service->sval[0], sizeof(szServiceName) - 1);
			szServiceName[sizeof(szServiceName) - 1] = '\0';
			CUtility::TrimQuotedWhitespcExtd(szServiceName);
			}
		if (szServiceName[0] == '\0')
			strcpy(szServiceName, "7869");


		char *pszMode;
		switch(PMode) {
			case 0:		
				pszMode = (char *)"Default";
				break;
			default:
				pszMode = (char *)"Unsupported";
			}

		gDiagnostics.DiagOut(eDLInfo, gszProcName,"Processing mode: '%s'",pszMode);
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Connect to Host name or IP: '%s'", szHostName);
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Connect to service name or Port: '%s'", szServiceName);
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Wait maximum of this long for service requester connection: %d minutes", MaxConnWait);
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Maximum number of service instances: %d", MaxServInsts);
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Maximum allocatable memory for service instances: %dGB", MaxServMemGB);

		gStopWatch.Start();
		CBKSProvider *pProvider;
		int Rslt = 0;

		do {
			pProvider = new CBKSProvider;

			Rslt = pProvider->Process(MaxConnWait,MaxServInsts, MaxServInsts, MaxServMemGB,szHostName, szServiceName); 

			delete pProvider;
			}
		while(Rslt == -1 || Rslt == 0);	// iterate whilst not timed out on connection or a fatal non-recoverable error; act as a psudeo service

		gStopWatch.Stop();
		
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Exit code: %d Total processing time: %s", Rslt, gStopWatch.Read());
		exit(Rslt);
	}
	return(1);
}


CPBSWService::CPBSWService()
{
}


CPBSWService::~CPBSWService()
{
}
