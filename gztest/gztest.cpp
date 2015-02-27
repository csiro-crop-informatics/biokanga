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
#else
#include <sys/mman.h>
#include <pthread.h>
#include "../libbiokanga/commhdrs.h"
#endif



const size_t cAllocInBuffer =   16000000;

int main(int argc, char* argv[])
{
char szgzFile[256];

INT64 TotLoaded;
size_t InBuffLen;
size_t AllocInBufferSize;
char *pszInBuffer;

printf("Starting..\n");
AllocInBufferSize = cAllocInBuffer;
pszInBuffer = (char *)malloc(AllocInBufferSize);

strcpy(szgzFile,"../../Experiments/SafflowerRNA_SN700819R/Trinity/SN700819R.fastq.gz");

gzFile RdsFile;
z_off_t FileOfs;

if((RdsFile = gzopen(szgzFile,"r"))==NULL)
	{
	printf("\nError: gzopen(%s) failed",szgzFile);
	free(pszInBuffer);
	exit(1);
	}
if(gzbuffer(RdsFile,cgzAllocInBuffer)!=0)
	{
	printf("\n gzbuffer(%d) failed",(int)cgzAllocInBuffer);
	free(pszInBuffer);
	exit(1);
	}

if(gzdirect(RdsFile))
	printf("\nLoading direct uncompressed");
else
	printf("\nLoading compressed");
TotLoaded = 0;
int Itrs = 0;
while((InBuffLen = (INT64)gzread(RdsFile,pszInBuffer,(int)AllocInBufferSize)) > 0)
	{
	TotLoaded += (INT64)InBuffLen;
	if(!(Itrs++ % 20))
		printf("\nTotLoaded: %lld",TotLoaded);

	FileOfs = gzoffset(RdsFile);
	printf("\nLoading compressed FileOfs: %ld TotLoaded: %lld ratio: %1.3f",FileOfs,TotLoaded,(double)FileOfs/TotLoaded);
	}

gzclose(RdsFile);
free(pszInBuffer);

printf("\nLoaded %lld chars",TotLoaded);
return 0;
}

