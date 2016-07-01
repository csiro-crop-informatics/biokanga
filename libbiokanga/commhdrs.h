#pragma once
// global headers in standard places used by this library
// local headers specific to this library are in "./conservlib.h" 
#include "stdafx.h"
#ifdef _WIN32
#include <io.h>
#else
// allows loading on older releases of Ubuntu back to 10.04
// __asm__(".symver memcpy,memcpy@GLIBC_2.11.1");
#define __USE_XOPEN_EXTENDED 1
#ifndef _LARGEFILE64_SOURCE
#define _LARGEFILE64_SOURCE 1
#endif
#ifndef _FILE_OFFSET_BITS
#define _FILE_OFFSET_BITS 64
#endif
#include <sys/time.h>
#include <sys/times.h>
#include <unistd.h>
#include <errno.h>
#include <regex.h>
#include <stdarg.h>
#include <ftw.h>
#include <fnmatch.h>
#endif
#include <sys/timeb.h>
#include <math.h>
#include <fcntl.h>      
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <search.h>
#include <limits.h>
#include <float.h>
#include <stddef.h>
#include <assert.h>
#include "./conservlib.h"

