// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once
#ifdef WIN32
#define _CRT_SECURE_NO_DEPRECATE 1
#define _CRT_NONSTDC_NO_DEPRECATE 1
#include <iostream>
#include <tchar.h>

// TODO: reference additional headers your program requires here

//&&Str Class Support
//Do not remove any lines from this block; refer to documentation
//  about conditional defines and their function
#define STR_WIN32
#define STR_THREAD_SAFE
#define STR_DEFINES_DONE
#define STR_USE_REGEX
#include "../libbiokanga/str.h"
//&&End of Str Class Support
#endif
