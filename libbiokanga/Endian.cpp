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

CEndian::CEndian(void)
{
INT16 EndianTest = 0x4321;
m_bIsBigEndian = (*(char *)&EndianTest == 0x21) ? false : true;
}

CEndian::~CEndian(void)
{
}

// Swaps 16bit (2byte) endian order
UINT16 
CEndian::SwapUI16Endians(UINT16 Val)
{
return (((Val&0x00FF)<< 8)+((Val&0xFF00)>>8));
}

// Swaps 32bit (4byte) endian order
UINT32 
CEndian::SwapUI32Endians(UINT32 Val)
{
return (((Val&0x000000FF)<<24)+((Val&0x0000FF00)<<8)+
   ((Val&0x00FF0000)>>8)+((Val&0xFF000000)>>24));
}

// Swaps 64bit (8byte) endian order
UINT64 
CEndian::SwapUI64Endians(UINT64 Val)
{
UINT32 LoInt = SwapUI32Endians((UINT32)(Val >> 32));
UINT32 HiInt = SwapUI32Endians((UINT32)Val);
return(((UINT64)HiInt << 32) + (UINT64)LoInt);
}

