#include "stdafx.h"
#ifdef _WIN32
#include "./commhdrs.h"
#else
#include "./commhdrs.h"
#endif

// Credits: the functions in this file were simply copied with only a couple
// of compiler specific warning fixes.
// The original code was sourced from files in randomc.zip - Agner Fog
//
// See http://www.agner.org/random/ for full details

/**************************   mersenne.cpp   **********************************
* Author:        Agner Fog
* Date created:  2001
* Last modified: 2008-11-16
* Project:       randomc.h
* Platform:      Any C++
* Description:
* Random Number generator of type 'Mersenne Twister'
*
* This random number generator is described in the article by
* M. Matsumoto & T. Nishimura, in:
* ACM Transactions on Modeling and Computer Simulation,
* vol. 8, no. 1, 1998, pp. 3-30.
* Details on the initialization scheme can be found at
* http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
*
* Further documentation:
* The file ran-instructions.pdf contains further documentation and 
* instructions.
*
* Copyright 2001-2008 by Agner Fog. 
* GNU General Public License http://www.gnu.org/licenses/gpl.html
*******************************************************************************/

typedef enum  TArch {LITTLE_ENDIAN1, BIG_ENDIAN1, NONIEEE, EXTENDEDPRECISIONLITTLEENDIAN, ENDIANUNKNOWN} eTarch;
static eTarch GetEndian(void)
{
eTarch Architecture;                  // conversion to float depends on computer architecture

  // detect computer architecture
union {double f; UINT32 i[3];} convert;

convert.i[2] = 0;
convert.f = 1.0;
if (convert.i[2] == 0x3FFF)				// must be 80bit long double 
	Architecture = EXTENDEDPRECISIONLITTLEENDIAN;
else
	if (convert.i[1] == 0x3FF00000)			// must be a 64bit double using little endian
		Architecture = LITTLE_ENDIAN1;
	else 
		if (convert.i[0] == 0x3FF00000)		// must be a 64bit double using big endian		
			Architecture = BIG_ENDIAN1;
		else 
			Architecture = NONIEEE;		// haven't a clue :-)
return(Architecture);
}

static double ToDouble(UINT32 x)
{
   // The fastest way to convert random bits to floating point is as follows:
  // Set the binary exponent of a floating point number to 1+bias and set
  // the mantissa to random bits. This will give a random number in the 
  // interval [1,2). Then subtract 1.0 to get a random number in the interval
  // [0,1). This procedure requires that we know how floating point numbers
  // are stored. The storing method is tested in function RandomInit and saved 
  // in the variable Architecture. The following switch statement can be
  // omitted if the architecture is known. (generally an Intel  PC running Windows or Linux uses
  // LITTLE_ENDIAN1 architecture):
union {double f; UINT32 i[3];} convert;
  convert.f = 1.0;
  switch (GetEndian()) {
	  case EXTENDEDPRECISIONLITTLEENDIAN:		// 80bit
		convert.i[1] = 0x80000000 | (x >> 1);
		convert.i[0] =  x << 31;
		return convert.f - 1.0;

	  case LITTLE_ENDIAN1:						// 64bit little endian
		convert.i[0] =  x << 20;
		convert.i[1] = (x >> 12) | convert.i[1];
		return convert.f - 1.0;

	  case BIG_ENDIAN1:
		convert.i[0] = (x >> 12) | convert.i[1];	// 64bit big endian
		convert.i[1] =  x << 20;
		return convert.f - 1.0;

		case NONIEEE: default:
			break;
	} 
  // This somewhat slower method works for all architectures, including 
  // non-IEEE floating point representation:
  return (double)x * (1.0/((double)(UINT32)(-1L)+1.0));
}

void CRandomMersenne::Init0(int seed) {
   // Seed generator
   const uint32_t factor = 1812433253UL;
   mt[0]= seed;
   for (mti=1; mti < MERS_N; mti++) {
      mt[mti] = (factor * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
   }
}

void CRandomMersenne::RandomInit(int seed) {
   // Initialize and seed
   Init0(seed);

   // Randomize some more
   for (int i = 0; i < 37; i++) BRandom();
}


void CRandomMersenne::RandomInitByArray(int const seeds[], int NumSeeds) {
   // Seed by more than 32 bits
   int i, j, k;

   // Initialize
   Init0(19650218);

   if (NumSeeds <= 0) return;

   // Randomize mt[] using whole seeds[] array
   i = 1;  j = 0;
   k = (MERS_N > NumSeeds ? MERS_N : NumSeeds);
   for (; k; k--) {
      mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL)) + (uint32_t)seeds[j] + j;
      i++; j++;
      if (i >= MERS_N) {mt[0] = mt[MERS_N-1]; i=1;}
      if (j >= NumSeeds) j=0;}
   for (k = MERS_N-1; k; k--) {
      mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL)) - i;
      if (++i >= MERS_N) {mt[0] = mt[MERS_N-1]; i=1;}}
   mt[0] = 0x80000000UL;  // MSB is 1; assuring non-zero initial array

   // Randomize some more
   mti = 0;
   for (int i = 0; i <= MERS_N; i++) BRandom();
}


uint32_t CRandomMersenne::BRandom() {
   // Generate 32 random bits
   uint32_t y;

   if (mti >= MERS_N) {
      // Generate MERS_N words at one time
      const uint32_t LOWER_MASK = (1LU << MERS_R) - 1;       // Lower MERS_R bits
      const uint32_t UPPER_MASK = 0xFFFFFFFF << MERS_R;      // Upper (32 - MERS_R) bits
      static const uint32_t mag01[2] = {0, MERS_A};

      int kk;
      for (kk=0; kk < MERS_N-MERS_M; kk++) {    
         y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
         mt[kk] = mt[kk+MERS_M] ^ (y >> 1) ^ mag01[y & 1];}

      for (; kk < MERS_N-1; kk++) {    
         y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
         mt[kk] = mt[kk+(MERS_M-MERS_N)] ^ (y >> 1) ^ mag01[y & 1];}      

      y = (mt[MERS_N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
      mt[MERS_N-1] = mt[MERS_M-1] ^ (y >> 1) ^ mag01[y & 1];
      mti = 0;
   }
   y = mt[mti++];

   // Tempering (May be omitted):
   y ^=  y >> MERS_U;
   y ^= (y << MERS_S) & MERS_B;
   y ^= (y << MERS_T) & MERS_C;
   y ^=  y >> MERS_L;

   return y;
}


double CRandomMersenne::Random() {
   // Output random float number in the interval 0 <= x < 1
   // Multiply by 2^(-32)
   return (double)BRandom() * (1./(65536.*65536.));
}


int CRandomMersenne::IRandom(int min, int max) {
   // Output random integer in the interval min <= x <= max
   // Relative error on frequencies < 2^-32
   if (max <= min) {
      if (max == min) return min; else return 0x80000000;
   }
   // Multiply interval with random and truncate
   int r = int((double)(uint32_t)(max - min + 1) * Random() + min); 
   if (r > max) r = max;
   return r;
}


int CRandomMersenne::IRandomX(int min, int max) {
   // Output random integer in the interval min <= x <= max
   // Each output value has exactly the same probability.
   // This is obtained by rejecting certain bit values so that the number
   // of possible bit values is divisible by the interval length
   if (max <= min) {
      if (max == min) return min; else return 0x80000000;
   }
#ifdef  INT64_SUPPORTED
   // 64 bit integers available. Use multiply and shift method
   uint32_t interval;                    // Length of interval
   uint64_t longran;                     // Random bits * interval
   uint32_t iran;                        // Longran / 2^32
   uint32_t remainder;                   // Longran % 2^32

   interval = uint32_t(max - min + 1);
   if (interval != LastInterval) {
      // Interval length has changed. Must calculate rejection limit
      // Reject when remainder >= 2^32 / interval * interval
      // RLimit will be 0 if interval is a power of 2. No rejection then
      RLimit = uint32_t(((uint64_t)1 << 32) / interval) * interval - 1;
      LastInterval = interval;
   }
   do { // Rejection loop
      longran  = (uint64_t)BRandom() * interval;
      iran = (uint32_t)(longran >> 32);
      remainder = (uint32_t)longran;
   } while (remainder > RLimit);
   // Convert back to signed and return result
   return (int32_t)iran + min;

#else
   // 64 bit integers not available. Use modulo method
   uint32_t interval;                    // Length of interval
   uint32_t bran;                        // Random bits
   uint32_t iran;                        // bran / interval
   uint32_t remainder;                   // bran % interval

   interval = uint32_t(max - min + 1);
   if (interval != LastInterval) {
      // Interval length has changed. Must calculate rejection limit
      // Reject when iran = 2^32 / interval
      // We can't make 2^32 so we use 2^32-1 and correct afterwards
      RLimit = (uint32_t)0xFFFFFFFF / interval;
      if ((uint32_t)0xFFFFFFFF % interval == interval - 1) RLimit++;
   }
   do { // Rejection loop
      bran = BRandom();
      iran = bran / interval;
      remainder = bran % interval;
   } while (iran >= RLimit);
   // Convert back to signed and return result
   return (int32_t)remainder + min;

#endif
}


/**************************   mother.cpp   ************************************
* Author:        Agner Fog
* Date created:  1999
* Last modified: 2008-11-16
* Project:       randomc.h
* Platform:      This implementation uses 64-bit integers for intermediate calculations.
*                Works only on compilers that support 64-bit integers.
* Description:
* Random Number generator of type 'Mother-Of-All generator'.
*
* This is a multiply-with-carry type of random number generator
* invented by George Marsaglia.  The algorithm is:             
* S = 2111111111*X[n-4] + 1492*X[n-3] + 1776*X[n-2] + 5115*X[n-1] + C
* X[n] = S modulo 2^32
* C = floor(S / 2^32)
*
* Further documentation:
* The file ran-instructions.pdf contains further documentation and 
* instructions.
*
* Copyright 1999-2008 by Agner Fog. 
* GNU General Public License http://www.gnu.org/licenses/gpl.html
******************************************************************************/

// Output random bits
uint32_t CRandomMother::BRandom() {
  uint64_t sum;
  sum = (uint64_t)2111111111UL * (uint64_t)x[3] +
     (uint64_t)1492 * (uint64_t)(x[2]) +
     (uint64_t)1776 * (uint64_t)(x[1]) +
     (uint64_t)5115 * (uint64_t)(x[0]) +
     (uint64_t)x[4];
  x[3] = x[2];  x[2] = x[1];  x[1] = x[0];
  x[4] = (uint32_t)(sum >> 32);                  // Carry
  x[0] = (uint32_t)sum;                          // Low 32 bits of sum
  return x[0];
} 


// returns a random number between 0 and 1:
double CRandomMother::Random() {
   return (double)BRandom() * (1./(65536.*65536.));
}


// returns integer random number in desired interval:
int CRandomMother::IRandom(int min, int max) {
   // Output random integer in the interval min <= x <= max
   // Relative error on frequencies < 2^-32
   if (max <= min) {
      if (max == min) return min; else return 0x80000000;
   }
   // Assume 64 bit integers supported. Use multiply and shift method
   uint32_t interval;                  // Length of interval
   uint64_t longran;                   // Random bits * interval
   uint32_t iran;                      // Longran / 2^32

   interval = (uint32_t)(max - min + 1);
   longran  = (uint64_t)BRandom() * interval;
   iran = (uint32_t)(longran >> 32);
   // Convert back to signed and return result
   return (int32_t)iran + min;
}


// this function initializes the random number generator:
void CRandomMother::RandomInit (int seed) {
  int i;
  uint32_t s = seed;
  // make random numbers and put them into the buffer
  for (i = 0; i < 5; i++) {
    s = s * 29943829 - 1;
    x[i] = s;
  }
  // randomize some more
  for (i=0; i<19; i++) BRandom();
}

