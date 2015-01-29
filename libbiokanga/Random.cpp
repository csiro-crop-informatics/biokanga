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

/************************** MERSENNE.CPP ******************** AgF 2001-10-18 *
*  Random Number generator 'Mersenne Twister'                                *
*                                                                            *
*  This random number generator is described in the article by               *
*  M. Matsumoto & T. Nishimura, in:                                          *
*  ACM Transactions on Modeling and Computer Simulation,                     *
*  vol. 8, no. 1, 1998, pp. 3-30.                                            *
*  Details on the initialization scheme can be found at                      *
*  http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html                  *
*                                                                            *
*  Experts consider this an excellent random number generator.               *
*                                                                            *
*  © 2001 - 2004 A. Fog.                                                     *
*  GNU General Public License www.gnu.org/copyleft/gpl.html                  *
*****************************************************************************/

void TRandomMersenne::RandomInit(uint32 seed) {
  // re-seed generator
  mt[0]= seed;
  for (mti=1; mti < MERS_N; mti++) {
    mt[mti] = (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);}

  // detect computer architecture
  union {double f; uint32 i[2];} convert;
  convert.f = 1.0;
  // Note: Old versions of the Gnu g++ compiler may make an error here,
  // compile with the option  -fenum-int-equiv  to fix the problem
  if (convert.i[1] == 0x3FF00000) Architecture = LITTLE_ENDIAN1;
  else if (convert.i[0] == 0x3FF00000) Architecture = BIG_ENDIAN1;
  else Architecture = NONIEEE;}

  
void TRandomMersenne::RandomInitByArray(uint32 seeds[], int length) {
  // seed by more than 32 bits
  int i, j, k;
  RandomInit(19650218UL);
  if (length <= 0) return;
  i = 1;  j = 0;
  k = (MERS_N > length ? MERS_N : length);
  for (; k; k--) {
    mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL)) + seeds[j] + j;
    i++; j++;
    if (i >= MERS_N) {mt[0] = mt[MERS_N-1]; i=1;}
    if (j >= length) j=0;}
  for (k = MERS_N-1; k; k--) {
    mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL)) - i;
    if (++i >= MERS_N) {mt[0] = mt[MERS_N-1]; i=1;}}
  mt[0] = 0x80000000UL;} // MSB is 1; assuring non-zero initial array

  
uint32 TRandomMersenne::BRandom() {
  // generate 32 random bits
  uint32 y;

  if (mti >= MERS_N) {
    // generate MERS_N words at one time
    const uint32 LOWER_MASK = (1LU << MERS_R) - 1; // lower MERS_R bits
    const uint32 UPPER_MASK = -1L  << MERS_R;      // upper (32 - MERS_R) bits
    static const uint32 mag01[2] = {0, MERS_A};
    
    int kk;
    for (kk=0; kk < MERS_N-MERS_M; kk++) {    
      y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
      mt[kk] = mt[kk+MERS_M] ^ (y >> 1) ^ mag01[y & 1];}

    for (; kk < MERS_N-1; kk++) {    
      y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
      mt[kk] = mt[kk+(MERS_M-MERS_N)] ^ (y >> 1) ^ mag01[y & 1];}      

    y = (mt[MERS_N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
    mt[MERS_N-1] = mt[MERS_M-1] ^ (y >> 1) ^ mag01[y & 1];
    mti = 0;}

  y = mt[mti++];

  // Tempering (May be omitted):
  y ^=  y >> MERS_U;
  y ^= (y << MERS_S) & MERS_B;
  y ^= (y << MERS_T) & MERS_C;
  y ^=  y >> MERS_L;
  return y;}

  
double TRandomMersenne::Random() {
  // output random float number in the interval 0 <= x < 1
  union {double f; uint32 i[2];} convert;
  uint32 r = BRandom(); // get 32 random bits
  // The fastest way to convert random bits to floating point is as follows:
  // Set the binary exponent of a floating point number to 1+bias and set
  // the mantissa to random bits. This will give a random number in the 
  // interval [1,2). Then subtract 1.0 to get a random number in the interval
  // [0,1). This procedure requires that we know how floating point numbers
  // are stored. The storing method is tested in function RandomInit and saved 
  // in the variable Architecture. The following switch statement can be
  // omitted if the architecture is known. (A PC running Windows or Linux uses
  // LITTLE_ENDIAN1 architecture):
  switch (Architecture) {
  case LITTLE_ENDIAN1:
    convert.i[0] =  r << 20;
    convert.i[1] = (r >> 12) | 0x3FF00000;
    return convert.f - 1.0;
  case BIG_ENDIAN1:
    convert.i[1] =  r << 20;
    convert.i[0] = (r >> 12) | 0x3FF00000;
    return convert.f - 1.0;
  case NONIEEE: default:
  ;} 
  // This somewhat slower method works for all architectures, including 
  // non-IEEE floating point representation:
  return (double)r * (1./((double)(uint32)(-1L)+1.));}

  
int TRandomMersenne::IRandom(int min, int max) {
  // output random integer in the interval min <= x <= max
  int r; 
  r = int((max - min + 1) * Random()) + min; // multiply interval with random and truncate
  if (r > max) r = max;
  if (max < min) return 0x80000000;
  return r;}

/************************** MOTHER.CPP ****************** AgF 1999-03-03 *
*  'Mother-of-All' random number generator                               *
*                                                                        *
*  This is a multiply-with-carry type of random number generator         *
*  invented by George Marsaglia.  The algorithm is:                      *
*  S = 2111111111*X[n-4] + 1492*X[n-3] + 1776*X[n-2] + 5115*X[n-1] + C   *
*  X[n] = S modulo 2^32                                                  *
*  C = floor(S / 2^32)                                                   *
*                                                                        *
*  IMPORTANT:
*  This implementation uses a long double for C. Note that not all       *
*  computers and compilers support the long double (80-bit) floating     *
*  point format. It is recommended to use a Borland or Gnu compiler on   *
*  a PC. The Microsoft compiler doesn't support the long double format.  *
*  You will get an error message if your system doesn't support this.    *
*                                                                        *
* © 2002 A. Fog. GNU General Public License www.gnu.org/copyleft/gpl.html*
*************************************************************************/

// constructor:
TRandomMotherOfAll::TRandomMotherOfAll(uint32 seed) {
  // Check that compiler supports 80-bit long double:
  assert(sizeof(long double)>9);
  // initialize
  RandomInit(seed);}


// returns a random number between 0 and 1:
double TRandomMotherOfAll::Random() {
  long double c;
  c = (long double)2111111111.0 * x[3] +
      1492.0 * (x[3] = x[2]) +
      1776.0 * (x[2] = x[1]) +
      5115.0 * (x[1] = x[0]) +
      x[4];

// sjs: on IBM p690 seems to support long doubles but __NO_LONG_DOUBLE_MATH has been defined????
// if long doubles truely not supported then there would have been an earlier assert when sizeof(long double) > 9 checked
//
#ifdef __NO_LONG_DOUBLE_MATH
  x[4] = floor(c);
#else
  x[4] = floorl(c);
#endif

  x[0] = c - x[4];
  x[4] = x[4] * (1./(65536.*65536.));
  return x[0];}


// returns integer random number in desired interval:
int TRandomMotherOfAll::IRandom(int min, int max) {
  int iinterval = max - min + 1;
  if (iinterval <= 0) return 0x80000000; // error
  int i = int(iinterval * Random());     // truncate
  if (i >= iinterval) i = iinterval-1;
  return min + i;}


// this function initializes the random number generator:
void TRandomMotherOfAll::RandomInit (uint32 seed) {
  int i;
  uint32 s = seed;
  // make random numbers and put them into the buffer
  for (i=0; i<5; i++) {
    s = s * 29943829 - 1;
    x[i] = s * (1./(65536.*65536.));}
  // randomize some more
  for (i=0; i<19; i++) Random();}


/************************* RANROTB.CPP ****************** AgF 1999-03-03 *
*  Random Number generator 'RANROT' type B                               *
*                                                                        *
*  This is a lagged-Fibonacci type of random number generator with       *
*  rotation of bits.  The algorithm is:                                  *
*  X[n] = ((X[n-j] rotl r1) + (X[n-k] rotl r2)) modulo 2^b               *
*                                                                        *
*  The last k values of X are stored in a circular buffer named          *
*  randbuffer.                                                           *
*  The code includes a self-test facility which will detect any          *
*  repetition of previous states.                                        *
*                                                                        *
*  The theory of the RANROT type of generators and the reason for the    *
*  self-test are described at www.agner.org/random                       *
*                                                                        *
* © 2002 A. Fog. GNU General Public License www.gnu.org/copyleft/gpl.html*
*************************************************************************/



// If your system doesn't have a rotate function for 32 bits integers,
// then use the definition below. If your system has the _lrotl function 
// then remove this.
#ifndef _WIN32
uint32 _lrotl (uint32 x, int r) {
	return (x << r) | (x >> (sizeof(x)*8-r));}
#endif

// constructor:
TRanrotBGenerator::TRanrotBGenerator(uint32 seed) {
  RandomInit(seed);
  // detect computer architecture
  union {double f; uint32 i[2];} convert;
  convert.f = 1.0;
  if (convert.i[1] == 0x3FF00000) Architecture = LITTLE_ENDIAN1;
  else if (convert.i[0] == 0x3FF00000) Architecture = BIG_ENDIAN1;
  else Architecture = NONIEEE;}


// returns a random number between 0 and 1:
double TRanrotBGenerator::Random() {
  uint32 x;
  // generate next random number
  x = randbuffer[p1] = _lrotl(randbuffer[p2], R1) + _lrotl(randbuffer[p1], R2);
  // rotate list pointers
  if (--p1 < 0) p1 = KK - 1;
  if (--p2 < 0) p2 = KK - 1;
  // perform self-test
  if (randbuffer[p1] == randbufcopy[0] &&
    memcmp(randbuffer, randbufcopy+KK-p1, KK*sizeof(uint32)) == 0) {
      // self-test failed
      if ((p2 + KK - p1) % KK != JJ) {
        // note: the way of printing error messages depends on system
        // In Windows you may use FatalAppExit
        printf("Random number generator not initialized");}
      else {
        printf("Random number generator returned to initial state");}
      exit(1);}
  // conversion to float:
  union {double f; uint32 i[2];} convert;
  switch (Architecture) {
  case LITTLE_ENDIAN1:
    convert.i[0] =  x << 20;
    convert.i[1] = (x >> 12) | 0x3FF00000;
    return convert.f - 1.0;
  case BIG_ENDIAN1:
    convert.i[1] =  x << 20;
    convert.i[0] = (x >> 12) | 0x3FF00000;
    return convert.f - 1.0;
  case NONIEEE: default:
  ;} 
  // This somewhat slower method works for all architectures, including 
  // non-IEEE floating point representation:
  return (double)x * (1./((double)(uint32)(-1L)+1.));}


// returns integer random number in desired interval:
int TRanrotBGenerator::IRandom(int min, int max) {
  int iinterval = max - min + 1;
  if (iinterval <= 0) return -1; // error
  int i = (int)(iinterval * Random()); // truncate
  if (i >= iinterval) i = iinterval-1;
  return min + i;}
  

void TRanrotBGenerator::RandomInit (uint32 seed) {
  // this function initializes the random number generator.
  int i;
  uint32 s = seed;

  // make random numbers and put them into the buffer
  for (i=0; i<KK; i++) {
    s = s * 2891336453 + 1;
    randbuffer[i] = s;}

  // check that the right data formats are used by compiler:
  union {
    double randp1;
    uint32 randbits[2];};
  randp1 = 1.5;
  assert(randbits[1]==0x3FF80000); // check that IEEE double precision float format used

  // initialize pointers to circular buffer
  p1 = 0;  p2 = JJ;
  // store state for self-test
  memmove (randbufcopy, randbuffer, KK*sizeof(uint32));
  memmove (randbufcopy+KK, randbuffer, KK*sizeof(uint32));
  // randomize some more
  for (i=0; i<9; i++) Random();
}

/************************* RANROTW.CPP ****************** AgF 1999-03-03 *
*  Random Number generator 'RANROT' type W                               *
*  This version is used when a resolution higher that 32 bits is desired.*
*                                                                        *
*  This is a lagged-Fibonacci type of random number generator with       *
*  rotation of bits.  The algorithm is:                                  *
*  Z[n] = (Y[n-j] + (Y[n-k] rotl r1)) modulo 2^(b/2)                     *
*  Y[n] = (Z[n-j] + (Z[n-k] rotl r2)) modulo 2^(b/2)                     *
*  X[n] = Y[n] + Z[n]*2^(b/2)                                            *
*                                                                        *
*  The last k values of Y and Z are stored in a circular buffer named    *
*  randbuffer.                                                           *
*  The code includes a self-test facility which will detect any          *
*  repetition of previous states.                                        *
*  The function uses a fast method for conversion to floating point.     *
*  This method relies on floating point numbers being stored in the      *
*  standard 64-bit IEEE format or the 80-bit long double format.         *
*                                                                        *
*  The theory of the RANROT type of generators and the reason for the    *
*  self-test are described at www.agner.org/random/theory                *
*                                                                        *
* ©2002 A. Fog. GNU General Public License www.gnu.org/copyleft/gpl.html *
*************************************************************************/

// If your system doesn't have a rotate function for 32 bits integers,
// then use the definition below. If your system has the _lrotl function 
// then remove this.
// uint32 _lrotl (uint32 x, int r) {
//   return (x << r) | (x >> (sizeof(x)*8-r));}


// constructor:
TRanrotWGenerator::TRanrotWGenerator(uint32 seed) {
  RandomInit(seed);
  // detect computer architecture
  randbits[2] = 0; randp1 = 1.0;
  if (randbits[2] == 0x3FFF) Architecture = EXTENDEDPRECISIONLITTLEENDIAN;
  else if (randbits[1] == 0x3FF00000) Architecture = LITTLE_ENDIAN1;
  else if (randbits[0] == 0x3FF00000) Architecture = BIG_ENDIAN1;
  else Architecture = NONIEEE;}


uint32 TRanrotWGenerator::BRandom() {
  // generate next random number
  uint32 y, z;
  // generate next number
  z = _lrotl(randbuffer[p1][0], R1) + randbuffer[p2][0];
  y = _lrotl(randbuffer[p1][1], R2) + randbuffer[p2][1];
  randbuffer[p1][0] = y; randbuffer[p1][1] = z;
  // rotate list pointers
  if (--p1 < 0) p1 = KK - 1;
  if (--p2 < 0) p2 = KK - 1;
  // perform self-test
  if (randbuffer[p1][0] == randbufcopy[0][0] &&
    memcmp(randbuffer, randbufcopy[KK-p1], 2*KK*sizeof(uint32)) == 0) {
      // self-test failed
      if ((p2 + KK - p1) % KK != JJ) {
        // note: the way of printing error messages depends on system
        // In Windows you may use FatalAppExit
        printf("Random number generator not initialized");}
      else {
        printf("Random number generator returned to initial state");}
      exit(1);}
  randbits[0] = y;
  randbits[1] = z;
  return y;}


long double TRanrotWGenerator::Random() {
  // returns a random number between 0 and 1.
  uint32 z = BRandom();  // generate 64 random bits
  switch (Architecture) {
  case EXTENDEDPRECISIONLITTLEENDIAN:
    // 80 bits floats = 63 bits resolution
    randbits[1] = z | 0x80000000;  break;
  case LITTLE_ENDIAN1:
    // 64 bits floats = 52 bits resolution
    randbits[1] = (z & 0x000FFFFF) | 0x3FF00000;  break;
  case BIG_ENDIAN1:
    // 64 bits floats = 52 bits resolution
    randbits[0] = (randbits[0] & 0x000FFFFF) | 0x3FF00000;  break;
  case NONIEEE: default:
    // not a recognized floating point format. 32 bits resolution
    return (double)z * (1./((double)(uint32)(-1L)+1.));}
  return randp1 - 1.0;}


int TRanrotWGenerator::IRandom(int min, int max) {
  // get integer random number in desired interval
  int iinterval = max - min + 1;
  if (iinterval <= 0) return 0x80000000;  // error
  int i = int(iinterval * Random());      // truncate
  if (i >= iinterval) i = iinterval-1;
  return min + i;}


void TRanrotWGenerator::RandomInit (uint32 seed) {
  // this function initializes the random number generator.
  int i, j;

  // make random numbers and put them into the buffer
  for (i=0; i<KK; i++) {
    for (j=0; j<2; j++) {
      seed = seed * 2891336453UL + 1;
      randbuffer[i][j] = seed;}}
  // set exponent of randp1
  randbits[2] = 0; randp1 = 1.0;
  // initialize pointers to circular buffer
  p1 = 0;  p2 = JJ;
  // store state for self-test
  memmove (randbufcopy, randbuffer, 2*KK*sizeof(uint32));
  memmove (randbufcopy[KK], randbuffer, 2*KK*sizeof(uint32));
  // randomize some more
  for (i=0; i<31; i++) BRandom();
}





