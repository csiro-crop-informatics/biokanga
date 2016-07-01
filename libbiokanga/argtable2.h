/*********************************************************************
This file is part of the argtable2 library.
Copyright (C) 1998,1999,2000,2001,2003 Stewart Heitmann
sheitmann@users.sourceforge.net

The argtable2 library is free software; you can redistribute it and/or
modify it under the terms of the GNU Library General Public License as
published by the Free Software Foundation; either version 2 of the
License, or (at your option) any later version.

This software is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Library General Public License for more details.

You should have received a copy of the GNU Library General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
USA.
**********************************************************************/
#pragma once

#include <stdio.h>    /* typedef FILE */

#ifdef __cplusplus
extern "C" {
#endif


/* bit masks for arg_hdr.flag */
enum
    {
    ARG_TERMINATOR=0x1,
    ARG_HASVALUE=0x2
    };

typedef void (arg_resetfn)(void *parent);
typedef int  (arg_scanfn)(void *parent, const char *argval);
typedef int  (arg_checkfn)(void *parent);
typedef void (arg_errorfn)(void *parent, FILE *fp, int error, const char *argval, const char *progname);


/*
 * The arg_hdr struct defines properties that are common to all arg_xxx structs.
 * The argtable library requires each arg_xxx struct to have an arg_hdr
 * struct as its first data member.
 * The argtable library functions then use this data to identify the
 * properties of the command line option, such as its option tags,
 * datatype string, and glossary strings, and so on.
 * Moreover, the arg_hdr struct contains pointers to custom functions that
 * are provided by each arg_xxx struct which perform the tasks of parsing
 * that particular arg_xxx arguments, performing post-parse checks, and
 * reporting errors.
 * These functions are private to the individual arg_xxx source code
 * and are the pointer to them are initiliased by that arg_xxx struct's
 * constructor function. The user could alter them after construction
 * if desired, but the original intention is for them to be set by the
 * constructor and left unaltered.
 */
struct arg_hdr
   {
   char         flag;        /* Modifier flags: ARG_TERMINATOR, ARG_HASVALUE. */
   const char  *shortopts;   /* String defining the short options */
   const char  *longopts;    /* String defiing the long options */
   const char  *datatype;    /* Description of the argument data type */
   const char  *glossary;    /* Description of the option as shown by arg_print_glossary function */
   int          mincount;    /* Minimum number of occurences of this option accepted */
   int          maxcount;    /* Maximum number of occurences if this option accepted */
   void        *parent;      /* Pointer to parent arg_xxx struct */
   arg_resetfn *resetfn;     /* Pointer to parent arg_xxx reset function */
   arg_scanfn  *scanfn;      /* Pointer to parent arg_xxx scan function */
   arg_checkfn *checkfn;     /* Pointer to parent arg_xxx check function */
   arg_errorfn *errorfn;     /* Pointer to parent arg_xxx error function */
   };

struct arg_rem
   {
   struct arg_hdr hdr;      /* The mandatory argtable header struct */
   };

struct arg_lit
   {
   struct arg_hdr hdr;      /* The mandatory argtable header struct */
   int count;               /* Number of occurences of this parsed */
   };

struct arg_int
   {
   struct arg_hdr hdr;      /* The mandatory argtable header struct */
   int count;               /* Number of occurences of this parsed */
   int *ival;               /* Storage for integer argument values */
   };

struct arg_dbl
   {
   struct arg_hdr hdr;      /* The mandatory argtable header struct */
   int count;               /* Number of occurences of this parsed */
   double *dval;            /* Storage for double argument values */
   };

struct arg_str
   {
   struct arg_hdr hdr;      /* The mandatory argtable header struct */
   int count;               /* Number of occurences of this parsed */
   const char **sval;       /* Storage for string argument pointers */
   };

struct arg_file
   {
   struct arg_hdr hdr;      /* The mandatory argtable header struct */
   int count;               /* Number of occurences of this parsed */
   const char **filename;   /* Storage for filename string pointers */
   const char **basename;   /* Storage for basename string pointers */
   const char **extension;  /* Storage for extension string pointers */
   };


enum {ARG_ELIMIT=1, ARG_EMALLOC, ARG_ENOMATCH, ARG_ELONGOPT, ARG_EMISSARG};
struct arg_end
   {
   struct arg_hdr hdr;      /* The mandatory argtable header struct */
   int count;               /* Number of occurences of this parsed */
   int *error;              /* Storage for error codes */
   void **parent;           /* Storage for pointers to parent arg_xxx */
   const char **argval;     /* Storage for pointers to offending command line string */
   };


/**** arg_xxx constructor functions *********************************/

struct arg_rem* arg_rem(const char* datatype, const char* glossary);

struct arg_lit* arg_lit0(const char* shortopts,
                         const char* longopts,
                         const char* glossary);
struct arg_lit* arg_lit1(const char* shortopts,
                         const char* longopts,
                         const char *glossary);
struct arg_lit* arg_litn(const char* shortopts,
                         const char* longopts,
                         int mincount,
                         int maxcount,
                         const char *glossary);

struct arg_int* arg_int0(const char* shortopts,
                         const char* longopts,
                         const char* datatype,
                         const char* glossary);
struct arg_int* arg_int1(const char* shortopts,
                         const char* longopts,
                         const char* datatype,
                         const char *glossary);
struct arg_int* arg_intn(const char* shortopts,
                         const char* longopts,
                         const char *datatype,
                         int mincount,
                         int maxcount,
                         const char *glossary);

struct arg_dbl* arg_dbl0(const char* shortopts,
                         const char* longopts,
                         const char* datatype,
                         const char* glossary);
struct arg_dbl* arg_dbl1(const char* shortopts,
                         const char* longopts,
                         const char* datatype,
                         const char *glossary);
struct arg_dbl* arg_dbln(const char* shortopts,
                         const char* longopts,
                         const char *datatype,
                         int mincount,
                         int maxcount,
                         const char *glossary);

struct arg_str* arg_str0(const char* shortopts,
                         const char* longopts,
                         const char* datatype,
                         const char* glossary);
struct arg_str* arg_str1(const char* shortopts,
                         const char* longopts,
                         const char* datatype,
                         const char *glossary);
struct arg_str* arg_strn(const char* shortopts,
                         const char* longopts,
                         const char* datatype,
                         int mincount,
                         int maxcount,
                         const char *glossary);

struct arg_file* arg_file0(const char* shortopts,
                           const char* longopts,
                           const char* datatype,
                           const char* glossary);
struct arg_file* arg_file1(const char* shortopts,
                           const char* longopts,
                           const char* datatype,
                           const char *glossary);
struct arg_file* arg_filen(const char* shortopts,
                           const char* longopts,
                           const char* datatype,
                           int mincount,
                           int maxcount,
                           const char *glossary);

struct arg_end* arg_end(int maxerrors);


/**** other functions *******************************************/
int arg_nullcheck(void **argtable);
int arg_parse(int argc, char **argv, void **argtable);
void arg_print_option(FILE *fp, const char *shortopts, const char *longopts, const char *datatype, const char *suffix);
void arg_print_syntax(FILE *fp, void **argtable, const char *suffix);
void arg_print_syntaxv(FILE *fp, void **argtable, const char *suffix);
void arg_print_glossary(FILE *fp, void **argtable, const char *format);
void arg_print_errors(FILE* fp, struct arg_end* end, const char* progname);
void arg_free(void **argtable);

#ifdef __cplusplus
}
#endif

/* Declarations for getopt.
   Copyright (C) 1989,90,91,92,93,94,96,97 Free Software Foundation, Inc.

   This file is part of the GNU C Library.  Its master source is NOT part of
   the C library, however.  The master source lives in /gd/gnu/lib.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Library General Public License as
   published by the Free Software Foundation; either version 2 of the
   License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Library General Public License for more details.

   You should have received a copy of the GNU Library General Public
   License along with the GNU C Library; see the file COPYING.LIB.  If not,
   write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330,
   Boston, MA 02111-1307, USA.  */

#pragma once

#ifdef	__cplusplus
extern "C"
{
#endif

/* For communication from `getopt' to the caller.
   When `getopt' finds an option that takes an argument,
   the argument value is returned here.
   Also, when `ordering' is RETURN_IN_ORDER,
   each non-option ARGV-element is returned here.  */

	extern char *optarg;

/* Index in ARGV of the next element to be scanned.
   This is used for communication to and from the caller
   and for communication between successive calls to `getopt'.

   On entry to `getopt', zero means this is the first call; initialize.

   When `getopt' returns -1, this is the index of the first of the
   non-option elements that the caller should itself scan.

   Otherwise, `optind' communicates from one call to the next
   how much of ARGV has been scanned so far.  */

	extern int optind;

/* Callers store zero here to inhibit the error message `getopt' prints
   for unrecognized options.  */

	extern int opterr;

/* Set to an option character which was unrecognized.  */

	extern int optopt;

/* Describe the long-named options requested by the application.
   The LONG_OPTIONS argument to getopt_long or getopt_long_only is a vector
   of `struct option' terminated by an element containing a name which is
   zero.

   The field `has_arg' is:
   no_argument          (or 0) if the option does not take an argument,
   required_argument    (or 1) if the option requires an argument,
   optional_argument    (or 2) if the option takes an optional argument.

   If the field `flag' is not NULL, it points to a variable that is set
   to the value given in the field `val' when the option is found, but
   left unchanged if the option is not found.

   To have a long-named option do something other than set an `int' to
   a compiled-in constant, such as set a value from `optarg', set the
   option's `flag' field to zero and its `val' field to a nonzero
   value (the equivalent single-letter option character, if there is
   one).  For long options that have a zero `flag' field, `getopt'
   returns the contents of the `val' field.  */

	struct option
	{
	const char *name;
		/* has_arg can't be an enum because some compilers complain about
		   type mismatches in all the code that assumes it is an int.  */
	int has_arg;
	int *flag;
		int val;
	};

/* Names for the values of the `has_arg' field of `struct option'.  */

#define	no_argument		0
#define required_argument	1
#define optional_argument	2

	extern int getopt(int argc, char *const *argv, const char *shortopts);
	extern int getopt_long(int argc, char *const *argv, const char *shortopts,
			       const struct option *longopts, int *longind);
	extern int getopt_long_only(int argc, char *const *argv,
				    const char *shortopts,
			       const struct option *longopts, int *longind);

/* Internal only.  Users should not call this directly.  */
	extern int _getopt_internal(int argc, char *const *argv,
						const char *shortopts,const struct option *longopts, int *longind,int long_only);

#ifdef	__cplusplus
}
#endif


