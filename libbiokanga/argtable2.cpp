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
#include "stdafx.h"
#ifdef _WIN32
#include "./commhdrs.h"
#else
#include "./commhdrs.h"
#endif


#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

/* local error codes */
enum {EMINCOUNT=1,EMAXCOUNT,EBADINT,EBADDOUBLE};

static
void arg_register_error(struct arg_end *end, void *parent, int error, const char *argval)
    {
    /* printf("arg_register_error(%p,%p,%d,%s)\n",end,parent,error,argval); */
    if (end->count < end->hdr.maxcount)
        {
        end->error[end->count] = error;
        end->parent[end->count] = parent;
        end->argval[end->count] = argval;
        end->count++;
        }
    else
        {
        end->error[end->hdr.maxcount-1]  = ARG_ELIMIT;
        end->parent[end->hdr.maxcount-1] = end;
        end->argval[end->hdr.maxcount-1] = NULL;
        }
    }


/*
 * Return index of first table entry with a matching short option
 * or -1 if no match was found.
 */
static
int find_shortoption(struct arg_hdr **table, char shortopt)
    {
    int tabindex;
    for(tabindex=0; !(table[tabindex]->flag&ARG_TERMINATOR); tabindex++)
        {
        if (table[tabindex]->shortopts && strchr(table[tabindex]->shortopts,shortopt))
            return tabindex;
        }
    return -1;
    }


struct longoptions
    {
    int getoptval;
    int noptions;
    struct option *options;
    };

#ifndef NDEBUG
static
void dump_longoptions(struct longoptions* longoptions)
    {
    int i;
    printf("getoptval = %d\n", longoptions->getoptval);
    printf("noptions  = %d\n", longoptions->noptions);
    for (i=0; i<longoptions->noptions; i++)
        {
        printf("options[%d].name    = \"%s\"\n", i, longoptions->options[i].name);
        printf("options[%d].has_arg = %d\n", i, longoptions->options[i].has_arg);
        printf("options[%d].flag    = %p\n", i, longoptions->options[i].flag);
        printf("options[%d].val     = %d\n", i, longoptions->options[i].val);
        }
    }
#endif

static
struct longoptions* alloc_longoptions(struct arg_hdr **table)
    {
    struct longoptions *result;
    size_t nbytes;
    int noptions = 1;
    size_t longoptlen = 0;
    int tabindex;

    /*
     * Determine the total number of option structs required
     * by counting the number of comma separated long options
     * in all table entries and return the count in noptions.
     * note: noptions starts at 1 not 0 because we getoptlong
     * requires a NULL option entry to terminate the option array.
     * While we are at it, count the number of chars required
     * to store private copies of all the longoption strings
     * and return that count in logoptlen.
     */
     tabindex=0;
     do
        {
        const char *longopts = table[tabindex]->longopts;
        longoptlen += (longopts?strlen(longopts):0) + 1;
        while (longopts)
            {
            noptions++;
            longopts = strchr(longopts+1,',');
            }
        }while(!(table[tabindex++]->flag&ARG_TERMINATOR));
    /*printf("%d long options consuming %d chars in total\n",noptions,longoptlen);*/


    /* allocate storage for return data structure as: */
    /* (struct longoptions) + (struct options)[noptions] + char[longoptlen] */
    nbytes = sizeof(struct longoptions)
           + sizeof(struct option)*noptions
           + longoptlen;
    result = (struct longoptions*)malloc(nbytes);
    if (result)
        {
        int option_index=0;
        char *store;

        result->getoptval=0;
        result->noptions = noptions;
        result->options = (struct option*)(result + 1);
        store = (char*)(result->options + noptions);

        for(tabindex=0; !(table[tabindex]->flag&ARG_TERMINATOR); tabindex++)
            {
            const char *longopts = table[tabindex]->longopts;

            while(longopts && *longopts)
                {
                char *storestart = store;

                /* copy progressive longopt strings into the store */
                while (*longopts!=0 && *longopts!=',')
                    *store++ = *longopts++;
                *store++ = 0;
                if (*longopts==',')
                    longopts++;
                /*fprintf(stderr,"storestart=\"%s\"\n",storestart);*/

                result->options[option_index].name    = storestart;
                result->options[option_index].has_arg = (table[tabindex]->flag&ARG_HASVALUE) ? 1 : 0;
                result->options[option_index].flag    = &(result->getoptval);
                result->options[option_index].val     = tabindex;
                option_index++;
                }
            }
        /* terminate the options array with a zero-filled entry */
        result->options[option_index].name    = 0;
        result->options[option_index].has_arg = 0;
        result->options[option_index].flag    = 0;
        result->options[option_index].val     = 0;
        }

    /*dump_longoptions(result);*/
    return result;
    }

static
char* alloc_shortoptions(struct arg_hdr **table)
   {
   char *result;
   size_t len = 2;
   int tabindex;

   /* determine the total number of option chars required */
   for(tabindex=0; !(table[tabindex]->flag&ARG_TERMINATOR); tabindex++)
       {
       struct arg_hdr *hdr = table[tabindex];
       len += 2 * (hdr->shortopts?strlen(hdr->shortopts):0);
       }

   result = (char *)malloc(len);
   if (result)
        {
        char *res = result;

        /* add a leading ':' so getopt return codes distinguish    */
        /* unrecognised option and options missing argument values */
        *res++=':';

        for(tabindex=0; !(table[tabindex]->flag&ARG_TERMINATOR); tabindex++)
            {
            struct arg_hdr *hdr = table[tabindex];
            const char *shortopts = hdr->shortopts;
            while(shortopts && *shortopts)
                {
                *res++ = *shortopts++;
                if (hdr->flag & ARG_HASVALUE)
                    *res++=':';
                }
            }
        /* null terminate the string */
        *res=0;
        }

   /*printf("alloc_shortoptions() returns \"%s\"\n",(result?result:"NULL"));*/
   return result;
   }


/* return index of the table terminator entry */
static
int arg_endindex(struct arg_hdr **table)
    {
    int tabindex=0;
    while (!(table[tabindex]->flag&ARG_TERMINATOR))
        tabindex++;
    return tabindex;
    }


static
void arg_parse_tagged(int argc, char **argv, struct arg_hdr **table, struct arg_end *endtable)
    {
    struct longoptions *longoptions;
    char *shortoptions;
    int copt;

    /*printf("arg_parse_tagged(%d,%p,%p,%p)\n",argc,argv,table,endtable);*/

    /* allocate short and long option arrays for the given opttable[].   */
    /* if the allocs fail then put an error msg in the last table entry. */
    longoptions  = alloc_longoptions(table);
    shortoptions = alloc_shortoptions(table);
    if (!longoptions || !shortoptions)
        {
        /* one or both memory allocs failed */
        arg_register_error(endtable,endtable,ARG_EMALLOC,NULL);
        /* free anything that was allocated (this is null safe) */
        free(shortoptions);
        free(longoptions);
        return;
        }

    /*dump_longoptions(longoptions);*/

    /* reset getopts internal option-index to zero, and disable error reporting */
    optind = 0;
    opterr = 0;

    /* fetch and process args using getopt_long */
    while( (copt=getopt_long(argc,argv,shortoptions,longoptions->options,NULL)) != -1)
        {
        /*
        printf("optarg=%s\n",optarg);
        printf("optind=%d\n",optind);
        printf("copt=%c\n",(char)copt);
        printf("optopt=%c\n",optopt);
        */
        switch(copt)
            {
            case 0:
                {
                int tabindex = longoptions->getoptval;
                void *parent  = table[tabindex]->parent;
                /*printf("long option detected from argtable[%d]\n", tabindex);*/
                if (table[tabindex]->scanfn)
                    {
                    int errorcode = table[tabindex]->scanfn(parent,optarg);
                    if (errorcode!=0)
                        arg_register_error(endtable,parent,errorcode,optarg);
                    }
                }
                break;

            case '?':
                /*
                * getarg_long() found an unrecognised short option.
                * if it was a short option its value is in optopt
                * if it was a long option then optopt=0
                */
                switch (optopt)
                    {
                    case 0:
                        /*printf("?0 unrecognised long option %s\n",argv[optind-1]);*/
                        arg_register_error(endtable,endtable,ARG_ELONGOPT,argv[optind-1]);
                        break;
                    default:
                        /*printf("?* unrecognised short option '%c'\n",optopt);*/
                        arg_register_error(endtable,endtable,optopt,NULL);
                        break;
                    }
                break;

            case':':
                /*
                * getarg_long() found an option with its argument missing
                */
                /*printf(": option %s requires an argument\n",argv[optind-1]);*/
                arg_register_error(endtable,endtable,ARG_EMISSARG,argv[optind-1]);
                break;

            default:
                {
                /* getopt_long() found a valid short option */
                int tabindex = find_shortoption(table,(char)copt);
                /*printf("short option detected from argtable[%d]\n", tabindex);*/
                if (tabindex==-1)
                    {
                    /* should never get here - but handle it just in case */
                    /*printf("unrecognised short option %d\n",copt);*/
                    arg_register_error(endtable,endtable,copt,NULL);
                    }
                else
                    {
                    if (table[tabindex]->scanfn)
                        {
					
                        void *parent  = table[tabindex]->parent;
                        int errorcode = table[tabindex]->scanfn(parent,optarg);
                        if (errorcode!=0)
                            arg_register_error(endtable,parent,errorcode,optarg);
                        }
                    }
                break;
                }
            }
        }

    free(shortoptions);
    free(longoptions);
    }


static
void arg_parse_untagged(int argc, char **argv, struct arg_hdr **table, struct arg_end *endtable)
    {
    int tabindex=0;

    /*printf("arg_parse_untagged(%d,%p,%p,%p)\n",argc,argv,table,endtable);*/
    while (!(table[tabindex]->flag&ARG_TERMINATOR))
        {
        void *parent;
        int errorcode;

        /* if we have exhausted our argv[optind] entries then we have finished */
        if (optind>=argc)
            {
            /*printf("arg_parse_untagged(): argv[] exhausted\n");*/
            return;
            }

        /* skip table entries with non-null long or short options (they are not untagged entries) */
        if (table[tabindex]->longopts || table[tabindex]->shortopts)
            {
            /*printf("arg_parse_untagged(): skipping argtable[%d] (tagged argument)\n",tabindex);*/
            tabindex++;
            continue;
            }

        /* skip table entries with NULL scanfn */
        if (!(table[tabindex]->scanfn))
            {
            /*printf("arg_parse_untagged(): skipping argtable[%d] (NULL scanfn)\n",tabindex);*/
            tabindex++;
            continue;
            }

        /* attempt to scan the current argv[optind] with the current     */
        /* table[tabindex] entry. If it succeeds then keep it, otherwise */
        /* try again with the next table[] entry.                        */
        parent = table[tabindex]->parent;
        errorcode = table[tabindex]->scanfn(parent,argv[optind]);
        if (errorcode==0)
            {
            /* success, move onto next argv[optind] but stay with same table[tabindex] */
            /*printf("arg_parse_untagged(): argtable[%d] successfully matched\n",tabindex);*/
            optind++;
            }
        else
            {
            /* failure, try same argv[optind] with next table[tabindex] entry */
            /*printf("arg_parse_untagged(): argtable[%d] failed match\n",tabindex);*/
            tabindex++;
            }

        }

    /* only get here when not all argv[] entries were consumed */
    /* register an error for each unused argv[] entry */
    while (optind<argc)
        {
        /*printf("arg_parse_untagged(): argv[%d]=\"%s\" not consumed\n",optind,argv[optind]);*/
        arg_register_error(endtable,endtable,ARG_ENOMATCH,argv[optind++]);
        }

    return;
    }


static
void arg_parse_check(struct arg_hdr **table, struct arg_end *endtable)
    {
    int tabindex=0;
    do
        {
        if (table[tabindex]->checkfn)
            {
            void *parent  = table[tabindex]->parent;
            int errorcode = table[tabindex]->checkfn(parent);
            if (errorcode!=0)
                arg_register_error(endtable,parent,errorcode,NULL);
            }
        }while(!(table[tabindex++]->flag&ARG_TERMINATOR));
    }


static
void arg_reset(void **argtable)
    {
    struct arg_hdr **table=(struct arg_hdr**)argtable;
    int tabindex=0;
    /*printf("arg_reset(%p)\n",argtable);*/
    do
        {
        if (table[tabindex]->resetfn)
            table[tabindex]->resetfn(table[tabindex]->parent);
        } while(!(table[tabindex++]->flag&ARG_TERMINATOR));
    }


    
int arg_parse(int argc, char **argv, void **argtable)
    {
    struct arg_hdr **table = (struct arg_hdr **)argtable;
    struct arg_end *endtable;
    int endindex;

    /*printf("arg_parse(%d,%p,%p)\n",argc,argv,argtable);*/

    /* reset any argtable data from previous invocations */
    arg_reset(argtable);

    /* locate the first end-of-table marker within the array */
    endindex = arg_endindex(table);
    endtable = (struct arg_end*)table[endindex];

    /* parse the command line for tagged options */
    arg_parse_tagged(argc,argv,table,endtable);

    /* parse the command line for untagged options */
    arg_parse_untagged(argc,argv,table,endtable);

    /* perform post-parse checks */
    arg_parse_check(table,endtable);

    return endtable->count;
    }


/*
 * Concatenate contents of src[] string onto *pdest[] string.
 * The *pdest pointer is altered to point to the end of the
 * target string and *pndest is decremented by the same number
 * of chars.
 * Does not append more than *pndest chars into *pdest[]
 * so as to prevent buffer overruns.
 * Its something like strncat() but more efficient for repeated
 * calls on the same destination string.
 * Example of use:
 *   char dest[30] = "good"
 *   size_t ndest = sizeof(dest);
 *   char *pdest = dest;
 *   arg_char(&pdest,"bye ",&ndest);
 *   arg_char(&pdest,"cruel ",&ndest);
 *   arg_char(&pdest,"world!",&ndest);
 * Results in:
 *   dest[] == "goodbye cruel world!"
 *   ndest  == 10
 */
static
void arg_cat(char **pdest, const char *src, size_t *pndest)
    {
    char *dest = *pdest;
    char *end  = dest + *pndest;

    /*locate null terminator of dest string */
    while(dest<end && *dest!=0)
        dest++;

    /* concat src string to dest string */
    while(dest<end && *src!=0)
        *dest++ = *src++;

    /* null terminate dest string */
    *dest=0;

    /* update *pdest and *pndest */
    *pndest = end - dest;
    *pdest  = dest;
    }


static
void arg_cat_option(char *dest, size_t ndest, const char *shortopts, const char *longopts, const char *datatype)
    {
    if (shortopts)
        {
        char option[3] = { '-', shortopts[0], 0 };
        arg_cat(&dest,option,&ndest);
        if (datatype)
            {
            arg_cat(&dest," ",&ndest);
            arg_cat(&dest,datatype,&ndest);
            }
        }
    else if (longopts)
        {
        size_t ncspn;

        /* add "--" tag prefix */
        arg_cat(&dest,"--",&ndest);

        /* add comma separated option tag */
        ncspn = strcspn(longopts,",");
        strncat(dest,longopts,(ncspn<ndest)?ncspn:ndest);

        if (datatype)
            {
            arg_cat(&dest,"=",&ndest);
            arg_cat(&dest,datatype,&ndest);
            }
        }
    else if (datatype)
        arg_cat(&dest,datatype,&ndest);
    }

static
void arg_cat_optionv(char *dest, size_t ndest, const char *shortopts, const char *longopts, const char *datatype, const char *separator)
    {
    separator = separator ? separator : "";

    if (shortopts)
        {
        const char *c = shortopts;
        while(*c)
            {
            /* "-a|-b|-c" */
            char shortopt[]={'-',*c,0};
            arg_cat(&dest,shortopt,&ndest);
            if (*++c)
                arg_cat(&dest,separator,&ndest);
            }
        }

    /* put separator between long opts and short opts */
    if (shortopts && longopts)
        arg_cat(&dest,separator,&ndest);

    if (longopts)
        {
        const char *c = longopts;
        while(*c)
            {
            size_t ncspn;

            /* add "--" tag prefix */
            arg_cat(&dest,"--",&ndest);

            /* add comma separated option tag */
            ncspn = strcspn(c,",");
            strncat(dest,c,(ncspn<ndest)?ncspn:ndest);
            c+=ncspn;

            /* add given separator in place of comma */
            if (*c==',')
                 {
                 arg_cat(&dest,separator,&ndest);
                 c++;
                 }
            }
        }

    if (datatype)
        {
        if (longopts)
            arg_cat(&dest,"=",&ndest);
        else if (shortopts)
            arg_cat(&dest," ",&ndest);
        arg_cat(&dest,datatype,&ndest);
        }
    }


void arg_print_option(FILE *fp, const char *shortopts, const char *longopts, const char *datatype, const char *suffix)
    {
    char syntax[200]="";
    suffix = suffix ? suffix : "";
    arg_cat_optionv(syntax,sizeof(syntax),shortopts,longopts,datatype,"|");
    fputs(syntax,fp);
    fputs(suffix,fp);
    }


/*
 * Print a GNU style [OPTION] string in which all short options that
 * do not take argument values are presented in abbreviated form, as
 * in: -xvfsd, or -xvf[sd], or [-xvsfd]
 */
static
void arg_print_gnuswitch(struct arg_hdr **table)
    {
    int tabindex;
    const char *format1="-%c";
    const char *format2="[-%c";
    const char *suffix="";

    /* print all mandatory switches that are without argument values */
    for(tabindex=0; table[tabindex] && !(table[tabindex]->flag&ARG_TERMINATOR); tabindex++)
        {
        /* skip optional options */
        if (table[tabindex]->mincount<1)
            continue;

        /* skip non-short options */
        if (table[tabindex]->shortopts==NULL)
            continue;

        /* skip options that take argument values */
        if (table[tabindex]->flag&ARG_HASVALUE)
            continue;

        /* print the short option (only the first short option char, ignore multiple choices)*/
        printf(format1,table[tabindex]->shortopts[0]);
        format1="%c";
        format2="[%c";
        }

    /* print all optional switches that are without argument values */
    for(tabindex=0; table[tabindex] && !(table[tabindex]->flag&ARG_TERMINATOR); tabindex++)
        {
        /* skip mandatory args */
        if (table[tabindex]->mincount>0)
            continue;

        /* skip args without short options */
        if (table[tabindex]->shortopts==NULL)
            continue;

        /* skip args with values */
        if (table[tabindex]->flag&ARG_HASVALUE)
            continue;

        /* print first short option */
        printf(format2,table[tabindex]->shortopts[0]);
        format2="%c";
        suffix="]";
        }

    printf("%s",suffix);
    }


void arg_print_syntax(FILE *fp, void **argtable, const char *suffix)
    {
    struct arg_hdr **table = (struct arg_hdr**)argtable;
    int i,tabindex;

    /* print GNU style [OPTION] string */
    arg_print_gnuswitch(table);

    /* print remaining options in abbreviated style */
    for(tabindex=0; table[tabindex] && !(table[tabindex]->flag&ARG_TERMINATOR); tabindex++)
        {
        char syntax[200]="";
        const char *shortopts, *longopts, *datatype;

        /* skip short options without arg values (they were printed by arg_print_gnu_switch) */
        if (table[tabindex]->shortopts && !(table[tabindex]->flag&ARG_HASVALUE))
            continue;

        shortopts = table[tabindex]->shortopts;
        longopts  = table[tabindex]->longopts;
        datatype  = table[tabindex]->datatype;
        arg_cat_option(syntax,sizeof(syntax),shortopts,longopts,datatype);

        if (strlen(syntax)>0)
            {
            /* print mandatory options */
            for (i=0; i<table[tabindex]->mincount; i++)
                printf(" %s",syntax);

            /* print optional args enclosed in "[..]" */
            switch ( table[tabindex]->maxcount - table[tabindex]->mincount )
                {
                case 0:
                    break;
                case 1:
                    printf(" [%s]",syntax);
                    break;
                case 2:
                    printf(" [%s] [%s]",syntax,syntax);
                    break;
                default:
                    printf(" [%s]...",syntax);
                    break;
                }
            }
        }

    if (suffix)
        printf("%s",suffix);
    }


void arg_print_syntaxv(FILE *fp, void **argtable, const char *suffix)
    {
    struct arg_hdr **table = (struct arg_hdr**)argtable;
    int i,tabindex;

    /* print remaining options in abbreviated style */
    for(tabindex=0; table[tabindex] && !(table[tabindex]->flag&ARG_TERMINATOR); tabindex++)
        {
        char syntax[200]="";
        const char *shortopts, *longopts, *datatype;

        shortopts = table[tabindex]->shortopts;
        longopts  = table[tabindex]->longopts;
        datatype  = table[tabindex]->datatype;
        arg_cat_optionv(syntax,sizeof(syntax),shortopts,longopts,datatype,"|");

        /* print mandatory options */
        for (i=0; i<table[tabindex]->mincount; i++)
            printf(" %s",syntax);

        /* print optional args enclosed in "[..]" */
        switch ( table[tabindex]->maxcount - table[tabindex]->mincount )
            {
            case 0:
                break;
            case 1:
                printf(" [%s]",syntax);
                break;
            case 2:
                printf(" [%s] [%s]",syntax,syntax);
                break;
            default:
                printf(" [%s]...",syntax);
                break;
            }
        }

    if (suffix)
        printf("%s",suffix);
    }


void arg_print_glossary(FILE *fp, void **argtable, const char *format)
    {
    struct arg_hdr **table = (struct arg_hdr**)argtable;
    int tabindex;

    format = format ? format : "  %-20s %s\n";
    for(tabindex=0; !(table[tabindex]->flag&ARG_TERMINATOR); tabindex++)
        {
        if (table[tabindex]->glossary)
            {
            char syntax[200]="";
            const char *shortopts = table[tabindex]->shortopts;
            const char *longopts  = table[tabindex]->longopts;
            const char *datatype  = table[tabindex]->datatype;
            const char *glossary  = table[tabindex]->glossary;
            arg_cat_optionv(syntax,sizeof(syntax),shortopts,longopts,datatype,", ");
            fprintf(fp,format,syntax,glossary);
            }
        }
    }


/**
 * Checks the argtable[] array for NULL entries and returns 1
 * if any are found, zero otherwise.
 */
int arg_nullcheck(void **argtable)
    {
    struct arg_hdr **table = (struct arg_hdr **)argtable;
    int tabindex;
    /*printf("arg_nullcheck(%p)\n",argtable);*/

    if (!table)
        return 1;

    for (tabindex=0; !(table[tabindex]->flag&ARG_TERMINATOR); tabindex++)
        {
        if (!table[tabindex])
            return 1;
        }

    return 0;
    }



void arg_free(void **argtable)
    {
    struct arg_hdr **table=(struct arg_hdr**)argtable;
    int tabindex=0;
    int flag;
    /*printf("arg_free(%p)\n",argtable);*/
    do
        {
        flag = table[tabindex]->flag;
        free(table[tabindex]);
        table[tabindex++]=NULL;
        } while(!(flag&ARG_TERMINATOR));
    }




static void resetfndbl(struct arg_dbl *parent)
    {
    /*printf("%s:resetfn(%p)\n",__FILE__,parent);*/
    parent->count=0;
    }

static int scanfndbl(struct arg_dbl *parent, const char *argval)
    {
    int errorcode = 0;

    if (parent->count < parent->hdr.maxcount )
        {
        double val;
        char *end;

        /* extract double from argval into val */
        val = strtod(argval,&end);

        /* if success then store result in parent->dval[] array otherwise return error*/
        if (*end==0)
            parent->dval[parent->count++] = val;
        else
            errorcode = EBADDOUBLE;
        }
    else
        errorcode = EMAXCOUNT;

    /*printf("%s:scanfn(%p) returns %d\n",__FILE__,parent,errorcode);*/
    return errorcode;
    }

static int checkfndbl(struct arg_dbl *parent)
    {
    int errorcode = (parent->count < parent->hdr.mincount) ? EMINCOUNT : 0;
    /*printf("%s:checkfn(%p) returns %d\n",__FILE__,parent,errorcode);*/
    return errorcode;
    }

static void errorfndbl(struct arg_dbl *parent, FILE *fp, int errorcode, const char *argval, const char *progname)
    {
    const char *shortopts = parent->hdr.shortopts;
    const char *longopts  = parent->hdr.longopts;
    const char *datatype  = parent->hdr.datatype;

    fprintf(fp,"%s: ",progname);
    switch(errorcode)
        {
        case EMINCOUNT:
            fputs("missing option ",fp);
            arg_print_option(fp,shortopts,longopts,datatype,"\n");
            break;

        case EMAXCOUNT:
            fputs("excess option ",fp);
            arg_print_option(fp,shortopts,longopts,argval,"\n");
            break;

        case EBADDOUBLE:
            fprintf(fp,"invalid argument \"%s\" to option ",argval);
            arg_print_option(fp,shortopts,longopts,datatype,"\n");
            break;
        }
    }


struct arg_dbl* arg_dbl0(const char* shortopts,
                               const char* longopts,
                               const char *datatype,
                               const char *glossary)
    {
    return arg_dbln(shortopts,longopts,datatype,0,1,glossary);
    }

struct arg_dbl* arg_dbl1(const char* shortopts,
                               const char* longopts,
                               const char *datatype,
                               const char *glossary)
    {
    return arg_dbln(shortopts,longopts,datatype,1,1,glossary);
    }


struct arg_dbl* arg_dbln(const char* shortopts,
                               const char* longopts,
                               const char *datatype,
                               int mincount,
                               int maxcount,
                               const char *glossary)
    {
    size_t nbytes;
    struct arg_dbl *result;

    nbytes = sizeof(struct arg_dbl)     /* storage for struct arg_dbl */
           + maxcount * sizeof(double);    /* storage for dval[maxcount] array */

    result = (struct arg_dbl*)malloc(nbytes);
    if (result)
        {
        /* init the arg_hdr struct */
        result->hdr.flag      = ARG_HASVALUE;
        result->hdr.shortopts = shortopts;
        result->hdr.longopts  = longopts;
        result->hdr.datatype  = datatype ? datatype : "<double>";
        result->hdr.glossary  = glossary;
        result->hdr.mincount  = mincount;
        result->hdr.maxcount  = maxcount;
        result->hdr.parent    = result;
        result->hdr.resetfn   = (arg_resetfn*)resetfndbl;
        result->hdr.scanfn    = (arg_scanfn*)scanfndbl;
        result->hdr.checkfn   = (arg_checkfn*)checkfndbl;
        result->hdr.errorfn   = (arg_errorfn*)errorfndbl;

        /* store the dval[maxcount] array immediately after the arg_dbl struct */
        result->dval  = (double*)(result+1);
        result->count = 0;
        }
    /*printf("arg_dbln() returns %p\n",result);*/
    return result;
    }


static void resetfn(struct arg_end *parent)
    {
    /*printf("%s:resetfn(%p)\n",__FILE__,parent);*/
    parent->count = 0;
    }

static void errorfn(void *parent, FILE *fp, int error, const char *argval, const char *progname)
    {
    progname = progname ? progname : "";
    argval = argval ? argval : "";

    fprintf(fp,"%s: ",progname);
    switch(error)
        {
        case ARG_ELIMIT:
            fputs("too many errors to display",fp);
            break;
        case ARG_EMALLOC:
            fputs("insufficent memory",fp);
            break;
        case ARG_ENOMATCH:
            fprintf(fp,"excess parameter \"%s\"",argval);
            break;
        case ARG_EMISSARG:
            fprintf(fp,"option \"%s\" requires an argument",argval);
            break;
        case ARG_ELONGOPT:
            fprintf(fp,"invalid option \"%s\"",argval);
            break;
        default:
            fprintf(fp,"invalid option \"-%c\"",error);
            break;
        }
    fputc('\n',fp);
    }


struct arg_end* arg_end(int maxcount)
    {
    size_t nbytes;
    struct arg_end *result;

    nbytes = sizeof(struct arg_end)
           + maxcount * sizeof(int)             /* storage for int error[maxcount] array*/
           + maxcount * sizeof(void*)           /* storage for void* parent[maxcount] array */
           + maxcount * sizeof(char*);          /* storage for char* argval[maxcount] array */

    result = (struct arg_end*)malloc(nbytes);
    if (result)
        {
		memset(result,0,nbytes);

        /* init the arg_hdr struct */
        result->hdr.flag      = ARG_TERMINATOR;
        result->hdr.shortopts = NULL;
        result->hdr.longopts  = NULL;
        result->hdr.datatype  = NULL;
        result->hdr.glossary  = NULL;
        result->hdr.mincount  = 1;
        result->hdr.maxcount  = maxcount;
        result->hdr.parent    = result;
        result->hdr.resetfn   = (arg_resetfn*)resetfn;
        result->hdr.scanfn    = NULL;
        result->hdr.checkfn   = NULL;
        result->hdr.errorfn   = errorfn;

        /* store error[maxcount] array immediately after struct arg_end */
        result->error = (int*)(result+1);

        /* store parent[maxcount] array immediately after error[] array */
        result->parent = (void**)(result->error + maxcount + 1);

        /* store argval[maxcount] array immediately after parent[] array */
        result->argval = (const char**)(result->parent + maxcount + 1);
        }

    /*printf("arg_end(%d) returns %p\n",maxcount,result);*/
    return result;
    }


void arg_print_errors(FILE* fp, struct arg_end* end, const char* progname)
    {
    int i;
    /*printf("arg_errors()\n");*/
    for (i=0; i<end->count; i++)
        {
        struct arg_hdr *errorparent = (struct arg_hdr *)(end->parent[i]);
        if (errorparent->errorfn)
            errorparent->errorfn(end->parent[i],fp,end->error[i],end->argval[i],progname);
        }
    }



#define FILESEPARATOR '\\'


static void resetfnfilen(struct arg_file *parent)
    {
    /*printf("%s:resetfn(%p)\n",__FILE__,parent);*/
    parent->count=0;
    }


/* Returns ptr to the base filename within *filename */
static const char* arg_basename(const char *filename)
    {
    const char *result = (filename ? strrchr(filename,FILESEPARATOR) : NULL);
    if (result)
        result++;
    else
        result = filename;
    return result;
    }


/* Returns ptr to the file extension within *filename */
static const char* arg_extension(const char *filename)
    {
    const char *result = (filename ? strrchr(filename,'.') : NULL);
    if (filename && !result)
        result = filename+strlen(filename);
    return result;
    }


static int scanfnfilen(struct arg_file *parent, const char *argval)
    {
    int errorcode = 0;

    if (parent->count < parent->hdr.maxcount )
        {
        parent->filename[parent->count]  = argval;
        parent->basename[parent->count]  = arg_basename(argval);
        parent->extension[parent->count] = arg_extension(argval);
        parent->count++;
        }
    else
        errorcode = EMAXCOUNT;

    /*printf("%s:scanfn(%p) returns %d\n",__FILE__,parent,errorcode);*/
    return errorcode;
    }


static int checkfnfilen(struct arg_file *parent)
    {
    int errorcode = (parent->count < parent->hdr.mincount) ? EMINCOUNT : 0;
    /*printf("%s:checkfn(%p) returns %d\n",__FILE__,parent,errorcode);*/
    return errorcode;
    }


static void errorfnfilen(struct arg_file *parent, FILE *fp, int errorcode, const char *argval, const char *progname)
    {
    const char *shortopts = parent->hdr.shortopts;
    const char *longopts  = parent->hdr.longopts;
    const char *datatype  = parent->hdr.datatype;

    fprintf(fp,"%s: ",progname);
    switch(errorcode)
        {
        case EMINCOUNT:
            fputs("missing option ",fp);
            arg_print_option(fp,shortopts,longopts,datatype,"\n");
            break;

        case EMAXCOUNT:
            fputs("excess option ",fp);
            arg_print_option(fp,shortopts,longopts,argval,"\n");
            break;

        default:
            fprintf(fp,"unknown error at \"%s\"\n",argval);
        }
    }


struct arg_file* arg_file0(const char* shortopts,
                           const char* longopts,
                           const char *datatype,
                           const char *glossary)
    {
    return arg_filen(shortopts,longopts,datatype,0,1,glossary);
    }


struct arg_file* arg_file1(const char* shortopts,
                           const char* longopts,
                           const char *datatype,
                           const char *glossary)
    {
    return arg_filen(shortopts,longopts,datatype,1,1,glossary);
    }


struct arg_file* arg_filen(const char* shortopts,
                           const char* longopts,
                           const char *datatype,
                           int mincount,
                           int maxcount,
                           const char *glossary)
    {
    size_t nbytes;
    struct arg_file *result;

    nbytes = sizeof(struct arg_file)     /* storage for struct arg_file */
           + sizeof(char*) * maxcount    /* storage for filename[maxcount] array */
           + sizeof(char*) * maxcount    /* storage for basename[maxcount] array */
           + sizeof(char*) * maxcount;   /* storage for extension[maxcount] array */

    result = (struct arg_file*)malloc(nbytes);
    if (result)
        {
        /* init the arg_hdr struct */
        result->hdr.flag      = ARG_HASVALUE;
        result->hdr.shortopts = shortopts;
        result->hdr.longopts  = longopts;
        result->hdr.glossary  = glossary;
        result->hdr.datatype  = datatype ? datatype : "<file>";
        result->hdr.mincount  = mincount;
        result->hdr.maxcount  = maxcount;
        result->hdr.parent    = result;
        result->hdr.resetfn   = (arg_resetfn*)resetfnfilen;
        result->hdr.scanfn    = (arg_scanfn*)scanfnfilen;
        result->hdr.checkfn   = (arg_checkfn*)checkfnfilen;
        result->hdr.errorfn   = (arg_errorfn*)errorfnfilen;

        /* store the filename,basename,extension arrays immediately after the arg_file struct */
        result->filename  = (const char**)(result+1);
        result->basename  = result->filename + maxcount;
        result->extension = result->basename + maxcount;
        result->count = 0;
        }
    /*printf("arg_filen() returns %p\n",result);*/
    return result;
    }


static void resetfnintn(struct arg_int *parent)
    {
    /*printf("%s:resetfn(%p)\n",__FILE__,parent);*/
    parent->count=0;
    }

static int scanfnintn(struct arg_int *parent, const char *argval)
    {
    int errorcode = 0;

    if (parent->count < parent->hdr.maxcount )
        {
        int val;
        char *end;

        /* extract base10 integer from argval into val */
        val = (int)strtol(argval,&end,10);

        /* if success then store result in parent->ival[] array otherwise return error*/
        if (*end==0)
            parent->ival[parent->count++] = val;
        else
            errorcode = EBADINT;
        }
    else
        errorcode = EMAXCOUNT;

    /*printf("%s:scanfn(%p) returns %d\n",__FILE__,parent,errorcode);*/
    return errorcode;
    }

static int checkfnintn(struct arg_int *parent)
    {
    int errorcode = (parent->count < parent->hdr.mincount) ? EMINCOUNT : 0;
    /*printf("%s:checkfn(%p) returns %d\n",__FILE__,parent,errorcode);*/
    return errorcode;
    }

static void errorfnintn(struct arg_int *parent, FILE *fp, int errorcode, const char *argval, const char *progname)
    {
    const char *shortopts = parent->hdr.shortopts;
    const char *longopts  = parent->hdr.longopts;
    const char *datatype  = parent->hdr.datatype;

    fprintf(fp,"%s: ",progname);
    switch(errorcode)
        {
        case EMINCOUNT:
            fputs("missing option ",fp);
            arg_print_option(fp,shortopts,longopts,datatype,"\n");
            break;

        case EMAXCOUNT:
            fputs("excess option ",fp);
            arg_print_option(fp,shortopts,longopts,argval,"\n");
            break;

        case EBADINT:
            fprintf(fp,"invalid argument \"%s\" to option ",argval);
            arg_print_option(fp,shortopts,longopts,datatype,"\n");
            break;
        }
    }


struct arg_int* arg_int0(const char* shortopts,
                         const char* longopts,
                         const char *datatype,
                         const char *glossary)
    {
    return arg_intn(shortopts,longopts,datatype,0,1,glossary);
    }

struct arg_int* arg_int1(const char* shortopts,
                         const char* longopts,
                         const char *datatype,
                         const char *glossary)
    {
    return arg_intn(shortopts,longopts,datatype,1,1,glossary);
    }


struct arg_int* arg_intn(const char* shortopts,
                         const char* longopts,
                         const char *datatype,
                         int mincount,
                         int maxcount,
                         const char *glossary)
    {
    size_t nbytes;
    struct arg_int *result;

    nbytes = sizeof(struct arg_int)     /* storage for struct arg_int */
           + maxcount * sizeof(int);    /* storage for ival[maxcount] array */

    result = (struct arg_int*)malloc(nbytes);
    if (result)
        {
        /* init the arg_hdr struct */
        result->hdr.flag      = ARG_HASVALUE;
        result->hdr.shortopts = shortopts;
        result->hdr.longopts  = longopts;
        result->hdr.datatype  = datatype ? datatype : "<int>";
        result->hdr.glossary  = glossary;
        result->hdr.mincount  = mincount;
        result->hdr.maxcount  = maxcount;
        result->hdr.parent    = result;
        result->hdr.resetfn   = (arg_resetfn*)resetfnintn;
        result->hdr.scanfn    = (arg_scanfn*)scanfnintn;
        result->hdr.checkfn   = (arg_checkfn*)checkfnintn;
        result->hdr.errorfn   = (arg_errorfn*)errorfnintn;

        /* store the ival[maxcount] array immediately after the arg_int struct */
        result->ival  = (int*)(result+1);
        result->count = 0;
        }
    /*printf("arg_intn() returns %p\n",result);*/
    return result;
    }

static void resetfnlitn(struct arg_lit *parent)
    {
    /*printf("%s:resetfn(%p)\n",__FILE__,parent);*/
    parent->count = 0;
    }

static int scanfnlitn(struct arg_lit *parent, const char *argval)
    {
    int errorcode = 0;
    if (parent->count < parent->hdr.maxcount )
        parent->count++;
    else
        errorcode = EMAXCOUNT;
    /*printf("%s:scanfn(%p,%s) returns %d\n",__FILE__,parent,argval,errorcode);*/
    return errorcode;
    }

static int checkfnlitn(struct arg_lit *parent)
    {
    int errorcode = (parent->count < parent->hdr.mincount) ? EMINCOUNT : 0;
    /*printf("%s:checkfn(%p) returns %d\n",__FILE__,parent,errorcode);*/
    return errorcode;
    }

static void errorfnlitn(struct arg_lit *parent, FILE *fp, int errorcode, const char *argval, const char *progname)
    {
    const char *shortopts = parent->hdr.shortopts;
    const char *longopts  = parent->hdr.longopts;
    const char *datatype  = parent->hdr.datatype;

    switch(errorcode)
        {
        case EMINCOUNT:
            fprintf(fp,"%s: missing option ",progname);
            arg_print_option(fp,shortopts,longopts,datatype,"\n");
            fprintf(fp,"\n");
            break;

        case EMAXCOUNT:
            fprintf(fp,"%s: extraneous option ",progname);
            arg_print_option(fp,shortopts,longopts,datatype,"\n");
            break;
        }
    }

struct arg_lit* arg_lit0(const char* shortopts,
                         const char* longopts,
                         const char* glossary)
    {return arg_litn(shortopts,longopts,0,1,glossary);}

struct arg_lit* arg_lit1(const char* shortopts,
                         const char* longopts,
                         const char* glossary)
    {return arg_litn(shortopts,longopts,1,1,glossary);}


struct arg_lit* arg_litn(const char* shortopts,
                         const char* longopts,
                         int mincount,
                         int maxcount,
                         const char *glossary)
    {
    struct arg_lit *result = (struct arg_lit*)malloc(sizeof(struct arg_lit));
    if (result)
        {
        /* init the arg_hdr struct */
        result->hdr.flag      = 0;
        result->hdr.shortopts = shortopts;
        result->hdr.longopts  = longopts;
        result->hdr.datatype  = NULL;
        result->hdr.glossary  = glossary;
        result->hdr.mincount  = mincount;
        result->hdr.maxcount  = maxcount;
        result->hdr.parent    = result;
        result->hdr.resetfn   = (arg_resetfn*)resetfnlitn;
        result->hdr.scanfn    = (arg_scanfn*)scanfnlitn;
        result->hdr.checkfn   = (arg_checkfn*)checkfnlitn;
        result->hdr.errorfn   = (arg_errorfn*)errorfnlitn;

        /* init local variables */
        result->count = 0;
        }
    /*printf("arg_litn() returns %p\n",result);*/
    return result;
    }


struct arg_rem* arg_rem(const char *datatype,
                        const char *glossary)
    {
    struct arg_rem *result = (struct arg_rem*)malloc(sizeof(struct arg_rem));
    if (result)
        {
        /* init the arg_hdr struct */
        result->hdr.flag      = 0;
        result->hdr.shortopts = NULL;
        result->hdr.longopts  = NULL;
        result->hdr.datatype  = datatype;
        result->hdr.glossary  = glossary;
        result->hdr.mincount  = 1;
        result->hdr.maxcount  = 1;
        result->hdr.parent    = result;
        result->hdr.resetfn   = NULL;
        result->hdr.scanfn    = NULL;
        result->hdr.checkfn   = NULL;
        result->hdr.errorfn   = NULL;
        }
    /*printf("arg_rem() returns %p\n",result);*/
    return result;
    }

static void resetfnstrn(struct arg_str *parent)
    {
    /*printf("%s:resetfn(%p)\n",__FILE__,parent);*/
    parent->count=0;
    }

static int scanfnstrn(struct arg_str *parent, const char *argval)
    {
    int errorcode = 0;

    if (parent->count < parent->hdr.maxcount )
        parent->sval[parent->count++] = argval;
    else
        errorcode = EMAXCOUNT;

    /*printf("%s:scanfn(%p) returns %d\n",__FILE__,parent,errorcode);*/
    return errorcode;
    }

static int checkfnstrn(struct arg_str *parent)
    {
    int errorcode = (parent->count < parent->hdr.mincount) ? EMINCOUNT : 0;
    /*printf("%s:checkfn(%p) returns %d\n",__FILE__,parent,errorcode);*/
    return errorcode;
    }

static void errorfnstrn(struct arg_str *parent, FILE *fp, int errorcode, const char *argval, const char *progname)
    {
    const char *shortopts = parent->hdr.shortopts;
    const char *longopts  = parent->hdr.longopts;
    const char *datatype  = parent->hdr.datatype;

    fprintf(fp,"%s: ",progname);
    switch(errorcode)
        {
        case EMINCOUNT:
            fputs("missing option ",fp);
            arg_print_option(fp,shortopts,longopts,datatype,"\n");
            break;

        case EMAXCOUNT:
            fputs("excess option ",fp);
            arg_print_option(fp,shortopts,longopts,argval,"\n");
            break;
        }
    }


struct arg_str* arg_str0(const char* shortopts,
                         const char* longopts,
                         const char *datatype,
                         const char *glossary)
    {
    return arg_strn(shortopts,longopts,datatype,0,1,glossary);
    }

struct arg_str* arg_str1(const char* shortopts,
                         const char* longopts,
                         const char *datatype,
                         const char *glossary)
    {
    return arg_strn(shortopts,longopts,datatype,1,1,glossary);
    }


struct arg_str* arg_strn(const char* shortopts,
                         const char* longopts,
                         const char *datatype,
                         int mincount,
                         int maxcount,
                         const char *glossary)
    {
    size_t nbytes;
    struct arg_str *result;

    nbytes = sizeof(struct arg_str)     /* storage for struct arg_str */
           + maxcount * sizeof(char*);  /* storage for sval[maxcount] array */

    result = (struct arg_str*)malloc(nbytes);
    if (result)
        {
        /* init the arg_hdr struct */
        result->hdr.flag      = ARG_HASVALUE;
        result->hdr.shortopts = shortopts;
        result->hdr.longopts  = longopts;
        result->hdr.datatype  = datatype ? datatype : "<string>";
        result->hdr.glossary  = glossary;
        result->hdr.mincount  = mincount;
        result->hdr.maxcount  = maxcount;
        result->hdr.parent    = result;
        result->hdr.resetfn   = (arg_resetfn*)resetfnstrn;
        result->hdr.scanfn    = (arg_scanfn*)scanfnstrn;
        result->hdr.checkfn   = (arg_checkfn*)checkfnstrn;
        result->hdr.errorfn   = (arg_errorfn*)errorfnstrn;

        /* store the sval[maxcount] array immediately after the arg_str struct */
        result->sval  = (const char**)(result+1);
        result->count = 0;
        }
    /*printf("arg_strn() returns %p\n",result);*/
    return result;
    }


/* Getopt for GNU.
   NOTE: getopt is now part of the C library, so if you don't know what
   "Keep this file name-space clean" means, talk to roland@gnu.ai.mit.edu
   before changing it!

   Copyright (C) 1987, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97
   Free Software Foundation, Inc.

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
/* This tells Alpha OSF/1 not to define a getopt prototype in <stdio.h>.
   Ditto for AIX 3.2 and <stdlib.h>.  */

/* This version of `getopt' appears to the caller like standard Unix `getopt'
   but it behaves differently for the user, since it allows the user
   to intersperse the options with the other arguments.

   As `getopt' works, it permutes the elements of ARGV so that,
   when it is done, all the options precede everything else.  Thus
   all application programs are extended to handle flexible argument order.

   Setting the environment variable POSIXLY_CORRECT disables permutation.
   Then the behavior is completely standard.

   GNU application programs can use a third alternative mode in which
   they can distinguish the relative order of options and other arguments.  */

/* For communication from `getopt' to the caller.
   When `getopt' finds an option that takes an argument,
   the argument value is returned here.
   Also, when `ordering' is RETURN_IN_ORDER,
   each non-option ARGV-element is returned here.  */

char *optarg = NULL;

/* Index in ARGV of the next element to be scanned.
   This is used for communication to and from the caller
   and for communication between successive calls to `getopt'.

   On entry to `getopt', zero means this is the first call; initialize.

   When `getopt' returns -1, this is the index of the first of the
   non-option elements that the caller should itself scan.

   Otherwise, `optind' communicates from one call to the next
   how much of ARGV has been scanned so far.  */

/* 1003.2 says this must be 1 before any call.  */
int optind = 1;

/* Formerly, initialization of getopt depended on optind==0, which
   causes problems with re-calling getopt as programs generally don't
   know that. */

int __getopt_initialized = 0;

/* The next char to be scanned in the option-element
   in which the last option character we returned was found.
   This allows us to pick up the scan where we left off.

   If this is zero, or a null string, it means resume the scan
   by advancing to the next ARGV-element.  */

static char *nextchar;

/* Callers store zero here to inhibit the error message
   for unrecognized options.  */

int opterr = 1;

/* Set to an option character which was unrecognized.
   This must be initialized on some systems to avoid linking in the
   system's own getopt implementation.  */

int optopt = '?';

/* Describe how to deal with options that follow non-option ARGV-elements.

   If the caller did not specify anything,
   the default is REQUIRE_ORDER if the environment variable
   POSIXLY_CORRECT is defined, PERMUTE otherwise.

   REQUIRE_ORDER means don't recognize them as options;
   stop option processing when the first non-option is seen.
   This is what Unix does.
   This mode of operation is selected by either setting the environment
   variable POSIXLY_CORRECT, or using `+' as the first character
   of the list of option characters.

   PERMUTE is the default.  We permute the contents of ARGV as we scan,
   so that eventually all the non-options are at the end.  This allows options
   to be given in any order, even with programs that were not written to
   expect this.

   RETURN_IN_ORDER is an option available to programs that were written
   to expect options and other ARGV-elements in any order and that care about
   the ordering of the two.  We describe each non-option ARGV-element
   as if it were the argument of an option with character code 1.
   Using `-' as the first character of the list of option characters
   selects this mode of operation.

   The special argument `--' forces an end of option-scanning regardless
   of the value of `ordering'.  In the case of RETURN_IN_ORDER, only
   `--' can cause `getopt' to return -1 with `optind' != ARGC.  */

static enum
{
	REQUIRE_ORDER, PERMUTE, RETURN_IN_ORDER
}
ordering;

/* Value of POSIXLY_CORRECT environment variable.  */
static char *posixly_correct = NULL;

#include <string.h>
#define	my_index	strchr

/* Handle permutation of arguments.  */

/* Describe the part of ARGV that contains non-options that have
   been skipped.  `first_nonopt' is the index in ARGV of the first of them;
   `last_nonopt' is the index after the last of them.  */

static int first_nonopt;
static int last_nonopt;

/* Exchange two adjacent subsequences of ARGV.
   One subsequence is elements [first_nonopt,last_nonopt)
   which contains all the non-options that have been skipped so far.
   The other is elements [last_nonopt,optind), which contains all
   the options processed since those non-options were skipped.

   `first_nonopt' and `last_nonopt' are relocated so that they describe
   the new indices of the non-options in ARGV after they are moved.  */


static void
 exchange(char **argv)
{
	int bottom = first_nonopt;
	int middle = last_nonopt;
	int top = optind;
	char *tem;

	/* Exchange the shorter segment with the far end of the longer segment.
	   That puts the shorter segment into the right place.
	   It leaves the longer segment in the right place overall,
	   but it consists of two parts that need to be swapped next.  */

	while (top > middle && middle > bottom)
	{
		if (top - middle > middle - bottom)
		{
			/* Bottom segment is the short one.  */
			int len = middle - bottom;
			register int i;

			/* Swap it with the top part of the top segment.  */
			for (i = 0; i < len; i++)
			{
				tem = argv[bottom + i];
				argv[bottom + i] = argv[top - (middle - bottom) + i];
				argv[top - (middle - bottom) + i] = tem;
			}
			/* Exclude the moved bottom segment from further swapping.  */
			top -= len;
		}
		else
		{
			/* Top segment is the short one.  */
			int len = top - middle;
			register int i;

			/* Swap it with the bottom part of the bottom segment.  */
			for (i = 0; i < len; i++)
			{
				tem = argv[bottom + i];
				argv[bottom + i] = argv[middle + i];
				argv[middle + i] = tem;
			}
			/* Exclude the moved top segment from further swapping.  */
			bottom += len;
		}
	}

	/* Update records for the slots the non-options now occupy.  */

	first_nonopt += (optind - last_nonopt);
	last_nonopt = optind;
}

/* Initialize the internal data when the first call is made.  */

static const char *
     _getopt_initialize(int argc,  char *const *argv, const char *optstring)
{
	/* Start processing options with ARGV-element 1 (since ARGV-element 0
	   is the program name); the sequence of previously skipped
	   non-option ARGV-elements is empty.  */

	first_nonopt = last_nonopt = optind = 1;

	nextchar = NULL;

	/* Determine how to handle the ordering of options and nonoptions.  */

	if (optstring[0] == '-')
	{
		ordering = RETURN_IN_ORDER;
		++optstring;
	}
	else if (optstring[0] == '+')
	{
		ordering = REQUIRE_ORDER;
		++optstring;
	}
	else if (posixly_correct != NULL)
		ordering = REQUIRE_ORDER;
	else
		ordering = PERMUTE;

	return optstring;
}

/* Scan elements of ARGV (whose length is ARGC) for option characters
   given in OPTSTRING.

   If an element of ARGV starts with '-', and is not exactly "-" or "--",
   then it is an option element.  The characters of this element
   (aside from the initial '-') are option characters.  If `getopt'
   is called repeatedly, it returns successively each of the option characters
   from each of the option elements.

   If `getopt' finds another option character, it returns that character,
   updating `optind' and `nextchar' so that the next call to `getopt' can
   resume the scan with the following option character or ARGV-element.

   If there are no more option characters, `getopt' returns -1.
   Then `optind' is the index in ARGV of the first ARGV-element
   that is not an option.  (The ARGV-elements have been permuted
   so that those that are not options now come last.)

   OPTSTRING is a string containing the legitimate option characters.
   If an option character is seen that is not listed in OPTSTRING,
   return '?' after printing an error message.  If you set `opterr' to
   zero, the error message is suppressed but we still return '?'.

   If a char in OPTSTRING is followed by a colon, that means it wants an arg,
   so the following text in the same ARGV-element, or the text of the following
   ARGV-element, is returned in `optarg'.  Two colons mean an option that
   wants an optional arg; if there is text in the current ARGV-element,
   it is returned in `optarg', otherwise `optarg' is set to zero.

   If OPTSTRING starts with `-' or `+', it requests different methods of
   handling the non-option ARGV-elements.
   See the comments about RETURN_IN_ORDER and REQUIRE_ORDER, above.

   Long-named options begin with `--' instead of `-'.
   Their names may be abbreviated as long as the abbreviation is unique
   or is an exact match for some defined option.  If they have an
   argument, it follows the option name in the same ARGV-element, separated
   from the option name by a `=', or else the in next ARGV-element.
   When `getopt' finds a long-named option, it returns 0 if that option's
   `flag' field is nonzero, the value of the option's `val' field
   if the `flag' field is zero.

   The elements of ARGV aren't really const, because we permute them.
   But we pretend they're const in the prototype to be compatible
   with other systems.

   LONGOPTS is a vector of `struct option' terminated by an
   element containing a name which is zero.

   LONGIND returns the index in LONGOPT of the long-named option found.
   It is only valid when a long-named option has been found by the most
   recent call.

   If LONG_ONLY is nonzero, '-' as well as '--' can introduce
   long-named options.  */

int
    _getopt_internal(int argc,  char *const *argv,  const char *optstring,  const struct option *longopts,
     int *longind, int long_only)
{
	optarg = NULL;

	if (!__getopt_initialized || optind == 0)
	{
		optstring = _getopt_initialize(argc, argv, optstring);
		optind = 1;	/* Don't scan ARGV[0], the program name.  */
		__getopt_initialized = 1;
	}

	/* Test whether ARGV[optind] points to a non-option argument.
	   Either it does not have option syntax, or there is an environment flag
	   from the shell indicating it is not an option.  The later information
	   is only used when the used in the GNU libc.  */

#define NONOPTION_P (argv[optind][0] != '-' || argv[optind][1] == '\0')

	if (nextchar == NULL || *nextchar == '\0')
	{
		/* Advance to the next ARGV-element.  */

		/* Give FIRST_NONOPT & LAST_NONOPT rational values if OPTIND has been
		   moved back by the user (who may also have changed the arguments).  */
		if (last_nonopt > optind)
			last_nonopt = optind;
		if (first_nonopt > optind)
			first_nonopt = optind;

		if (ordering == PERMUTE)
		{
			/* If we have just processed some options following some non-options,
			   exchange them so that the options come first.  */

			if (first_nonopt != last_nonopt && last_nonopt != optind)
				exchange((char **) argv);
			else if (last_nonopt != optind)
				first_nonopt = optind;

			/* Skip any additional non-options
			   and extend the range of non-options previously skipped.  */

			while (optind < argc && NONOPTION_P)
				optind++;
			last_nonopt = optind;
		}

		/* The special ARGV-element `--' means premature end of options.
		   Skip it like a null option,
		   then exchange with previous non-options as if it were an option,
		   then skip everything else like a non-option.  */

		if (optind != argc && !strcmp(argv[optind], "--"))
		{
			optind++;

			if (first_nonopt != last_nonopt && last_nonopt != optind)
				exchange((char **) argv);
			else if (first_nonopt == last_nonopt)
				first_nonopt = optind;
			last_nonopt = argc;

			optind = argc;
		}

		/* If we have done all the ARGV-elements, stop the scan
		   and back over any non-options that we skipped and permuted.  */

		if (optind == argc)
		{
			/* Set the next-arg-index to point at the non-options
			   that we previously skipped, so the caller will digest them.  */
			if (first_nonopt != last_nonopt)
				optind = first_nonopt;
			return -1;
		}

		/* If we have come to a non-option and did not permute it,
		   either stop the scan or describe it to the caller and pass it by.  */

		if (NONOPTION_P)
		{
			if (ordering == REQUIRE_ORDER)
				return -1;
			optarg = argv[optind++];
			return 1;
		}

		/* We have found another option-ARGV-element.
		   Skip the initial punctuation.  */

		nextchar = (argv[optind] + 1
			    + (longopts != NULL && argv[optind][1] == '-'));
	}

	/* Decode the current option-ARGV-element.  */

	/* Check whether the ARGV-element is a long option.

	   If long_only and the ARGV-element has the form "-f", where f is
	   a valid short option, don't consider it an abbreviated form of
	   a long option that starts with f.  Otherwise there would be no
	   way to give the -f short option.

	   On the other hand, if there's a long option "fubar" and
	   the ARGV-element is "-fu", do consider that an abbreviation of
	   the long option, just like "--fu", and not "-f" with arg "u".

	   This distinction seems to be the most useful approach.  */

	if (longopts != NULL
	    && (argv[optind][1] == '-'
		|| (long_only && (argv[optind][2] || !my_index(optstring, argv[optind][1])))))
	{
		char *nameend;
		const struct option *p;
		const struct option *pfound = NULL;
		int exact = 0;
		int ambig = 0;
		int indfound = -1;
		int option_index;

		for (nameend = nextchar; *nameend && *nameend != '='; nameend++)
			/* Do nothing.  */ ;

		/* Test all long options for either exact match
		   or abbreviated matches.  */
		for (p = longopts, option_index = 0; p->name; p++, option_index++)
			if (!strncmp(p->name, nextchar, nameend - nextchar))
			{
				if ((unsigned int) (nameend - nextchar)
				    == (unsigned int) strlen(p->name))
				{
					/* Exact match found.  */
					pfound = p;
					indfound = option_index;
					exact = 1;
					break;
				}
				else if (pfound == NULL)
				{
					/* First nonexact match found.  */
					pfound = p;
					indfound = option_index;
				}
				else
					/* Second or later nonexact match found.  */
					ambig = 1;
			}

		if (ambig && !exact)
		{
			if (opterr)
				fprintf(stderr, "%s: option `%s' is ambiguous\n",
					argv[0], argv[optind]);
			nextchar += strlen(nextchar);
			optind++;
			optopt = 0;
			return '?';
		}

		if (pfound != NULL)
		{
			option_index = indfound;
			optind++;
			if (*nameend)
			{
				/* Don't test has_arg with >, because some C compilers don't
				   allow it to be used on enums.  */
				if (pfound->has_arg)
					optarg = nameend + 1;
				else
				{
					if (opterr)
					{
						if (argv[optind - 1][1] == '-')
							/* --option */
							fprintf(stderr,
								"%s: option `--%s' doesn't allow an argument\n",
								argv[0], pfound->name);
						else
							/* +option or -option */
							fprintf(stderr,
								"%s: option `%c%s' doesn't allow an argument\n",
								argv[0], argv[optind - 1][0], pfound->name);
					}

					nextchar += strlen(nextchar);

					optopt = pfound->val;
					return '?';
				}
			}
			else if (pfound->has_arg == 1)
			{
				if (optind < argc)
					optarg = argv[optind++];
				else
				{
					if (opterr)
						fprintf(stderr,
							"%s: option `%s' requires an argument\n",
						 argv[0], argv[optind - 1]);
					nextchar += strlen(nextchar);
					optopt = pfound->val;
					return optstring[0] == ':' ? ':' : '?';
				}
			}
			nextchar += strlen(nextchar);
			if (longind != NULL)
				*longind = option_index;
			if (pfound->flag)
			{
				*(pfound->flag) = pfound->val;
				return 0;
			}
			return pfound->val;
		}

		/* Can't find it as a long option.  If this is not getopt_long_only,
		   or the option starts with '--' or is not a valid short
		   option, then it's an error.
		   Otherwise interpret it as a short option.  */
		if (!long_only || argv[optind][1] == '-'
		    || my_index(optstring, *nextchar) == NULL)
		{
			if (opterr)
			{
				if (argv[optind][1] == '-')
					/* --option */
					fprintf(stderr, "%s: unrecognized option `--%s'\n",
						argv[0], nextchar);
				else
					/* +option or -option */
					fprintf(stderr, "%s: unrecognized option `%c%s'\n",
					argv[0], argv[optind][0], nextchar);
			}
			nextchar = (char *) "";
			optind++;
			optopt = 0;
			return '?';
		}
	}

	/* Look at and handle the next short option-character.  */

	{
		char c = *nextchar++;
		char *temp = my_index((char *)optstring, c);

		/* Increment `optind' when we start to process its last character.  */
		if (*nextchar == '\0')
			++optind;

		if (temp == NULL || c == ':')
		{
			if (opterr)
			{
				if (posixly_correct)
					/* 1003.2 specifies the format of this message.  */
					fprintf(stderr, "%s: illegal option -- %c\n",
						argv[0], c);
				else
					fprintf(stderr, "%s: invalid option -- %c\n",
						argv[0], c);
			}
			optopt = c;
			return '?';
		}
		/* Convenience. Treat POSIX -W foo same as long option --foo */
		if (temp[0] == 'W' && temp[1] == ';')
		{
			char *nameend;
			const struct option *p;
			const struct option *pfound = NULL;
			int exact = 0;
			int ambig = 0;
			int indfound = 0;
			int option_index;

			/* This is an option that requires an argument.  */
			if (*nextchar != '\0')
			{
				optarg = nextchar;
				/* If we end this ARGV-element by taking the rest as an arg,
				   we must advance to the next element now.  */
				optind++;
			}
			else if (optind == argc)
			{
				if (opterr)
				{
					/* 1003.2 specifies the format of this message.  */
					fprintf(stderr, "%s: option requires an argument -- %c\n",
						argv[0], c);
				}
				optopt = c;
				if (optstring[0] == ':')
					c = ':';
				else
					c = '?';
				return c;
			}
			else
				/* We already incremented `optind' once;
				   increment it again when taking next ARGV-elt as argument.  */
				optarg = argv[optind++];

			/* optarg is now the argument, see if it's in the
			   table of longopts.  */

			for (nextchar = nameend = optarg; *nameend && *nameend != '='; nameend++)
				/* Do nothing.  */ ;

			/* Test all long options for either exact match
			   or abbreviated matches.  */
			for (p = longopts, option_index = 0; p->name; p++, option_index++)
				if (!strncmp(p->name, nextchar, nameend - nextchar))
				{
					if ((unsigned int) (nameend - nextchar) == strlen(p->name))
					{
						/* Exact match found.  */
						pfound = p;
						indfound = option_index;
						exact = 1;
						break;
					}
					else if (pfound == NULL)
					{
						/* First nonexact match found.  */
						pfound = p;
						indfound = option_index;
					}
					else
						/* Second or later nonexact match found.  */
						ambig = 1;
				}
			if (ambig && !exact)
			{
				if (opterr)
					fprintf(stderr, "%s: option `-W %s' is ambiguous\n",
						argv[0], argv[optind]);
				nextchar += strlen(nextchar);
				optind++;
				return '?';
			}
			if (pfound != NULL)
			{
				option_index = indfound;
				if (*nameend)
				{
					/* Don't test has_arg with >, because some C compilers don't
					   allow it to be used on enums.  */
					if (pfound->has_arg)
						optarg = nameend + 1;
					else
					{
						if (opterr)
							fprintf(stderr, "%s: option `-W %s' doesn't allow an argument\n",
								argv[0], pfound->name);

						nextchar += strlen(nextchar);
						return '?';
					}
				}
				else if (pfound->has_arg == 1)
				{
					if (optind < argc)
						optarg = argv[optind++];
					else
					{
						if (opterr)
							fprintf(stderr,
								"%s: option `%s' requires an argument\n",
								argv[0], argv[optind - 1]);
						nextchar += strlen(nextchar);
						return optstring[0] == ':' ? ':' : '?';
					}
				}
				nextchar += strlen(nextchar);
				if (longind != NULL)
					*longind = option_index;
				if (pfound->flag)
				{
					*(pfound->flag) = pfound->val;
					return 0;
				}
				return pfound->val;
			}
			nextchar = NULL;
			return 'W';	/* Let the application handle it.   */
		}
		if (temp[1] == ':')
		{
			if (temp[2] == ':')
			{
				/* This is an option that accepts an argument optionally.  */
				if (*nextchar != '\0')
				{
					optarg = nextchar;
					optind++;
				}
				else
					optarg = NULL;
				nextchar = NULL;
			}
			else
			{
				/* This is an option that requires an argument.  */
				if (*nextchar != '\0')
				{
					optarg = nextchar;
					/* If we end this ARGV-element by taking the rest as an arg,
					   we must advance to the next element now.  */
					optind++;
				}
				else if (optind == argc)
				{
					if (opterr)
					{
						/* 1003.2 specifies the format of this message.  */
						fprintf(stderr,
							"%s: option requires an argument -- %c\n",
							argv[0], c);
					}
					optopt = c;
					if (optstring[0] == ':')
						c = ':';
					else
						c = '?';
				}
				else
					/* We already incremented `optind' once;
					   increment it again when taking next ARGV-elt as argument.  */
					optarg = argv[optind++];
				nextchar = NULL;
			}
		}
		return c;
	}
}

int
    getopt(int argc, char *const *argv, const char *optstring)
{
	return _getopt_internal(argc, argv, optstring,
				(const struct option *) 0,
				(int *) 0,
				0);
}

/* getopt_long and getopt_long_only entry points for GNU getopt.
   Copyright (C) 1987,88,89,90,91,92,93,94,96,97 Free Software Foundation, Inc.

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

int
    getopt_long(int argc, char *const *argv,  const char *options,
     const struct option *long_options, int *opt_index)
{
	return _getopt_internal(argc, argv, options, long_options, opt_index, 0);
}

/* Like getopt_long, but '-' as well as '--' can indicate a long option.
   If an option that starts with '-' (not '--') doesn't match a long option,
   but does match a short option, it is parsed as a short option
   instead.  */

int
    getopt_long_only(int argc, char *const *argv, const char *options,
     const struct option *long_options, int *opt_index)
{
	return _getopt_internal(argc, argv, options, long_options, opt_index, 1);
}








