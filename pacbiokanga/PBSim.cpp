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

// The code is this module is copied, with relatively minor modifications, from the source code (PBSIM 1.0.3) released in support
// of the following paper and full attribution is hereby given to the authors:
// "PBSIM: PacBio reads simulator - toward accurate genome assembly" 
// Yukiteru Ono, Kiyoshi Asai and Michiaki Hamada 
// Bioinformatics (2013) 29 (1) 119-121.  
// doi: 10.1093/bioinformatics/bts649 

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

#include "pacbiokanga.h"

#include <math.h>
#define _USE_MATH_DEFINES // for C++
#include <cmath>
#include "PBSim.h"


CPBSim PBSim;

void PBSim_print_help(void);

/////////////////////////////////////////
// Main                                //
/////////////////////////////////////////

#ifdef _WIN32
int ProcPBSim(int argc, char* argv[])
{
// determine my process name
_splitpath(argv[0],NULL,NULL,gszProcName,NULL);
#else
int
ProcPBSim(int argc, char** argv)
{
// determine my process name
CUtility::splitpath((char *)argv[0],NULL,gszProcName);
#endif

	// initialise diagnostics log system to a default log file
	// if can't open then will ignore the error and continue with log messages to screen only
	if(!gDiagnostics.Open("pacbiokanga_PBSim.log",eDLInfo,eDLInfo,true))
		gDiagnostics.Open(NULL,eDLInfo,eDLNone,true);
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Subprocess %s Version %s starting",gpszSubProcess->pszName,cpszProgVer);

	// give credit where credit is due!
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"This module integrates code derived from PBSIM 1.0.3 and full attribution is hereby given to the authors:\n \"PBSIM: PacBio reads simulator - toward accurate genome assembly\",Yukiteru Ono, Kiyoshi Asai and Michiaki Hamada, Bioinformatics (2013) 29 (1) 119-121.");

  // show user current resource limits
#ifndef _WIN32
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Resources: %s",CUtility::ReportResourceLimits());
#endif


  struct sim_t sim;
  char *tp, *tmp_buf;
  long num;
  long ratio;


  memset(&sim, 0, sizeof(sim));

  sim.seed = (unsigned int)time(NULL);
  sim.len_min = FASTQ_LEN_MIN;

  // Variables for Option
  int opt, option_index;
  struct option long_options[] = {
    {"sample-fastq", 1, NULL, 0},
    {"data-type", 1, NULL, 0},
    {"depth", 1, NULL, 0},
    {"length-mean", 1, NULL, 0},
    {"length-sd", 1, NULL, 0},
    {"length-min", 1, NULL, 0},
    {"length-max", 1, NULL, 0},
    {"accuracy-mean", 1, NULL, 0},
    {"accuracy-sd", 1, NULL, 0},
    {"accuracy-min", 1, NULL, 0},
    {"accuracy-max", 1, NULL, 0},
    {"difference-ratio", 1, NULL, 0},
    {"model_qc", 1, NULL, 0},
    {"prefix", 1, NULL, 0},
    {"sample-profile-id", 1, NULL, 0},
    {"seed", 1, NULL, 0},
    {0, 0, 0, 0}
  };

  // Option parsing
  option_index = 0;
  while ((opt = getopt_long(argc, argv, "", long_options, &option_index)) != -1) {
    switch (opt) {
    case 0:
      sim.set_flg[option_index] = 1;

      switch (option_index) {
      case 0:
		strcpy(sim.szFastqFileName,optarg);
        break;

      case 1:
        if (strncmp(optarg, "CLR", 3) == 0) {
          sim.data_type = DATA_TYPE_CLR;
        } else if (strncmp(optarg, "CCS", 3) == 0) {
          sim.data_type = DATA_TYPE_CCS;
        } else {
          gDiagnostics.DiagOut(eDLFatal, gszProcName,  "ERROR (data-type: %s): Acceptable value is CLR or CCS.\n", optarg);
          exit(-1);
        }
        break;

      case 2:
        sim.depth = atof(optarg);
        if (sim.depth <= 0.0) {
          gDiagnostics.DiagOut(eDLFatal, gszProcName,  "ERROR (depth: %s): Acceptable range is more than 0.\n", optarg);
          exit(-1);
        }
        break;

      case 3:
        sim.len_mean = atof(optarg);
        if ((sim.len_mean < FASTQ_LEN_MIN) || (sim.len_mean > FASTQ_LEN_MAX)) {
          gDiagnostics.DiagOut(eDLFatal, gszProcName,  "ERROR (length-mean: %s): Acceptable range is %d-%d.\n", optarg,FASTQ_LEN_MIN, FASTQ_LEN_MAX);
          exit(-1);
        }
        break;

      case 4:
        sim.len_sd = atof(optarg);
        if ((sim.len_sd < 0) || (sim.len_sd > FASTQ_LEN_MAX)) {
          gDiagnostics.DiagOut(eDLFatal, gszProcName, "ERROR (length-sd: %s): Acceptable range is 0-%d.\n",
            optarg, FASTQ_LEN_MAX);
          exit(-1);
        }
        break;

      case 5:
        if (strlen(optarg) >= 8) {
          gDiagnostics.DiagOut(eDLFatal, gszProcName,  "ERROR (length-min: %s): Acceptable range is %d-%d.\n", optarg, FASTQ_LEN_MIN,FASTQ_LEN_MAX);
          exit(-1);
        }
        sim.len_min = atoi(optarg);
        if ((sim.len_min < FASTQ_LEN_MIN) || (sim.len_min > FASTQ_LEN_MAX)) {
          gDiagnostics.DiagOut(eDLFatal, gszProcName,  "ERROR (length-min: %s): Acceptable range is %d-%d.\n", optarg, FASTQ_LEN_MIN,FASTQ_LEN_MAX);
          exit(-1);
        }
        break;

      case 6:
        if (strlen(optarg) >= 8) {
          gDiagnostics.DiagOut(eDLFatal, gszProcName,  "ERROR (length-max: %s): Acceptable range is %d-%d.\n", optarg,FASTQ_LEN_MIN, FASTQ_LEN_MAX);
          exit(-1);
        }
        sim.len_max = atoi(optarg);
        if ((sim.len_max < 1) || (sim.len_max > FASTQ_LEN_MAX)) {
          gDiagnostics.DiagOut(eDLFatal, gszProcName,  "ERROR (length-max: %s): Acceptable range is %d-%d.\n", optarg, FASTQ_LEN_MIN,FASTQ_LEN_MAX);
          exit(-1);
        }
        break;

      case 7:
        sim.accuracy_mean = atof(optarg);
        if ((sim.accuracy_mean < 0.0) || (sim.accuracy_mean > 1.0)) {
          gDiagnostics.DiagOut(eDLFatal, gszProcName, "ERROR (accuracy-mean: %s): Acceptable range is 0.0-1.0.\n", optarg);
          exit(-1);
        }
        break;

      case 8:
        sim.accuracy_sd = atof(optarg);
        if ((sim.accuracy_sd < 0.0) || (sim.accuracy_sd > 1.0)) {
          gDiagnostics.DiagOut(eDLFatal, gszProcName,  "ERROR (accuracy-sd: %s): Acceptable range is 0.0-1.0.\n", optarg);
          exit(-1);
        }
        break;

      case 9:
        sim.accuracy_min = atof(optarg);
        if ((sim.accuracy_min < 0.0) || (sim.accuracy_min > 1.0)) {
          gDiagnostics.DiagOut(eDLFatal, gszProcName, "ERROR (accuracy-min: %s): Acceptable range is 0.0-1.0.\n", optarg);
          exit(-1);
        }
        break;

      case 10:
        sim.accuracy_max = atof(optarg);
        if ((sim.accuracy_max < 0.0) || (sim.accuracy_max > 1.0)) {
          gDiagnostics.DiagOut(eDLFatal, gszProcName,  "ERROR (accuracy-max: %s): Acceptable range is 0.0-1.0.\n", optarg);
          exit(-1);
        }
        break;

      case 11:
        if ((tmp_buf = (char *)malloc(strlen(optarg) + 1)) == 0) {
          gDiagnostics.DiagOut(eDLFatal, gszProcName,  "ERROR: Cannot allocate memory.\n");
          exit(-1);
        }
        strcpy(tmp_buf, optarg);
        num = 0;
        tp = strtok(tmp_buf, ":");
        while (num < 3) {
          if (tp == NULL) {
            gDiagnostics.DiagOut(eDLFatal, gszProcName,  "ERROR (difference-ratio: %s): Format is sub:ins:del.\n", optarg);
            exit(-1);
          }
          if (strlen(tp) >= 5) {
            gDiagnostics.DiagOut(eDLFatal, gszProcName, "ERROR (difference-ratio: %s): Acceptable range is 0-%d.\n", optarg, RATIO_MAX);
            exit(-1);
          }
          ratio = atoi(tp);
          if ((ratio < 0) || (ratio > RATIO_MAX)) {
            gDiagnostics.DiagOut(eDLFatal, gszProcName, "ERROR (difference-ratio: %s): Acceptable range is 0-%d.\n", optarg, RATIO_MAX);
            exit(-1);
          }
          if (num == 0) {
            sim.sub_ratio = ratio;
          } else if (num == 1) {
            sim.ins_ratio = ratio;
          } else if (num == 2) {
            sim.del_ratio = ratio;
          }
          num ++;
          tp = strtok(NULL, ":");
        }
        break;

      case 12:
        if ((sim.model_qc_file = (char *)malloc(strlen(optarg) + 1)) == 0) {
          gDiagnostics.DiagOut(eDLFatal, gszProcName,  "ERROR: Cannot allocate memory.\n");
          exit(-1);
        }
        strcpy(sim.model_qc_file, optarg);
        break;

      case 13:
        if ((sim.prefix = (char *)malloc(strlen(optarg) + 1)) == 0) {
          gDiagnostics.DiagOut(eDLFatal, gszProcName,  "ERROR: Cannot allocate memory.\n");
          exit(-1);
        }
        strcpy(sim.prefix, optarg);
        break;

      case 14:
        if ((sim.profile_id = (char *)malloc(strlen(optarg) + 1)) == 0) {
          gDiagnostics.DiagOut(eDLFatal, gszProcName, "ERROR: Cannot allocate memory.\n");
          exit(-1);
        }
        strcpy(sim.profile_id, optarg);
        break;

      case 15:
        sim.seed = (unsigned int)atoi(optarg);
        break;

      default:
        break;
      }

    default:
      break;
    }
  }


  if (argv[optind] == '\0') {
    PBSim_print_help();
    exit(-1);
  }

  strcpy(sim.szref_file, argv[optind]);

return(PBSim.Process(&sim));
}

CPBSim::CPBSim()
{
fp_filtered=NULL;
fp_stats = NULL;
fp_ref = NULL;
fp_fq = NULL;
fp_maf = NULL;
memset(&sim,0,sizeof(sim));

memset(&fastq,0,sizeof(fastq));
memset(&ref,0,sizeof(ref));
memset(&mut,0,sizeof(mut));

memset(model_qc,0,sizeof(model_qc));
memset(freq_len,0,sizeof(freq_len));
memset(freq_accuracy,0,sizeof(freq_accuracy));
}

CPBSim::~CPBSim()
{


}


TRandomCombined<CRandomMother,CRandomMersenne> RGseeds((int)time(0));

int CPBSim::Process(struct sim_t *psim)
{
  long i;
  long t1, t2;
  sim = *psim;

  if ((fastq.file = (char *)malloc(strlen(sim.szFastqFileName) + 1)) == 0) {
          gDiagnostics.DiagOut(eDLFatal, gszProcName, "ERROR: Cannot allocate memory.\n");
          exit(-1);
        }
  strcpy(fastq.file, sim.szFastqFileName);
  t1 = get_time();
 RGseeds.RandomInit((int)sim.seed);

  // Quality code to error probability
  for (i=0; i<=93; i++) {
    qc[i].prob = pow(10, (double)i / -10);
    qc[i].character = (char)(i+33);
  }

  // Setting of simulation parameters     
  if (set_sim_param() == PBFAILED) {
    exit(-1);
  }
  print_sim_param();

  // FASTQ
  if (sim.process == PROCESS_SAMPLING) {
    if ((fp_filtered = tmpfile()) == NULL) {
      gDiagnostics.DiagOut(eDLFatal, gszProcName, "ERROR: Cannot open temporary file\n");
      return PBFAILED;
    }

    if (get_fastq_inf() == PBFAILED) {
      exit(-1);
    }
    print_fastq_stats();
  } else if (sim.process == PROCESS_SAMPLING_STORE) {
    if ((fp_filtered = fopen(sim.profile_fq, "w+")) == NULL) {
      gDiagnostics.DiagOut(eDLFatal, gszProcName, "ERROR: Cannot open sample_profile\n");
      return PBFAILED;
    }
    if ((fp_stats = fopen(sim.profile_stats, "w+")) == NULL) {
      gDiagnostics.DiagOut(eDLFatal, gszProcName,  "ERROR: Cannot open sample_profile\n");
      return PBFAILED;
    }

    if (get_fastq_inf() == PBFAILED) {
      exit(-1);
    }
    print_fastq_stats();
  } else if (sim.process == PROCESS_SAMPLING_REUSE) {
    if ((fp_filtered = fopen(sim.profile_fq, "r")) == NULL) {
      gDiagnostics.DiagOut(eDLFatal, gszProcName, "ERROR: Cannot open sample_profile\n");
      return PBFAILED;
    }
    if ((fp_stats = fopen(sim.profile_stats, "r")) == NULL) {
      gDiagnostics.DiagOut(eDLFatal, gszProcName, "ERROR: Cannot open sample_profile\n");
      return PBFAILED;
    }

    if (get_fastq_inf() == PBFAILED) {
      exit(-1);
    }
    print_fastq_stats();
  }

  // Quality code model
  if (sim.process == PROCESS_MODEL) {
    if (set_model_qc() == PBFAILED) {
      exit(-1);
    }
  }

  // Reference sequence
  if ((ref.file = (char *)malloc(strlen(sim.szref_file) + 1)) == 0) {
    gDiagnostics.DiagOut(eDLFatal, gszProcName, "ERROR: Cannot allocate memory.\n");
    exit(-1);
  }
  strcpy(ref.file, sim.szref_file);

  if (get_ref_inf() == PBFAILED) {
    exit(-1);
  }

  // Set mutation parameters and varianeces
  if (set_mut() == PBFAILED) {
    exit(-1);
  }

  // Creating simulated reads
  for (ref.num=1; ref.num<=ref.num_seq; ref.num++) {
    if (get_ref_seq() == PBFAILED) {
      exit(-1);
    }

    init_sim_res();

    sprintf(sim.outfile_fq, "%s_%04ld.fastq", sim.prefix, ref.num);
    if ((fp_fq = fopen(sim.outfile_fq, "w")) == NULL) {
      gDiagnostics.DiagOut(eDLFatal, gszProcName, "ERROR: Cannot open output file: %s\n", sim.outfile_fq);
      return PBFAILED;
    }

    sprintf(sim.outfile_maf, "%s_%04ld.maf.txt", sim.prefix, ref.num);
    if ((fp_maf = fopen(sim.outfile_maf, "w")) == NULL) {
      gDiagnostics.DiagOut(eDLFatal, gszProcName, "ERROR: Cannot open output file: %s\n", sim.outfile_maf);
      return PBFAILED;
    }

    sim.len_quota = (long long)(sim.depth * ref.len);

    if (sim.process == PROCESS_MODEL) {
      if (simulate_by_model() == PBFAILED) {
        exit(-1);
      }
    } else {
      if (simulate_by_sampling() == PBFAILED) {
        exit(-1);
      }
    }

    print_simulation_stats();

    fclose(fp_fq);
    fclose(fp_maf);
  }

  if ((sim.process == PROCESS_SAMPLING_STORE) || (sim.process == PROCESS_SAMPLING_REUSE)) {
    fclose(fp_filtered);
    fclose(fp_stats);
  }

  t2 = get_time();

 gDiagnostics.DiagOut(eDLInfo, gszProcName,  ":::: System utilization ::::");
  gDiagnostics.DiagOut(eDLInfo, gszProcName,  "Elapsed time(s) : %ld", t2 - t1);

  return(0);
}

///////////////////////////////////////
// Function: trim - Remove "\n"      //
///////////////////////////////////////

int CPBSim::trim(char *line) {
  int end_pos = (int)strlen(line) - 1;

  if (line[end_pos] == '\n') {
    line[end_pos] = '\0';
    return 1;
  }
  return 0;
}

///////////////////////////////////////////////////////
// Function: get_ref_inf - Get reference information //
///////////////////////////////////////////////////////

int CPBSim::get_ref_inf(void) {
  FILE *fp;
  char *pLine;
  int ret;
  long max_len = 0;

  pLine = new char [BUF_SIZE+1];
  gDiagnostics.DiagOut(eDLInfo, gszProcName, ":::: Reference stats ::::");
	gDiagnostics.DiagOut(eDLInfo, gszProcName,"file name : %s\n", ref.file);


  if ((fp = fopen(ref.file, "r")) == NULL) {
    gDiagnostics.DiagOut(eDLFatal, gszProcName,"ERROR: Cannot open file: %s", ref.file);
	delete pLine;
    return PBFAILED;
  }

  ref.num_seq = 0;
  ref.len = 0;

  while (fgets(pLine, BUF_SIZE, fp) != NULL) {
    ret = trim(pLine);

    if (pLine[0] == '>') {
      if (ref.num_seq != 0) {
        if (ref.len < REF_SEQ_LEN_MIN) {
          gDiagnostics.DiagOut(eDLFatal, gszProcName, "ERROR: Reference is too short. Acceptable length >= %ld", REF_SEQ_LEN_MIN);
		  delete pLine;
          return PBFAILED;
        }
        gDiagnostics.DiagOut(eDLInfo, gszProcName, "ref.%d (len:%d) : %s", ref.num_seq, ref.len, ref.id);
        fclose(fp_ref);
        if (ref.len > max_len) {
          max_len = ref.len;
        }
      }

      ref.num_seq ++;
      if (ref.num_seq > REF_SEQ_NUM_MAX) {
        gDiagnostics.DiagOut(eDLFatal, gszProcName,"ERROR: References are too many. Max number of reference is %ld", REF_SEQ_NUM_MAX);
        delete pLine;
        return PBFAILED;
      }

      strncpy(ref.id, pLine + 1, REF_ID_LEN_MAX);
      ref.id[REF_ID_LEN_MAX] = '\0';

      sprintf(sim.outfile_ref, "%s_%04ld.ref", sim.prefix, ref.num_seq);
      if ((fp_ref = fopen(sim.outfile_ref, "w")) == NULL) {
        gDiagnostics.DiagOut(eDLFatal, gszProcName, "ERROR: Cannot open output file: %s", sim.outfile_ref);
		delete pLine;
        return PBFAILED;
      }

      ref.len = 0;

      while (ret != EXISTS_LINE_FEED) {
        if (fgets(pLine, BUF_SIZE, fp) == NULL) {
          break;
        }
        ret = trim(pLine);
      }

      fprintf(fp_ref, ">%s\n", ref.id);
    } else {
      ref.len += (long)strlen(pLine);

      if (ref.len > REF_SEQ_LEN_MAX) {
        gDiagnostics.DiagOut(eDLFatal, gszProcName,"ERROR: Reference is too long. Acceptable length <= %ld", REF_SEQ_LEN_MAX);
		delete pLine;
        return PBFAILED;
      }

      fprintf(fp_ref, "%s\n", pLine);
    }
  }
  delete pLine;
  fclose(fp);
  
  if (ref.len < REF_SEQ_LEN_MIN) {
    gDiagnostics.DiagOut(eDLFatal, gszProcName, "ERROR: Reference is too short. Acceptable length >= %ld", REF_SEQ_LEN_MIN);
    return PBFAILED;
  }
  gDiagnostics.DiagOut(eDLFatal, gszProcName, "ref.%d (len:%d) : %s", ref.num_seq, ref.len, ref.id);
  fclose(fp_ref);
  if (ref.len > max_len) {
    max_len = ref.len;
  }

  if ((ref.seq = (char *)malloc(max_len + 1)) == 0) {
    gDiagnostics.DiagOut(eDLFatal, gszProcName, "ERROR: Cannot allocate memory");
    return PBFAILED;
  }

  return PBSUCCEEDED;
}

////////////////////////////////////////////////////
// Function: get_ref_seq - Get reference sequence //
////////////////////////////////////////////////////

int CPBSim::get_ref_seq(void) {
  FILE *fp;
  char *pLine;
  long offset = 0;
  long copy_size;
  int ret;
  

  sprintf(sim.outfile_ref, "%s_%04ld.ref", sim.prefix, ref.num);

  if ((fp = fopen(sim.outfile_ref, "r")) == NULL) {
    gDiagnostics.DiagOut(eDLFatal, gszProcName,"ERROR: Cannot open file: %s", sim.outfile_ref);
    return PBFAILED;
  }

  pLine = new char [BUF_SIZE];

  while (fgets(pLine, BUF_SIZE, fp) != NULL) {
    ret = trim(pLine);

    if (pLine[0] == '>') {
      while (ret != EXISTS_LINE_FEED) {
        if (fgets(pLine, BUF_SIZE, fp) == NULL) {
          break;
        }
        ret = trim(pLine);
      }
    } else {
      copy_size = (long)strlen(pLine);
      memcpy(ref.seq + offset, pLine, copy_size);
      offset += copy_size;
    }
  }
  fclose(fp);
  delete pLine;

  ref.seq[offset] = '\0';
  ref.len = (long)strlen(ref.seq);

  return PBSUCCEEDED;
}

/////////////////////////////////////////////////////
// Function: get_fastq_inf - Get FASTQ information //
/////////////////////////////////////////////////////

int CPBSim::get_fastq_inf(void) {
  FILE *fp;
  char *tp, *item;
  char *pLine;
  char *pqc_tmp;
  long len;
  double prob;
  double accuracy;
  double accuracy_total = 0;
  long value;
  double variance;
  long i;
  int line_num;

  pLine = new char [BUF_SIZE];
  memset(pLine,0,BUF_SIZE);
  pqc_tmp = new char [FASTQ_LEN_MAX];
  memset(pqc_tmp,0,FASTQ_LEN_MAX);
  memset(freq_len,0,sizeof(freq_len));
  memset(freq_accuracy,0,sizeof(freq_accuracy));
 
  fastq.num = 0;
  fastq.len_min = LONG_MAX;
  fastq.len_max = 0;
  fastq.len_total = 0;
  fastq.num_filtered = 0;
  fastq.len_min_filtered = LONG_MAX;
  fastq.len_max_filtered = 0;
  fastq.len_total_filtered = 0;

  if (sim.process == PROCESS_SAMPLING_REUSE) {
    while (fgets(pLine, BUF_SIZE, fp_stats) != NULL) {
      trim(pLine);
      tp = strtok(pLine, "\t");
      item = tp;
      tp = strtok(NULL, "\t");

      if (strcmp(item, "num") == 0) {
        fastq.num_filtered = atol(tp);
      } else if (strcmp(item, "len_total") == 0) {
        fastq.len_total_filtered = atol(tp);
      } else if (strcmp(item, "len_min") == 0) {
        fastq.len_min_filtered = atol(tp);
      } else if (strcmp(item, "len_max") == 0) {
        fastq.len_max_filtered = atol(tp);
      } else if (strcmp(item, "len_mean") == 0) {
        fastq.len_mean_filtered = atof(tp);
      } else if (strcmp(item, "len_sd") == 0) {
        fastq.len_sd_filtered = atof(tp);
      } else if (strcmp(item, "accuracy_mean") == 0) {
        fastq.accuracy_mean_filtered = atof(tp);
      } else if (strcmp(item, "accuracy_sd") == 0) {
        fastq.accuracy_sd_filtered = atof(tp);
      }
    }
  } else {
    if ((fp = fopen(fastq.file, "r")) == NULL) {
      gDiagnostics.DiagOut(eDLFatal, gszProcName,"ERROR: Cannot open file: %s", fastq.file);
	   delete pLine;
		delete pqc_tmp;
      return PBFAILED;
    }

    pqc_tmp[0] = '\0';
    len = 0;
    line_num = 0;

    while (fgets(pLine, BUF_SIZE, fp) != NULL) {
      if (trim(pLine) == EXISTS_LINE_FEED) {
        line_num ++;

        if (line_num == 4) {
          len += (long)strlen(pLine);

          if (len > FASTQ_LEN_MAX) {
            gDiagnostics.DiagOut(eDLFatal, gszProcName, "ERROR: fastq is too long. Max acceptable length is %ld", FASTQ_LEN_MAX);
             delete pLine;
			delete pqc_tmp;
            return PBFAILED;
          }

          fastq.num ++;
          fastq.len_total += len;

          if (fastq.num > FASTQ_NUM_MAX) {
            gDiagnostics.DiagOut(eDLFatal, gszProcName,"ERROR: fastq is too many. Max acceptable number is %ld", FASTQ_NUM_MAX);
            delete pLine;
			delete pqc_tmp;
            return PBFAILED;
          }

          if (len > fastq.len_max) {
            fastq.len_max = len;
          }
          if (len < fastq.len_min) {
            fastq.len_min = len;
          }

          if ((len >= sim.len_min) && (len <= sim.len_max)) {
            strcat(pqc_tmp, pLine);
            prob = 0.0;
            for (i=0; i<len; i++) {
              prob += qc[(int)pqc_tmp[i] - 33].prob;
            }
            accuracy = 1.0 - (prob / len);
			if(accuracy < 0.0)
				accuracy = 0.0;
		    else
				if(accuracy > 1.0)
					accuracy = 1.0;
            if ((accuracy >= sim.accuracy_min) && (accuracy <= sim.accuracy_max)) {
              accuracy_total += accuracy;
              fastq.num_filtered ++;
              fastq.len_total_filtered += len;

              freq_len[len] ++;
              value = (int)(accuracy * 100000); 
              freq_accuracy[value] ++;

              fprintf(fp_filtered, "%s\n", pqc_tmp);

              if (len > fastq.len_max_filtered) {
                fastq.len_max_filtered = len;
              }
              if (len < fastq.len_min_filtered) {
                fastq.len_min_filtered = len;
              }
            }
          }

          line_num = 0;
          pqc_tmp[0] = '\0';
          len = 0;
        }
      } else {
        if (line_num == 3) {
          len += (long)strlen(pLine);
          if (len > FASTQ_LEN_MAX) {
            gDiagnostics.DiagOut(eDLFatal, gszProcName,"ERROR: fastq is too long. Max acceptable length is %ld", FASTQ_LEN_MAX);
			 delete pLine;
			 delete pqc_tmp;
            return PBFAILED;
          }
          strcat(pqc_tmp, pLine);
        }
      }
    }

    fclose(fp);

    if (fastq.num_filtered < 1) {
      gDiagnostics.DiagOut(eDLFatal, gszProcName,"ERROR: there is no sample-fastq in the valid range of length and accuracy");
	 delete pLine;
	 delete pqc_tmp;
      return PBFAILED;
    }

    fastq.len_mean_filtered = (double)fastq.len_total_filtered / fastq.num_filtered;
    fastq.accuracy_mean_filtered = accuracy_total / fastq.num_filtered;

    variance = 0.0;
    for (i=0; i<=sim.len_max; i++) {
      if (freq_len[i] > 0) { 
        variance += pow((fastq.len_mean_filtered - i), 2) * freq_len[i];
      }
    }
    fastq.len_sd_filtered = sqrt(variance / fastq.num_filtered);

    variance = 0.0;
    for (i=0; i<=100000; i++) {
      if (freq_accuracy[i] > 0) { 
        variance += pow((fastq.accuracy_mean_filtered - i * 0.00001), 2) * freq_accuracy[i];
      }
    }
    fastq.accuracy_sd_filtered = sqrt(variance / fastq.num_filtered);

    if (sim.process == PROCESS_SAMPLING_STORE) {
      fprintf(fp_stats, "num\t%ld\n", fastq.num_filtered);
      fprintf(fp_stats, "len_total\t%lld\n", fastq.len_total_filtered);
      fprintf(fp_stats, "len_min\t%ld\n", fastq.len_min_filtered);
      fprintf(fp_stats, "len_max\t%ld\n", fastq.len_max_filtered);
      fprintf(fp_stats, "len_mean\t%f\n", fastq.len_mean_filtered);
      fprintf(fp_stats, "len_sd\t%f\n", fastq.len_sd_filtered);
      fprintf(fp_stats, "accuracy_mean\t%f\n", fastq.accuracy_mean_filtered);
      fprintf(fp_stats, "accuracy_sd\t%f\n", fastq.accuracy_sd_filtered);
    }
  }

 delete pLine;
 delete pqc_tmp;
 return PBSUCCEEDED;
}

/////////////////////////////////////////////////////
// Function: print_fastq_stats - Print FASTQ stats //
/////////////////////////////////////////////////////

void CPBSim::print_fastq_stats(void) {
  gDiagnostics.DiagOut(eDLInfo, gszProcName, ":::: FASTQ stats ::::");

  if (sim.process == PROCESS_SAMPLING_REUSE) {
    gDiagnostics.DiagOut(eDLInfo, gszProcName, "file name : %s", sim.profile_fq);
  } else {
    gDiagnostics.DiagOut(eDLInfo, gszProcName, "file name : %s", fastq.file);
    gDiagnostics.DiagOut(eDLInfo, gszProcName, ":: all reads ::");
    gDiagnostics.DiagOut(eDLInfo, gszProcName, "read num. : %ld", fastq.num);
    gDiagnostics.DiagOut(eDLInfo, gszProcName, "read total length : %lld", fastq.len_total);
    gDiagnostics.DiagOut(eDLInfo, gszProcName, "read min length : %ld", fastq.len_min);
    gDiagnostics.DiagOut(eDLInfo, gszProcName, "read max length : %ld", fastq.len_max);
  }

  gDiagnostics.DiagOut(eDLInfo, gszProcName, ":: filtered reads ::");
  gDiagnostics.DiagOut(eDLInfo, gszProcName, "read num. : %ld", fastq.num_filtered);
  gDiagnostics.DiagOut(eDLInfo, gszProcName, "read total length : %lld", fastq.len_total_filtered);
  gDiagnostics.DiagOut(eDLInfo, gszProcName, "read min length : %ld", fastq.len_min_filtered);
  gDiagnostics.DiagOut(eDLInfo, gszProcName, "read max length : %ld", fastq.len_max_filtered);
  gDiagnostics.DiagOut(eDLInfo, gszProcName, "read length mean (SD) : %f (%f)",
    fastq.len_mean_filtered, fastq.len_sd_filtered);
  gDiagnostics.DiagOut(eDLInfo, gszProcName, "read accuracy mean (SD) : %f (%f)",
    fastq.accuracy_mean_filtered, fastq.accuracy_sd_filtered);
}

//////////////////////////////////////////////////////////
// Function: init_sim_res - Initiate simulation results //
//////////////////////////////////////////////////////////

void CPBSim::init_sim_res(void) {
  sim.res_num = 0;
  sim.res_len_total = 0;
  sim.res_sub_num = 0;
  sim.res_ins_num = 0;
  sim.res_del_num = 0;
  sim.res_len_min = LONG_MAX;
  sim.res_len_max = 0;
}

/////////////////////////////////////////////////////////
// Function: set_sim_param - Set simulation parameters //
/////////////////////////////////////////////////////////

int CPBSim::set_sim_param(void) {
  FILE *fp;
  long sum;

  // data-type
  if (!(sim.set_flg[1])) {
    sim.data_type = DATA_TYPE_CLR;
  }

  // depth
  if (!(sim.set_flg[2])) {
    sim.depth = (sim.data_type == DATA_TYPE_CLR) ? 20.0 : 50.0;
  }

  // length-mean
  if (!(sim.set_flg[3])) {
    sim.len_mean = (sim.data_type == DATA_TYPE_CLR) ? 3000 : 450;
  }

  // length-sd
  if (!(sim.set_flg[4])) {
    sim.len_sd = (sim.data_type == DATA_TYPE_CLR) ? 2300 : 170;
  }

  // length-min
  if (!(sim.set_flg[5])) {
    sim.len_min = (sim.data_type == DATA_TYPE_CLR) ? 100 : 100;
  }

  // length-max
  if (!(sim.set_flg[6])) {
    sim.len_max = (sim.data_type == DATA_TYPE_CLR) ? 25000 : 2500;
  }

  // accuracy-mean
  if (sim.data_type == DATA_TYPE_CLR) {
    if (sim.set_flg[7]) {
      sim.accuracy_mean = int(sim.accuracy_mean * 100) * 0.01;
    } else {
      sim.accuracy_mean = 0.78;
    }
  } else {
    sim.accuracy_mean = 0.98;
  }

  // accuracy-sd
  if (sim.data_type == DATA_TYPE_CLR) {
    if (sim.set_flg[8]) {
      sim.accuracy_sd = int(sim.accuracy_sd * 100) * 0.01;
    } else {
      sim.accuracy_sd = 0.02;
    }
  } else {
    sim.accuracy_sd = 0.02;
  }

  // accuracy-min
  if (sim.data_type == DATA_TYPE_CLR) {
    if (sim.set_flg[9]) {
      sim.accuracy_min = int(sim.accuracy_min * 100) * 0.01;
    } else {
      sim.accuracy_min = 0.75;
    }
  } else {
    sim.accuracy_min = 0.75;
  }

  // accuracy-max
  if (sim.data_type == DATA_TYPE_CLR) {
    if (sim.set_flg[10]) {
      sim.accuracy_max = int(sim.accuracy_max * 100) * 0.01;
    } else {
      sim.accuracy_max = 1.0;
    }
  } else {
    sim.accuracy_max = 1.0;
  }

  // difference-ratio
  if (!(sim.set_flg[11])) {
    if (sim.data_type == DATA_TYPE_CLR) {
      sim.sub_ratio = 10;
      sim.ins_ratio = 60;
      sim.del_ratio = 30;
    } else {
      sim.sub_ratio = 6;
      sim.ins_ratio = 21;
      sim.del_ratio = 73;
    }
  }

  sum = sim.sub_ratio + sim.ins_ratio + sim.del_ratio;
  sim.sub_rate = (double)sim.sub_ratio / sum;
  sim.ins_rate = (double)sim.ins_ratio / sum;
  sim.del_rate = (double)sim.del_ratio / sum;

  // prefix and outfile
  if (!(sim.set_flg[13])) {
    if ((sim.prefix = (char *)malloc(10)) == 0) {
      gDiagnostics.DiagOut(eDLFatal, gszProcName,  "ERROR: Cannot allocate memory");
      exit(-1);
    }
    strcpy(sim.prefix, "sd");
  }

  if ((sim.outfile_ref = (char *)malloc(strlen(sim.prefix) + 50)) == 0) {
    gDiagnostics.DiagOut(eDLFatal, gszProcName, "ERROR: Cannot allocate memory");
    return PBFAILED;
  }

  if ((sim.outfile_fq = (char *)malloc(strlen(sim.prefix) + 50)) == 0) {
    gDiagnostics.DiagOut(eDLFatal, gszProcName, "ERROR: Cannot allocate memory");
    return PBFAILED;
  }

  if ((sim.outfile_maf = (char *)malloc(strlen(sim.prefix) + 50)) == 0) {
    gDiagnostics.DiagOut(eDLFatal, gszProcName, "ERROR: Cannot allocate memory");
    return PBFAILED;
  }

  // profile
  if (sim.set_flg[14]) {
    if ((sim.profile_fq = (char *)malloc(strlen(sim.profile_id) + 50)) == 0) {
      gDiagnostics.DiagOut(eDLFatal, gszProcName,  "ERROR: Cannot allocate memory.\n");
      return PBFAILED;
    }

    if ((sim.profile_stats = (char *)malloc(strlen(sim.profile_id) + 50)) == 0) {
      gDiagnostics.DiagOut(eDLFatal, gszProcName, "ERROR: Cannot allocate memory.\n");
      return PBFAILED;
    }

    sprintf(sim.profile_fq, "sample_profile_%s.fastq", sim.profile_id);
    sprintf(sim.profile_stats, "sample_profile_%s.stats", sim.profile_id);
  }

  // length and accuracy
  if (sim.len_min > sim.len_max) {
    gDiagnostics.DiagOut(eDLFatal, gszProcName, "ERROR: length min(%ld) is greater than max(%ld).\n", sim.len_min, sim.len_max);
    return PBFAILED;
  }
  if (sim.accuracy_min > sim.accuracy_max) {
    gDiagnostics.DiagOut(eDLFatal, gszProcName, "ERROR: accuracy min(%f) is greater than max(%f).\n", sim.accuracy_min, sim.accuracy_max);
    return PBFAILED;
  }

  // process
  if (sim.set_flg[12]) {
    if ((sim.set_flg[0]) || (sim.set_flg[14])) {
      gDiagnostics.DiagOut(eDLFatal, gszProcName, "ERROR: either --sample-fastq(and/or --sample-profile-id)(sampling-based) or --model_qc(model-based) should be set.\n");
      return PBFAILED;
    }
    sim.process = PROCESS_MODEL;
  } else {
    if (sim.set_flg[0]) {
      if (sim.set_flg[14]) {
        sim.process = PROCESS_SAMPLING_STORE;
      } else {
        sim.process = PROCESS_SAMPLING;
      }
    } else {
      if (sim.set_flg[14]) {
        sim.process = PROCESS_SAMPLING_REUSE;
      } else {
        gDiagnostics.DiagOut(eDLFatal, gszProcName,  "ERROR: either --sample-fastq(and/or --sample-profile-id)(sampling-based) or --model_qc(model-based) should be set.\n");
        return PBFAILED;
      }
    }
  }

  // sample profile
  if (sim.process == PROCESS_SAMPLING_STORE) {
    if ((fp = fopen(sim.profile_fq, "r")) != NULL) {
      gDiagnostics.DiagOut(eDLFatal, gszProcName,  "ERROR: %s exists.\n", sim.profile_fq);
      fclose(fp);
      return PBFAILED;
    }
    if ((fp = fopen(sim.profile_stats, "r")) != NULL) {
      gDiagnostics.DiagOut(eDLFatal, gszProcName,  "ERROR: %s exists.\n", sim.profile_stats);
      fclose(fp);
      return PBFAILED;
    }
  }

  if (sim.process == PROCESS_SAMPLING_REUSE) {
    if ((fp = fopen(sim.profile_fq, "r")) == NULL) {
      gDiagnostics.DiagOut(eDLFatal, gszProcName,  "ERROR: %s does not exist.\n", sim.profile_fq);
      return PBFAILED;
    }
    fclose(fp);
    if ((fp = fopen(sim.profile_stats, "r")) == NULL) {
      gDiagnostics.DiagOut(eDLFatal, gszProcName, "ERROR: %s does not exist.\n", sim.profile_stats);
      return PBFAILED;
    }
    fclose(fp);
  }

  return PBSUCCEEDED;
}

////////////////////////////////////////////////////////
// Function: simulate_by_sampling - Simulate by model //
////////////////////////////////////////////////////////

int CPBSim::simulate_by_sampling(void) {
  long len;
  long long len_total = 0;
  long sampling_num, sampling_interval, sampling_value, sampling_residue;
  long num;
  long i, j;
  long value;
  double accuracy, accuracy_total = 0.0;
  double prob, variance;
  char id[128];
  int digit_num1[4], digit_num2[4], digit_num[4];

  memset(freq_len,0,sizeof(freq_len));
  memset(freq_accuracy,0,sizeof(freq_accuracy));

  for (i=0; i<=93; i++) {
    mut.sub_thre[i] = int((qc[i].prob * sim.sub_rate) * 1000000);
    mut.ins_thre[i] = int((qc[i].prob * (sim.sub_rate + sim.ins_rate)) * 1000000);
  }
  mut.del_thre = int((1.0 - fastq.accuracy_mean_filtered) * sim.del_rate * 1000000);

  sampling_num = (long)(sim.len_quota / fastq.len_total_filtered);
  sampling_residue = (long)(sim.len_quota % fastq.len_total_filtered);
  if (sampling_residue == 0) {
    sampling_interval = 1;
  } else {
    sampling_interval = (long)((double)(fastq.len_total_filtered / sampling_residue) * 2);
    if (sampling_interval > (long)(fastq.num_filtered/2)) {
      sampling_interval = (long)(fastq.num_filtered/2);
    }
  }

  // Make simulation data
  while (len_total < sim.len_quota) {
    rewind(fp_filtered);
	
	if(fastq.num_filtered >= 2)
		sampling_value =  RGseeds.IRandom(0,fastq.num_filtered-1);
	else
		sampling_value = 0;
    while (fgets(mut.qc, fastq.len_max_filtered + 2, fp_filtered) != NULL) {
      if (len_total >= sim.len_quota) {
        break;
      }

      trim(mut.qc);

      if (sampling_value % sampling_interval == 0) {
        num = sampling_num + 1;
      } else {
        num = sampling_num;
      }
      sampling_value ++;

      for (i=0; i<num; i++) {
        if (len_total >= sim.len_quota) {
          break;
        }

        mut.tmp_len_max = (long)(sim.len_quota - len_total);
        if (mut.tmp_len_max < sim.len_min) {
          mut.tmp_len_max = sim.len_min;
        }

        if (mutate() == PBFAILED) {
          return PBFAILED;
        }

        sim.res_num ++;
        len = (long)strlen(mut.new_seq);
        sim.res_len_total += len;
        len_total += len;
        freq_len[len] ++;

        if (len > sim.res_len_max) {
          sim.res_len_max = len;
        }
        if (len < sim.res_len_min) {
          sim.res_len_min = len;
        }

        prob = 0.0;
        for (j=0; j<len; j++) {
          prob += qc[(int)mut.new_qc[j] - 33].prob;
        }
        accuracy = 1.0 - (prob / len);
        accuracy_total += accuracy;
        value = (int)(accuracy * 100000);
        freq_accuracy[value] ++;

        sprintf(id, "S%ld_%ld", ref.num, sim.res_num);
        fprintf(fp_fq, "@%s\n%s\n+%s\n%s\n", id, mut.new_seq, id, mut.new_qc);

        digit_num1[0] = 3;
        digit_num2[0] = 1 + count_digit(sim.res_num);
        digit_num[0] = (digit_num1[0] >= digit_num2[0]) ? digit_num1[0] : digit_num2[0];

        digit_num1[1] = count_digit((mut.seq_left - 1));
        digit_num2[1] = 1;
        digit_num[1] = (digit_num1[1] >= digit_num2[1]) ? digit_num1[1] : digit_num2[1];

        digit_num1[2] = count_digit((mut.seq_right - mut.seq_left + 1));
        digit_num2[2] = count_digit(len);
        digit_num[2] = (digit_num1[2] >= digit_num2[2]) ? digit_num1[2] : digit_num2[2];

        digit_num1[3] = count_digit(ref.len);
        digit_num2[3] = count_digit(len);
        digit_num[3] = (digit_num1[3] >= digit_num2[3]) ? digit_num1[3] : digit_num2[3];

        fprintf(fp_maf, "a\ns ref"); 
        while (digit_num1[0] ++ < digit_num[0]) {
          fprintf(fp_maf, " ");
        }
        while (digit_num1[1] ++ < digit_num[1]) {
          fprintf(fp_maf, " ");
        }
        fprintf(fp_maf, " %ld", mut.seq_left - 1);
        while (digit_num1[2] ++ < digit_num[2]) {
          fprintf(fp_maf, " ");
        }
        fprintf(fp_maf, " %ld +", mut.seq_right - mut.seq_left + 1);
        while (digit_num1[3] ++ < digit_num[3]) {
          fprintf(fp_maf, " ");
        }
        fprintf(fp_maf, " %ld %s\n", ref.len, mut.maf_ref_seq);
        fprintf(fp_maf, "s %s", id); 
        while (digit_num2[0] ++ < digit_num[0]) {
          fprintf(fp_maf, " ");
        }
        while (digit_num2[1] ++ < digit_num[1]) {
          fprintf(fp_maf, " ");
        }
        fprintf(fp_maf, " %d", 0);
        while (digit_num2[2] ++ < digit_num[2]) {
          fprintf(fp_maf, " ");
        }
        fprintf(fp_maf, " %ld %c", len, mut.seq_strand);
        while (digit_num2[3] ++ < digit_num[3]) {
          fprintf(fp_maf, " ");
        }
        fprintf(fp_maf, " %ld %s\n\n", len, mut.maf_seq);
      }
    }

    sampling_num = 0;
  }

  sim.res_len_mean = (double)sim.res_len_total / sim.res_num;
  sim.res_accuracy_mean = accuracy_total / sim.res_num;

  if (sim.res_num == 1) {
    sim.res_len_sd = 0.0;
    sim.res_accuracy_sd = 0.0;
  } else {
    variance = 0.0;
    for (i=0; i<=sim.len_max; i++) {
      if (freq_len[i] > 0) {
        variance += pow((sim.res_len_mean - i), 2) * freq_len[i];
      }
    }
    sim.res_len_sd = sqrt(variance / sim.res_num);

    variance = 0.0;
    for (i=0; i<=100000; i++) {
      if (freq_accuracy[i] > 0) {
        variance += pow((sim.res_accuracy_mean - i * 0.00001), 2) * freq_accuracy[i];
      }
    }
    sim.res_accuracy_sd = sqrt(variance / sim.res_num);
  }

  return PBSUCCEEDED;
}

/////////////////////////////////////////////////////
// Function: simulate_by_model - Simulate by Model //
/////////////////////////////////////////////////////

int CPBSim::simulate_by_model(void) {
  long len;
  long long len_total = 0;
  long num;
  long i, j, k;
  double prob, mean, variance, sd;
  double len_prob_total, accuracy_prob_total, qc_prob_total, value, sum;
  double accuracy_total = 0.0;
  int accuracy;

long *pprob2len;
long *pprob2accuracy;
long *pprob2qc;

  long len_rand_value, accuracy_rand_value, qc_rand_value[101];
  long start_wk, end_wk;
  long index;
  long accuracy_min, accuracy_max;
  char id[128];
  int digit_num1[4], digit_num2[4], digit_num[4];

//  long prob2len[100001], prob2accuracy[100001], prob2qc[101][1001];
pprob2len = new long [100002];			// not too confident that pbsim code will not reference past end so allow additional
memset(pprob2len,0,sizeof(long)*100002);
pprob2accuracy = new long [100002];
memset(pprob2accuracy,0,sizeof(long)*100002);
pprob2qc = new long [102*1002];
memset(pprob2qc,0,sizeof(long)*102*1002);
memset(qc_rand_value,0,sizeof(qc_rand_value));

memset(freq_len,0,sizeof(freq_len));
memset(freq_accuracy,0,sizeof(freq_accuracy));
 
  for (i=0; i<=93; i++) {
    mut.sub_thre[i] = int((qc[i].prob * sim.sub_rate) * 1000000);
    mut.ins_thre[i] = int((qc[i].prob * (sim.sub_rate + sim.ins_rate)) * 1000000);
  }
  mut.del_thre = int((1.0 - sim.accuracy_mean) * sim.del_rate * 1000000);

  accuracy_min = (long)(sim.accuracy_min * 100);
  accuracy_max = (long)(sim.accuracy_max * 100);

  // length distribution
  variance = log(1 + pow((sim.len_sd / sim.len_mean) ,2));
  mean = log(sim.len_mean) - variance * 0.5;
  sd = sqrt(variance);

  if (sim.len_sd == 0.0) {
    pprob2len[1] = int(sim.len_mean + 0.5);
    len_rand_value = 1;
  } else {
    start_wk = 1; 
    len_prob_total = 0.0;
    for (i=sim.len_min; i<=sim.len_max; i++) {
      len_prob_total += exp(-1 * pow((log(i)-mean), 2) / 2 / variance) / sqrt(2*PBS_PI) / sd / i;
      end_wk = int(len_prob_total * 100000);
      if (end_wk > 100000) {
        end_wk = 100000;
      }

      for (j=start_wk; j<=end_wk; j++) {
        pprob2len[j] = i;
      }

      if (end_wk >= 100000) {
        break;
      }
      start_wk = end_wk + 1;
    }
    len_rand_value = end_wk;
  }

  if (len_rand_value < 1) {
    gDiagnostics.DiagOut(eDLFatal, gszProcName,  "ERROR: length parameters are not appropriate.\n");
	if(pprob2len != NULL)
		delete pprob2len;
	if(pprob2accuracy != NULL)
		delete pprob2accuracy;
	if(pprob2qc != NULL)
		delete pprob2qc;
    return PBFAILED;
  }

  // accuracy distribution
  if (sim.data_type == DATA_TYPE_CLR) {
    mean = sim.accuracy_mean * 100;
    sd = sim.accuracy_sd * 100;
    variance = pow(sd, 2);

    if (sd == 0.0) {
      pprob2accuracy[1] = int(mean + 0.5);
      accuracy_rand_value = 1;
    } else {
      start_wk = 1; 
      accuracy_prob_total = 0.0;
      for (i=accuracy_min; i<=accuracy_max; i++) {
        accuracy_prob_total += exp(-1 * pow(i - mean, 2) / 2 / variance) / sqrt(2 * PBS_PI) / sd;
        end_wk = int(accuracy_prob_total * 100000);
        if (end_wk > 100000) {
          end_wk = 100000;
        }

        for (j=start_wk; j<=end_wk; j++) {
          pprob2accuracy[j] = i;
        }

        if (end_wk >= 100000) {
          break;
        }
        start_wk = end_wk + 1;
      }
      accuracy_rand_value = end_wk;
    }
  } else {
    sum = 0;
    for (i=accuracy_min; i<=accuracy_max; i++) {
      sum += exp(0.5 * (i - 75));
    }

    start_wk = 1; 
    accuracy_prob_total = 0.0;
    for (i=accuracy_min; i<=accuracy_max; i++) {
      accuracy_prob_total += exp(0.5 * (i - 75)) / sum;
      end_wk = int(accuracy_prob_total * 100000);
      if (end_wk > 100000) {
        end_wk = 100000;
      }

      for (j=start_wk; j<=end_wk; j++) {
        pprob2accuracy[j] = i;
      }

      if (end_wk >= 100000) {
        break;
      }
      start_wk = end_wk + 1;
    }
    accuracy_rand_value = end_wk;
  }

  if (accuracy_rand_value < 1) {
    gDiagnostics.DiagOut(eDLFatal, gszProcName,  "ERROR: accuracy parameters are not appropriate.\n");
	if(pprob2len != NULL)
		delete pprob2len;
	if(pprob2accuracy != NULL)
		delete pprob2accuracy;
	if(pprob2qc != NULL)
		delete pprob2qc;
    return PBFAILED;
  }

  // quality code distributiin
  for (i=accuracy_min; i<=accuracy_max; i++) {
    start_wk = 1; 
    qc_prob_total = 0.0;

    for (j=model_qc[i].min; j<=model_qc[i].max; j++) {
      qc_prob_total += model_qc[i].prob[j];
      end_wk = int(qc_prob_total * 1000);
      if (end_wk > 1000) {
        end_wk = 1000;
      }

      for (k=start_wk; k<=end_wk; k++) {
        pprob2qc[(i*1001)+k] = j;
      }

      if (end_wk >= 1000) {
        break;
      }
      start_wk = end_wk + 1;
    }
    qc_rand_value[i] = end_wk;
  }

  // simulation
  while (len_total < sim.len_quota) {
	if(len_rand_value > 1)
		index = RGseeds.IRandom(0,len_rand_value - 1);
	else
		index = 0; 
    len = pprob2len[index];
    if (len_total + len > sim.len_quota) {
      len = (long)(sim.len_quota - len_total);

      if (len < sim.len_min) {
        len = sim.len_min;
      }
    }

    mut.tmp_len_max = len;
	if(accuracy_rand_value > 1)
		index = RGseeds.IRandom(0,accuracy_rand_value - 1);
	else
		index = 0;
    accuracy = pprob2accuracy[index];

    num = 0;
	long indexR;
    while (num < len) {
      index = 0;
	  if(qc_rand_value[accuracy] > 1)
		indexR =  RGseeds.IRandom(0,qc_rand_value[accuracy] - 1);
	  else
		indexR = 0;
	  index	= pprob2qc[(accuracy*1001)+indexR];
      mut.qc[num ++] = qc[index].character;
      if (num >= len) {
        break;
      }
    }
    mut.qc[num] = '\0';

    if (mutate() == PBFAILED) {
	if(pprob2len != NULL)
		delete pprob2len;
	if(pprob2accuracy != NULL)
		delete pprob2accuracy;
	if(pprob2qc != NULL)
		delete pprob2qc;
      return PBFAILED;
    }

    len = (long)strlen(mut.new_seq);
    sim.res_len_total += len;
    len_total += len;
    freq_len[len] ++;
    sim.res_num ++;

    if (len > sim.res_len_max) {
      sim.res_len_max = len;
    }
    if (len < sim.res_len_min) {
      sim.res_len_min = len;
    }

    prob = 0.0;
    for (i=0; i<len; i++) {
      prob += qc[(int)mut.new_qc[i] - 33].prob;
    }
    value = 1.0 - (prob / len);
	if(value < 0.0)
		value = 0.0;
	else
		if(value > 1.0)
			value = 1.0;

    accuracy_total += value;
    accuracy = (int)(value * 100000);
    freq_accuracy[accuracy] ++;

    sprintf(id, "S%ld_%ld", ref.num, sim.res_num);
    fprintf(fp_fq, "@%s\n%s\n+%s\n%s\n", id, mut.new_seq, id, mut.new_qc);

    digit_num1[0] = 3;
    digit_num2[0] = 1 + count_digit(sim.res_num);
    digit_num[0] = (digit_num1[0] >= digit_num2[0]) ? digit_num1[0] : digit_num2[0];

    digit_num1[1] = count_digit((mut.seq_left - 1));
    digit_num2[1] = 1;
    digit_num[1] = (digit_num1[1] >= digit_num2[1]) ? digit_num1[1] : digit_num2[1];

    digit_num1[2] = count_digit((mut.seq_right - mut.seq_left + 1));
    digit_num2[2] = count_digit(len);
    digit_num[2] = (digit_num1[2] >= digit_num2[2]) ? digit_num1[2] : digit_num2[2];

    digit_num1[3] = count_digit(ref.len);
    digit_num2[3] = count_digit(len);
    digit_num[3] = (digit_num1[3] >= digit_num2[3]) ? digit_num1[3] : digit_num2[3];

    fprintf(fp_maf, "a\ns ref");
    while (digit_num1[0] ++ < digit_num[0]) {
      fprintf(fp_maf, " ");
    }
    while (digit_num1[1] ++ < digit_num[1]) {
      fprintf(fp_maf, " ");
    }
    fprintf(fp_maf, " %ld", mut.seq_left - 1);
    while (digit_num1[2] ++ < digit_num[2]) {
      fprintf(fp_maf, " ");
    }
    fprintf(fp_maf, " %ld +", mut.seq_right - mut.seq_left + 1);
    while (digit_num1[3] ++ < digit_num[3]) {
      fprintf(fp_maf, " ");
    }
    fprintf(fp_maf, " %ld %s\n", ref.len, mut.maf_ref_seq);
    fprintf(fp_maf, "s %s", id);
    while (digit_num2[0] ++ < digit_num[0]) {
      fprintf(fp_maf, " ");
    }
    while (digit_num2[1] ++ < digit_num[1]) {
      fprintf(fp_maf, " ");
    }
    fprintf(fp_maf, " %d", 0);
    while (digit_num2[2] ++ < digit_num[2]) {
      fprintf(fp_maf, " ");
    }
    fprintf(fp_maf, " %ld %c", len, mut.seq_strand);
    while (digit_num2[3] ++ < digit_num[3]) {
      fprintf(fp_maf, " ");
    }
    fprintf(fp_maf, " %ld %s\n\n", len, mut.maf_seq);
  }

  sim.res_len_mean = (double)sim.res_len_total / sim.res_num;
  sim.res_accuracy_mean = accuracy_total / sim.res_num;

  if (sim.res_num == 1) {
    sim.res_len_sd = 0.0;
    sim.res_accuracy_sd = 0.0;
  } else {
    variance = 0.0;
    for (i=0; i<=sim.len_max; i++) {
      if (freq_len[i] > 0) {
        variance += pow((sim.res_len_mean - i), 2) * freq_len[i];
      }
    }
    sim.res_len_sd = sqrt(variance / sim.res_num);

    variance = 0.0;
    for (i=0; i<=100000; i++) {
      if (freq_accuracy[i] > 0) {
        variance += pow((sim.res_accuracy_mean - i * 0.00001), 2) * freq_accuracy[i];
      }
    }
    sim.res_accuracy_sd = sqrt(variance / sim.res_num);
  }
if(pprob2len != NULL)
	delete pprob2len;
if(pprob2accuracy != NULL)
	delete pprob2accuracy;
if(pprob2qc != NULL)
	delete pprob2qc;
  return PBSUCCEEDED;
}

/////////////////////////////////////////////////////////////
// Function: print_sim_param - Print simulation parameters //
/////////////////////////////////////////////////////////////

void CPBSim::print_sim_param(void) {
  gDiagnostics.DiagOut(eDLInfo, gszProcName,  ":::: Simulation parameters :::");

  if (sim.process == PROCESS_MODEL) {
    gDiagnostics.DiagOut(eDLInfo, gszProcName, "Simulated by stochastic model");
  } else {
    gDiagnostics.DiagOut(eDLInfo, gszProcName, "Simulated by fastq sampling");
  }

  gDiagnostics.DiagOut(eDLInfo, gszProcName, "prefix : %s", sim.prefix);
  if (sim.set_flg[14]) {
    gDiagnostics.DiagOut(eDLInfo, gszProcName, "sample_profile_id : %s", sim.profile_id);
  }

  if (sim.data_type == DATA_TYPE_CLR) {
    gDiagnostics.DiagOut(eDLInfo, gszProcName, "data-type : CLR");
  } else {
    gDiagnostics.DiagOut(eDLInfo, gszProcName, "data-type : CCS");
  }

  gDiagnostics.DiagOut(eDLInfo, gszProcName, "depth : %lf", sim.depth);

  if (sim.set_flg[0]) {
    gDiagnostics.DiagOut(eDLInfo, gszProcName, "length-mean : (sampling FASTQ)");
    gDiagnostics.DiagOut(eDLInfo, gszProcName, "length-sd : (sampling FASTQ)");
  } else {
    gDiagnostics.DiagOut(eDLInfo, gszProcName, "length-mean : %f", sim.len_mean);
    gDiagnostics.DiagOut(eDLInfo, gszProcName, "length-sd : %f", sim.len_sd);
  }
  gDiagnostics.DiagOut(eDLInfo, gszProcName,  "length-min : %ld", sim.len_min);
  gDiagnostics.DiagOut(eDLInfo, gszProcName, "length-max : %ld", sim.len_max);

  if (sim.set_flg[0]) {
    gDiagnostics.DiagOut(eDLInfo, gszProcName,  "accuracy-mean : (sampling FASTQ)");
    gDiagnostics.DiagOut(eDLInfo, gszProcName,  "accuracy-sd : (sampling FASTQ)");
  } else {
    gDiagnostics.DiagOut(eDLInfo, gszProcName, "accuracy-mean : %f", sim.accuracy_mean);
    gDiagnostics.DiagOut(eDLInfo, gszProcName,  "accuracy-sd : %f", sim.accuracy_sd);
  }
  gDiagnostics.DiagOut(eDLInfo, gszProcName, "accuracy-min : %f", sim.accuracy_min);
  gDiagnostics.DiagOut(eDLInfo, gszProcName, "accuracy-max : %f", sim.accuracy_max);

 gDiagnostics.DiagOut(eDLInfo, gszProcName, "difference-ratio : %d:%d:%d",
    sim.sub_ratio, sim.ins_ratio, sim.del_ratio);
}

/////////////////////////////////////////////////////////////////
// Function: set_mut - Set mutation parameters and varianeces  //
/////////////////////////////////////////////////////////////////

int CPBSim::set_mut(void) {
  mut.sub_nt_a = (char *)"TGC";
  mut.sub_nt_t = (char *)"AGC";
  mut.sub_nt_g = (char *)"ATC";
  mut.sub_nt_c = (char *)"ATG";
  mut.sub_nt_n = (char *)"ATGC";
  mut.ins_nt = (char *)"ATGC";

  if ((mut.qc = (char *)malloc(sim.len_max + 10)) == 0) {
    gDiagnostics.DiagOut(eDLFatal, gszProcName, "ERROR: Cannot allocate memory");
    return PBFAILED;
  }

  if ((mut.new_qc = (char *)malloc(sim.len_max * 2 + 10)) == 0) {
    gDiagnostics.DiagOut(eDLFatal, gszProcName,"ERROR: Cannot allocate memory");
    return PBFAILED;
  }

  if ((mut.tmp_qc = (char *)malloc(sim.len_max * 2 + 10)) == 0) {
    gDiagnostics.DiagOut(eDLFatal, gszProcName,"ERROR: Cannot allocate memory");
    return PBFAILED;
  }

  if ((mut.seq = (char *)malloc(sim.len_max * 2 + 10)) == 0) {
    gDiagnostics.DiagOut(eDLFatal, gszProcName,"ERROR: Cannot allocate memory");
    return PBFAILED;
  }

  if ((mut.new_seq = (char *)malloc(sim.len_max * 2 + 10)) == 0) {
    gDiagnostics.DiagOut(eDLFatal, gszProcName, "ERROR: Cannot allocate memory");
    return PBFAILED;
  }

  if ((mut.maf_seq = (char *)malloc(sim.len_max * 2 + 10)) == 0) {
    gDiagnostics.DiagOut(eDLFatal, gszProcName,"ERROR: Cannot allocate memory");
    return PBFAILED;
  }

  if ((mut.maf_ref_seq = (char *)malloc(sim.len_max * 2 + 10)) == 0) {
    gDiagnostics.DiagOut(eDLFatal, gszProcName, "ERROR: Cannot allocate memory");
    return PBFAILED;
  }

  return PBSUCCEEDED;
}

////////////////////////////////////
// Function: mutate - Mutate read //
////////////////////////////////////

int CPBSim::mutate(void) {
  char nt;
  long i;
  long index;
  long rand_value;
  long qc_value;
  long len;
  long offset, seq_offset, maf_offset;

  len = (long)strlen(mut.qc);
  if (mut.tmp_len_max < len) {
    len = mut.tmp_len_max;
  }

  // Place deletions
  offset = 0;
  for (i=0; i<len-1; i++) {
    mut.tmp_qc[offset ++] = mut.qc[i];
    if ((rand_value =  RGseeds.IRandom(0,1000000-1)) < mut.del_thre) {
      mut.tmp_qc[offset ++] = ' ';
      sim.res_del_num ++;
    }
  }
  mut.tmp_qc[offset ++] = mut.qc[len - 1];
  mut.tmp_qc[offset] = '\0';

  len = (long)strlen(mut.tmp_qc);

  if (len >= ref.len) {
    offset = 0;
    len = ref.len;
  } else {
    offset =  RGseeds.IRandom(0,(ref.len - len));
  }

  mut.seq_left = offset + 1;
  mut.seq_right = offset + len;

  if (sim.res_num % 2 == 0) {
    mut.seq_strand = '+';

    for (i=0; i<len; i++) {
      nt = toupper(ref.seq[offset + i]);
      mut.seq[i] = nt;
    }
  } else {
    mut.seq_strand = '-';

    for (i=0; i<len; i++) {
      nt = toupper(ref.seq[offset + i]);

      if (nt == 'A') {
        mut.seq[len-1-i] = 'T';
      } else if (nt == 'T') {
        mut.seq[len-1-i] = 'A';
      } else if (nt == 'G') {
        mut.seq[len-1-i] = 'C';
      } else if (nt == 'C') {
        mut.seq[len-1-i] = 'G';
      } else {
        mut.seq[len-1-i] = nt;
      }
    }
  }
  mut.seq[len] = '\0';

  // Place substitutions and insertions
  offset = 0;
  seq_offset = 0;
  maf_offset = 0;
  for (i=0; i<len; i++) {
    nt = mut.seq[seq_offset ++];

    if (mut.tmp_qc[i] == ' ') {
      mut.maf_seq[maf_offset] = '-';
      mut.maf_ref_seq[maf_offset] = nt;
      maf_offset ++;
      continue;
    }

    mut.new_qc[offset] = mut.tmp_qc[i];

    rand_value = RGseeds.IRandom(0, 1000000-1);
    qc_value = (int)mut.tmp_qc[i] - 33;

    if (rand_value < mut.sub_thre[qc_value]) {
      sim.res_sub_num ++;
      index =  RGseeds.IRandom(0,2);
      if (nt == 'A') {
        mut.new_seq[offset] = mut.sub_nt_a[index];
      } else if (nt == 'T') {
        mut.new_seq[offset] = mut.sub_nt_t[index];
      } else if (nt == 'G') {
        mut.new_seq[offset] = mut.sub_nt_g[index];
      } else if (nt == 'C') {
        mut.new_seq[offset] = mut.sub_nt_c[index];
      } else {
        index = RGseeds.IRandom(0,3);
        mut.new_seq[offset] = mut.sub_nt_n[index];
      }
      mut.maf_ref_seq[maf_offset] = nt;
    } else if (rand_value < mut.ins_thre[qc_value]) {
      sim.res_ins_num ++;
      index = RGseeds.IRandom(0, 7);
      if (index >= 4) {
        mut.new_seq[offset] = nt;
      } else {
        mut.new_seq[offset] = mut.ins_nt[index];
      }
      seq_offset --;
      if (mut.seq_strand == '+') {
        mut.seq_right --;
      } else {
        mut.seq_left ++;
      }
      mut.maf_ref_seq[maf_offset] = '-';
    } else {
      mut.new_seq[offset] = nt;
      mut.maf_ref_seq[maf_offset] = nt;
    }
    mut.maf_seq[maf_offset] = mut.new_seq[offset];
    maf_offset ++;
    offset ++;
  }
  mut.new_qc[offset] = '\0';
  mut.new_seq[offset] = '\0';
  mut.maf_seq[maf_offset] = '\0';
  mut.maf_ref_seq[maf_offset] = '\0';

  if (mut.seq_strand == '-') {
    revcomp(mut.maf_seq);
    revcomp(mut.maf_ref_seq);
  }

  return PBSUCCEEDED;
}

////////////////////////////////////////////////////////////////
// Function: print_simulation_stats - Print Simulation Stats. //
////////////////////////////////////////////////////////////////

void CPBSim::print_simulation_stats(void) {
  sim.res_depth = (double)sim.res_len_total / ref.len;
  sim.res_sub_rate = (double)sim.res_sub_num / sim.res_len_total;
  sim.res_ins_rate = (double)sim.res_ins_num / sim.res_len_total;
  sim.res_del_rate = (double)sim.res_del_num / sim.res_len_total;

  gDiagnostics.DiagOut(eDLInfo, gszProcName, ":::: Simulation stats (ref.%d) ::::", ref.num);
  gDiagnostics.DiagOut(eDLInfo, gszProcName, "read num. : %ld", sim.res_num);
  gDiagnostics.DiagOut(eDLInfo, gszProcName,"depth : %lf", sim.res_depth);
  gDiagnostics.DiagOut(eDLInfo, gszProcName, "read length mean (SD) : %f (%f)",
    sim.res_len_mean, sim.res_len_sd);
  gDiagnostics.DiagOut(eDLInfo, gszProcName, "read length min : %ld", sim.res_len_min);
  gDiagnostics.DiagOut(eDLInfo, gszProcName, "read length max : %ld", sim.res_len_max);
  gDiagnostics.DiagOut(eDLInfo, gszProcName,"read accuracy mean (SD) : %f (%f)",
    sim.res_accuracy_mean, sim.res_accuracy_sd);
  gDiagnostics.DiagOut(eDLInfo, gszProcName,"substitution rate. : %f", sim.res_sub_rate);
  gDiagnostics.DiagOut(eDLInfo, gszProcName, "insertion rate. : %f", sim.res_ins_rate);
  gDiagnostics.DiagOut(eDLInfo, gszProcName,"deletion rate. : %f", sim.res_del_rate);
 }

///////////////////////////////////////////////////////
// Function: set_model_qc - Set quality code model   //
///////////////////////////////////////////////////////

int CPBSim::set_model_qc(void) {
  FILE *fp;
  char line[BUF_SIZE+1];
  char *tp;
  long accuracy;
  int num;
  int i, j;

  if ((fp = fopen(sim.model_qc_file, "r")) == NULL) {
    gDiagnostics.DiagOut(eDLFatal, gszProcName, "ERROR: Cannot open file: %s", sim.model_qc_file);
    return PBFAILED;
  }

  for (i=0; i<=100; i++) {
    for (j=0; j<=93; j++) {
      model_qc[i].prob[j] = 0.0;
    }
  }

  while (fgets(line, BUF_SIZE, fp) != NULL) {
    trim(line);

    tp = strtok(line, "\t");
    accuracy = atoi(tp);

    num = 0;
    tp = strtok(NULL, "\t");
    while (tp != NULL) {
      model_qc[accuracy].prob[num] = atof(tp);
      num ++;
      tp = strtok(NULL, "\t");
    }
  }
  fclose(fp);

  for (i=0; i<=100; i++) {
    model_qc[i].min = 0;
    model_qc[i].max = 93;

    for (j=0; j<=93; j++) {
      if (model_qc[i].prob[j] > 0.0) {
        model_qc[i].min = j;
        break;
      }
    }
    for (j=93; j>=0; j--) {
      if (model_qc[i].prob[j] > 0.0) {
        model_qc[i].max = j;
        break;
      }
    }
  }

  return PBSUCCEEDED;
}


///////////////////////////////////////////////////////
// Function: get_time - Get time                     //
///////////////////////////////////////////////////////

long CPBSim::get_time(void) {
time_t aclock;
time(&aclock );   // Get time in seconds
return (long)aclock;
}

///////////////////////////////////////
// Function: print_help - Print help //
///////////////////////////////////////

void PBSim_print_help(void) {
  fprintf(stderr, "\n");
  fprintf(stderr, "USAGE: pbsim [options] <reference>\n\n");
  fprintf(stderr, " <reference>           FASTA format file.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, " [general options]\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  --prefix             prefix of output files (sd).\n");
  fprintf(stderr, "  --data-type          data type. CLR or CCS (CLR).\n");
  fprintf(stderr, "  --depth              depth of coverage (CLR: 20.0, CCS: 50.0).\n");
  fprintf(stderr, "  --length-min         minimum length (100).\n");
  fprintf(stderr, "  --length-max         maximum length (CLR: 25000, CCS: 2500).\n");
  fprintf(stderr, "  --accuracy-min       minimum accuracy.\n");
  fprintf(stderr, "                       (CLR: 0.75, CCS: fixed as 0.75).\n");
  fprintf(stderr, "                       this option can be used only in case of CLR.\n");
  fprintf(stderr, "  --accuracy-max       maximum accuracy.\n");
  fprintf(stderr, "                       (CLR: 1.00, CCS: fixed as 1.00).\n");
  fprintf(stderr, "                       this option can be used only in case of CLR.\n");
  fprintf(stderr, "  --difference-ratio   ratio of differences. substitution:insertion:deletion.\n");
  fprintf(stderr, "                       each value up to 1000 (CLR: 10:60:30, CCS:6:21:73).\n");
  fprintf(stderr, "  --seed               for a pseudorandom number generator (Unix time).\n");
  fprintf(stderr, "\n");
  fprintf(stderr, " [options of sampling-based simulation]\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  --sample-fastq       FASTQ format file to sample.\n");
  fprintf(stderr, "  --sample-profile-id  sample-fastq (filtered) profile ID.\n");
  fprintf(stderr, "                       when using --sample-fastq, profile is stored.\n");
  fprintf(stderr, "                       'sample_profile_<ID>.fastq', and\n");
  fprintf(stderr, "                       'sample_profile_<ID>.stats' are created.\n");
  fprintf(stderr, "                       when not using --sample-fastq, profile is re-used.\n");
  fprintf(stderr, "                       Note that when profile is used, --length-min,max,\n");
  fprintf(stderr, "                       --accuracy-min,max would be the same as the profile.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, " [options of model-based simulation].\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  --model_qc           model of quality code.\n");
  fprintf(stderr, "  --length-mean        mean of length model (CLR: 3000.0, CCS:450.0).\n");
  fprintf(stderr, "  --length-sd          standard deviation of length model.\n");
  fprintf(stderr, "                       (CLR: 2300.0, CCS: 170.0).\n");
  fprintf(stderr, "  --accuracy-mean      mean of accuracy model.\n");
  fprintf(stderr, "                       (CLR: 0.78, CCS: fixed as 0.98).\n");
  fprintf(stderr, "                       this option can be used only in case of CLR.\n");
  fprintf(stderr, "  --accuracy-sd        standard deviation of accuracy model.\n");
  fprintf(stderr, "                       (CLR: 0.02, CCS: fixed as 0.02).\n");
  fprintf(stderr, "                       this option can be used only in case of CLR.\n");
  fprintf(stderr, "\n");
}

/////////////////////////////////////////
// Function: count_digit - count digit //
/////////////////////////////////////////

int CPBSim::count_digit(long num) {
  int digit = 1;
  int quotient;

  quotient = int(num / 10);

  while (quotient != 0) {
    digit ++;
    quotient = int(quotient / 10);
  }

  return digit;
}  

//////////////////////////////////////////////////////
// Function: revcomp - convert to reverse complement//
//////////////////////////////////////////////////////

void CPBSim::revcomp(char* str) {
  int i, len;
  char c;

  len = (long)strlen(str);

  for(i=0; i<len/2; i++) {
    c = str[i];
    str[i] = str[len-i-1];
    str[len-i-1] = c;
  }

  for(i=0; i<len; i++) {
    if (str[i] == 'A') {
      str[i] = 'T';
    } else if (str[i] == 'T') {
      str[i] = 'A';
    } else if (str[i] == 'G') {
      str[i] = 'C';
    } else if (str[i] == 'C') {
      str[i] = 'G';
    }
  }
}  
