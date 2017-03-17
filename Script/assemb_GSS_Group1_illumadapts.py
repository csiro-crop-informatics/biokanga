#!/usr/bin/env python
# Full processing lifecycle, from raw reads quality checking through to scaffolding, of a targeted set of NGS readsets
# This version is trimming reads containing adaptor sequences 
# This version is trimming the 5' end of reads by 0bp 
# Scaffolding is allowing 0 subs
# This version allows for the option of reassembling contigs only with a higher error rate of 5
# first phase is 'biokanga ngsqc'
# next phase is 'biokangs filter'
# next phase is 'biokanga assemb'
# final phase is 'biokanga scaffold'

from __future__ import print_function
from subprocess import *
import os
import platform
import errno
import logging

# specify the sample sets to be processed as a dictionary with the filenames as the values
# if PE then ensure the files are ordered as PE1 followed by PE2
# the resultant resultsets are prefixed with the dictionary name

SampleSets= {"GSS_1A_1B_1D" :
        ["1AL_001_5_1_sequence.fastq.gz","1AL_001_5_2_sequence.fastq.gz",
        "1AL_002_2_1_sequence.fastq.gz","1AL_002_2_2_sequence.fastq.gz",
        "1AL_003_4_1_sequence.fastq.gz","1AL_003_4_2_sequence.fastq.gz",
        "1AL_004_3_1_sequence.fastq.gz","1AL_004_3_2_sequence.fastq.gz",
        "1AL_005_4_1_sequence.fastq.gz","1AL_005_4_2_sequence.fastq.gz",
        "1AL_006_L003_R1.fastq.gz","1AL_006_L003_R2.fastq.gz",
        "1AL_007_L004_R1.fastq.gz","1AL_007_L004_R2.fastq.gz",
        "1AS_001_1_1_sequence.fastq.gz","1AS_001_1_2_sequence.fastq.gz",
        "1AS_002_2_1_sequence.fastq.gz","1AS_002_2_2_sequence.fastq.gz",
        "1AS_003_L001_R1.fastq.gz","1AS_003_L001_R2.fastq.gz",
        "1AS_004_L002_R1.fastq.gz","1AS_004_L002_R2.fastq.gz",
        "1BL_001_7_1_sequence.fastq.gz","1BL_001_7_2_sequence.fastq.gz",
        "1BL_002_3_1_sequence.fastq.gz","1BL_002_3_2_sequence.fastq.gz",
        "1BL_003_4_1_sequence.fastq.gz","1BL_003_4_2_sequence.fastq.gz",
        "1BL_004_L001_R1.fastq.gz","1BL_004_L001_R2.fastq.gz",
        "1BS_001_2_1_sequence.fastq.gz","1BS_001_2_2_sequence.fastq.gz",
        "1BS_002_3_1_sequence.fastq.gz","1BS_002_3_2_sequence.fastq.gz",
        "1BS_003_6_1_sequence.fastq.gz","1BS_003_6_2_sequence.fastq.gz",
        "1DL_001_L006_R1.fastq.gz","1DL_001_L006_R2.fastq.gz",
        "1DL_002_L007_R1.fastq.gz","1DL_002_L007_R2.fastq.gz",
        "1DL_003_L008_R1.fastq.gz","1DL_003_L008_R2.fastq.gz",
        "1DS_001_6_1_sequence.fastq.gz","1DS_001_6_2_sequence.fastq.gz",
        "1DS_002_6_1_sequence.fastq.gz","1DS_002_6_2_sequence.fastq.gz",
        "1DS_003_7_1_sequence.fastq.gz","1DS_003_7_2_sequence.fastq.gz"],
    "GSS_1AL" :
        ["1AL_001_5_1_sequence.fastq.gz","1AL_001_5_2_sequence.fastq.gz",
        "1AL_002_2_1_sequence.fastq.gz","1AL_002_2_2_sequence.fastq.gz",
        "1AL_003_4_1_sequence.fastq.gz","1AL_003_4_2_sequence.fastq.gz",
        "1AL_004_3_1_sequence.fastq.gz","1AL_004_3_2_sequence.fastq.gz",
        "1AL_005_4_1_sequence.fastq.gz","1AL_005_4_2_sequence.fastq.gz",
        "1AL_006_L003_R1.fastq.gz","1AL_006_L003_R2.fastq.gz",
        "1AL_007_L004_R1.fastq.gz","1AL_007_L004_R2.fastq.gz"],
    "GSS_1AS" :
        ["1AS_001_1_1_sequence.fastq.gz","1AS_001_1_2_sequence.fastq.gz",
        "1AS_002_2_1_sequence.fastq.gz","1AS_002_2_2_sequence.fastq.gz",
        "1AS_003_L001_R1.fastq.gz","1AS_003_L001_R2.fastq.gz",
        "1AS_004_L002_R1.fastq.gz","1AS_004_L002_R2.fastq.gz"],
    "GSS_1BL" :
        ["1BL_001_7_1_sequence.fastq.gz","1BL_001_7_2_sequence.fastq.gz",
        "1BL_002_3_1_sequence.fastq.gz","1BL_002_3_2_sequence.fastq.gz",
        "1BL_003_4_1_sequence.fastq.gz","1BL_003_4_2_sequence.fastq.gz",
        "1BL_004_L001_R1.fastq.gz","1BL_004_L001_R2.fastq.gz"],
    "GSS_1BS" :
        ["1BS_001_2_1_sequence.fastq.gz","1BS_001_2_2_sequence.fastq.gz",
        "1BS_002_3_1_sequence.fastq.gz","1BS_002_3_2_sequence.fastq.gz",
        "1BS_003_6_1_sequence.fastq.gz","1BS_003_6_2_sequence.fastq.gz"],
    "GSS_1DL" :
        ["1DL_001_L006_R1.fastq.gz","1DL_001_L006_R2.fastq.gz",
        "1DL_002_L007_R1.fastq.gz","1DL_002_L007_R2.fastq.gz",
        "1DL_003_L008_R1.fastq.gz","1DL_003_L008_R2.fastq.gz"],
    "GSS_1DS" :
        ["1DS_001_6_1_sequence.fastq.gz","1DS_001_6_2_sequence.fastq.gz",
        "1DS_002_6_1_sequence.fastq.gz","1DS_002_6_2_sequence.fastq.gz",
        "1DS_003_7_1_sequence.fastq.gz","1DS_003_7_2_sequence.fastq.gz"]}

# the resultant resultsets are prefixed with the dictionary name

Trim5 = "0"      # trim 5' ends of reads by 0bp
QCMinPhred = "0"  # when QC'ing then do not accept reads with Phred score at any base offset with Phred less than this
ScaffMinInsertSize = "110"    # when scaffolding then need to know the range of insert sizes for the scaffolding reads
ScaffMaxInsertSize = "1500"   # in this experiment, the filtered PE's used in the assembly are being used for scaffolding
 
# filenames which have following suffixes are assumed to be PE1 of a paired end
PE1Sfxs=["_1_sequence.fastq.gz","_R1.fastq.gz"]
# filenames which have following suffixes are assumed to be PE2 of a paired end
PE2Sfxs=["_2_sequence.fastq.gz","_R2.fastq.gz"]

# specify which processing phases are to be executed - True to execute, False to skip
bDoNGSqc = True          # set True for reads quality check reporting with 'biokanga ngsqc'
bDoFilter = True        # set True for filtering with 'biokanga filter' processing
bDoFilterAdaptors = True # set True for trimming adaptors from reads in filtering phase
bDoAssemb = True         # set True for de Novo assembly of filtered reads with 'biokanga filter'
bDoScaffold = True      # set True for scaffolding of de Novo assembled contigs with 'biokanga scaffold'
bDoReassembly = True     # set True to repeat the assembly and scaffolding phase with assembly of contigs only allowing higher mismatch rates
 
# max number of processing threads to be utilised
NumThreads = "48"

#make all directories on the given path and don't treat as error if any path components already exist
def mkdirs(newdir,mode=0777):
    try: os.makedirs(newdir,mode)
    except OSError, err:
        if(err.errno != errno.EEXIST or not os.path.isdir(newdir)):
            raise

#platform dependency on the rootpath
if(platform.platform(terse=1)[:3] != 'Win'):
    RootPath="/data1/Wheat"
else:
    RootPath="c:\\wheat"

#rootpath now known, construct dependent paths
#note that the individual phases write their primary outputs into directories with naming convention representing the phase
ExperimentsPath=RootPath+"/GSS_AssembReads"
SrcReadsPath=ExperimentsPath+"/Data"                   # directory containing the raw reads
NGSqcPath=ExperimentsPath+"/NGSqc"                         # directory to write ngsqc result files into
FilterPath=ExperimentsPath+"/Filtered"                     # directory to write filtered result files into
AssembPath=ExperimentsPath+"/Assembled"                    # directory to write assembled contig result files into
ScaffoldPath=ExperimentsPath+"/Scaffolded"                 # directory to write scaffolded assembly result files into
LogsPath=ExperimentsPath+"/Logs"
AdaptorSeqs=ExperimentsPath+"/Resources/IlluminaAdaptersAllpXENpQ108_nodup.fasta"
ReassembAssembPath=AssembPath+"Reassemb"                    # directory to write reassembled contig result files into
ReassembScaffoldPath=ScaffoldPath+"Reassemb"               # directory to write scaffolded reassembly result files into

if(os.path.isfile(LogsPath)):
    print("LogsPath '%s' exists but is a file" % LogsPath)
    exit(1)
if(not os.path.isdir(LogsPath)):
    print("LogsPath '%s' does not exist, creating path" % LogsPath)
    mkdirs(LogsPath)

thisScript = os.path.basename(__file__)
logger = logging.getLogger(thisScript)	# get logging instance using this script name as the identifier
logger.setLevel(logging.DEBUG)
# create file (all levels) and console handler (info and above) to log script generated messages
fh = logging.FileHandler(LogsPath+"/"+thisScript+".log")
fh.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
# create formatter and add it to the handlers
formatter = logging.Formatter('%(asctime)s , %(name)s , %(levelname)s , %(message)s')
fh.setFormatter(formatter)
ch.setFormatter(formatter)
# add the handlers to the logger
logger.addHandler(fh)
logger.addHandler(ch)
logger.info("Startup")


#check that dependent paths do actuallly exist and are directories, not files
if(os.path.isfile(ExperimentsPath)):
    print("ExperimentsPath '%s' exists but is a file" % ExperimentsPath)
    exit(1)
if(not os.path.isdir(ExperimentsPath)):
    print("ExperimentsPath '%s' does not exist, creating path" % ExperimentsPath)
    mkdirs(ExperimentsPath)

if(os.path.isfile(NGSqcPath)):
    print("NGSqcPath '%s' exists but is a file" % NGSqcPath)
    exit(1)
if(not os.path.isdir(NGSqcPath)):
    print("NGSqcPath '%s' does not exist, creating path" % NGSqcPath)
    mkdirs(NGSqcPath)

if(bDoFilterAdaptors and not os.path.isfile(AdaptorSeqs)):
    print("AdaptorSeqs '%s' is not a file" % AdaptorSeqs)
    exit(1)

if(os.path.isfile(FilterPath)):
    print("FilterPath '%s' exists but is a file" % FilterPath)
    exit(1)
if(not os.path.isdir(FilterPath)):
    print("FilterPath '%s' does not exist, creating path" % FilterPath)
    mkdirs(FilterPath)

if(os.path.isfile(AssembPath)):
    print("AssembPath '%s' exists but is a file" % AssembPath)
    exit(1)
if(not os.path.isdir(AssembPath)):
    print("AssembPath '%s' does not exist, creating path" % AssembPath)
    mkdirs(AssembPath)

if(os.path.isfile(ScaffoldPath)):
    print("ScaffoldPath '%s' exists but is a file" % ScaffoldPath)
    exit(1)
if(not os.path.isdir(ScaffoldPath)):
    print("ScaffoldPath '%s' does not exist, creating path" % ScaffoldPath)
    mkdirs(ScaffoldPath)

if(os.path.isfile(ReassembAssembPath)):
    print("ReassembAssembPath '%s' exists but is a file" % ReassembAssembPath)
    exit(1)
if(not os.path.isdir(ReassembAssembPath)):
    print("ReassembAssembPath '%s' does not exist, creating path" % ReassembAssembPath)
    mkdirs(ReassembAssembPath)

if(os.path.isfile(ReassembScaffoldPath)):
    print("ReassembScaffoldPath '%s' exists but is a file" % ReassembScaffoldPath)
    exit(1)
if(not os.path.isdir(ReassembScaffoldPath)):
    print("ReassembScaffoldPath '%s' does not exist, creating path" % ReassembScaffoldPath)
    mkdirs(ReassembScaffoldPath)

if(os.path.isfile(SrcReadsPath)):
    print("SrcReadsPath '%s' exists but is a file" % SrcReadsPath)
    exit(1)
if(not os.path.isdir(SrcReadsPath)):
    print("SrcReadsPath '%s' does not exist" % SrcReadsPath)
    exit(1)

if(os.path.isfile(LogsPath)):
    print("LogsPath '%s' exists but is a file" % LogsPath)
    exit(1)
if(not os.path.isdir(LogsPath)):
    print("LogsPath '%s' does not exist, creating path" % LogsPath)
    mkdirs(LogsPath)

# better to check now if source input files exist now rather than aborting processing because a specific file is missing after perhaps days of processing...

if(not os.path.isdir(SrcReadsPath)):
    print("SrcReadsPath '%s' does not exist or is a file" % SrcReadsPath)
    exit(1)
for SampleKey, ReadsFiles in SampleSets.iteritems():
    for ReadsFile in ReadsFiles:
        ReadsFileName=SrcReadsPath+"/"+ReadsFile
        if(not os.path.isfile(ReadsFileName)):
            print("Reads file '%s' for sample set '%s' does not exist or is a directory" % (ReadsFileName,SampleKey))
            exit(1)

# Seems that all raw read files do exist, start processing....

# 1st phase is for quality checks on the raw reads
if(bDoNGSqc):  
    for SampleKey, ReadsFiles in SampleSets.iteritems():
        logger.info("Starting to quality check reads in '%s'" % (SampleKey))
        Process2Exe = "biokanga"
        ParamList = [Process2Exe,
                         "ngsqc",
                         "-T"+NumThreads,
                         "-p"+QCMinPhred,
                         "-F"+LogsPath+"/"+Process2Exe+"_"+SampleKey+".log",
                         "-o"+NGSqcPath+"/"+SampleKey+".ngsqc"]
        if (bDoFilterAdaptors):
            ParamList.append("-c"+AdaptorSeqs)
        for ReadsFile in ReadsFiles:
                # if file name suffix does not match _R2.fasta.gz then assume it's a PE1
            if ReadsFile.endswith(tuple(PE1Sfxs)):
                ReadsFileNamePE1=SrcReadsPath+"/"+ReadsFile
                ParamList.append("-i"+ReadsFileNamePE1)
            else:
                if ReadsFile.endswith(tuple(PE2Sfxs)):
                   ReadsFileNamePE2=SrcReadsPath+"/"+ReadsFile
                   ParamList.append("-u"+ReadsFileNamePE2)
                else:
                   logger.info("Input reads file '%s' does not have PE1 (%s) or PE2 (%s) suffix" % (ReadsFile,PE1Sfx,PE2Sfx))
                   exit(1)
        logger.debug(ParamList)
        Process = Popen(ParamList,shell=False)
        sts = Process.wait()
        logger.info("Completed quality checking reads for '%s'" % (SampleKey))
    logger.info("All quality checking completed")

# next phase is to filter out adaptor sequences, also trim 5' reads by Trim5
if(bDoFilter): 
    for SampleKey, ReadsFiles in SampleSets.iteritems():
        logger.info("Starting to filter PE reads in '%s'" % (SampleKey))
        Process2Exe = "biokanga"
        ParamList = [Process2Exe,
                         "filter",
                         "-T"+NumThreads,
                         "-x"+Trim5,
                         "-F"+LogsPath+"/"+Process2Exe+"_"+SampleKey+".log",
                         "-O"+FilterPath+"/"+SampleKey+".x"+Trim5+".dists.csv",
                         "-o"+FilterPath+"/"+SampleKey+".x"+Trim5+".Filtered"]
        if (bDoFilterAdaptors):
            ParamList.append("-c"+AdaptorSeqs)
        for ReadsFile in ReadsFiles:
                # if file name suffix does not match _R2.fasta.gz then assume it's a PE1
            if ReadsFile.endswith(tuple(PE1Sfxs)):
               ReadsFileNamePE1=SrcReadsPath+"/"+ReadsFile
               ParamList.append("-i"+ReadsFileNamePE1)
            else:
               if ReadsFile.endswith(tuple(PE2Sfxs)):
                  ReadsFileNamePE2=SrcReadsPath+"/"+ReadsFile
                  ParamList.append("-I"+ReadsFileNamePE2)
               else:
                  logger.info("Input reads file '%s' does not have PE1 (%s) or PE2 (%s) suffix" % (ReadsFile,PE1Sfx,PE2Sfx))
                  exit(1)
        logger.debug(ParamList)
        Process = Popen(ParamList,shell=False)
        sts = Process.wait()
        logger.info("Filtering '%s' completed" % (SampleKey))
    logger.info("All filtering completed")

# next phase is to de Novo assemble
if(bDoAssemb): 
    for SampleKey, ReadsFiles in SampleSets.iteritems():
        logger.info("Starting to assemble filtered PE reads from '%s'" % (SampleKey))
        Process2Exe = "biokanga"
        ParamList = [Process2Exe,
                         "assemb",
                         "-T"+NumThreads,
                         "-F"+LogsPath+"/"+Process2Exe+"_"+SampleKey+".log",
                         "-a"+FilterPath+"/"+SampleKey+".x"+Trim5+".Filtered.R1.fasta",
                         "-A"+FilterPath+"/"+SampleKey+".x"+Trim5+".Filtered.R2.fasta",
                         "-o"+AssembPath+"/"+SampleKey+".x"+Trim5+".Filtered.assemb"]
        logger.debug(ParamList)
        Process = Popen(ParamList,shell=False)
        sts = Process.wait()
        logger.info("Assembly of '%s' completed" % (SampleKey))
    logger.info("All assembly completed")

# next phase is to scaffold the assembled contigs
if(bDoScaffold): 
    for SampleKey, ReadsFiles in SampleSets.iteritems():
        logger.info("Starting to scaffold contigs from '%s'" % (SampleKey))
        Process2Exe = "biokanga"
        ParamList = [Process2Exe,
                         "scaffold",
			 "-s0",
                         "-T"+NumThreads,
                         "-p"+ScaffMinInsertSize,
                         "-P"+ScaffMaxInsertSize,
                         "-F"+LogsPath+"/"+Process2Exe+"_"+SampleKey+".log",
                         "-a"+FilterPath+"/"+SampleKey+".x"+Trim5+".Filtered.R1.fasta",
                         "-A"+FilterPath+"/"+SampleKey+".x"+Trim5+".Filtered.R2.fasta",
                         "-c"+AssembPath+"/"+SampleKey+".x"+Trim5+".Filtered.assemb.SE.fasta",
                         "-o"+ScaffoldPath+"/"+SampleKey+".x"+Trim5+".Filtered.assemb.scaffold"]
        logger.debug(ParamList)
        Process = Popen(ParamList,shell=False)
        sts = Process.wait()
        logger.info("Scaffolding of '%s' completed" % (SampleKey))
    logger.info("All Scaffolding completed")

#perhaps the user has requested a reassembly of the contigs only but allowing a higher substitution rate
# next phase is to de Novo assemble
if(bDoReassembly):
    for SampleKey, ReadsFiles  in SampleSets.iteritems():
        logger.info("Starting to reassemble contigs previously generated for '%s'" % (SampleKey))
        Process2Exe = "biokanga"
        ParamList = [Process2Exe,
                         "assemb",
                         "-T"+NumThreads,
                         "-F"+LogsPath+"/"+Process2Exe+"_"+SampleKey+".log",
                         "-s5", 
                         "-c"+AssembPath+"/"+SampleKey+".x"+Trim5+".Filtered.assemb.SE.fasta",
                         "-o"+ReassembAssembPath+"/"+SampleKey+".x"+Trim5+".Filtered.assemb"]
        logger.debug(ParamList)
        Process = Popen(ParamList,shell=False)
        sts = Process.wait()
        logger.info("Reassembly of '%s' completed" % (SampleKey))
    logger.info("All assembly completed")

# next phase is to scaffold the assembled contigs
if(bDoReassembly):
    for SampleKey, ReadsFiles in SampleSets.iteritems():
        logger.info("Starting to scaffold reassembled contigs from '%s'" % (SampleKey))
        Process2Exe = "biokanga"
        ParamList = [Process2Exe,
                         "scaffold",
                         "-s0",
                         "-T"+NumThreads,
                         "-p"+ScaffMinInsertSize,
                         "-P"+ScaffMaxInsertSize,
                         "-F"+LogsPath+"/"+Process2Exe+"_"+SampleKey+".log",
                         "-a"+FilterPath+"/"+SampleKey+".x"+Trim5+".Filtered.R1.fasta",
                         "-A"+FilterPath+"/"+SampleKey+".x"+Trim5+".Filtered.R2.fasta",
                         "-c"+ReassembAssembPath+"/"+SampleKey+".x"+Trim5+".Filtered.assemb.SE.fasta",
                         "-o"+ReassembScaffoldPath+"/"+SampleKey+".x"+Trim5+".Filtered.assemb.scaffold"]
        logger.debug(ParamList)
        Process = Popen(ParamList,shell=False)
        sts = Process.wait()
        logger.info("Scaffolding of reassembled '%s' completed" % (SampleKey))
    logger.info("All Scaffolding completed")

logger.info("All completed")

