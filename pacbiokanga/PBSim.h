#pragma once

#define PBSUCCEEDED 1
#define PBFAILED 0
#define PBTRUE 1
#define PBFALSE 0
#define BUF_SIZE 10240
#define REF_ID_LEN_MAX 128
#define REF_SEQ_NUM_MAX 9999
#define REF_SEQ_LEN_MAX 1000000000
#define REF_SEQ_LEN_MIN 100
#define FASTQ_NUM_MAX 100000000
#define FASTQ_LEN_MAX 1000000
#define FASTQ_LEN_MIN 100
#define EXISTS_LINE_FEED 1
#define DEPTH_MAX 1000
#define RATIO_MAX 1000
#define DATA_TYPE_CLR 0
#define DATA_TYPE_CCS 1
#define PROCESS_MODEL 1
#define PROCESS_SAMPLING 2
#define PROCESS_SAMPLING_REUSE 3
#define PROCESS_SAMPLING_STORE 4

#define PBS_PI       3.14159265358979323846

#pragma pack(1)
/////////////////////////////////////////
// Definitions of structures           //
/////////////////////////////////////////

// simulation parameters
struct sim_t {
  int data_type;
  int process;
  double depth;
  double accuracy_mean, accuracy_sd, accuracy_max, accuracy_min;
  long len_min, len_max; 
  double len_mean, len_sd; 
  long long len_quota;
  long sub_ratio, ins_ratio, del_ratio;
  double sub_rate, ins_rate, del_rate;
  int set_flg[20];
  long res_num;
  long long res_len_total; 
  double res_depth;
  double res_accuracy_mean, res_accuracy_sd;
  long res_len_min, res_len_max; 
  double res_len_mean, res_len_sd; 
  long res_sub_num, res_ins_num, res_del_num;
  double res_sub_rate, res_ins_rate, res_del_rate;
  char *prefix, *outfile_ref, *outfile_fq, *outfile_maf;
  char *model_qc_file;
  char *profile_id, *profile_fq, *profile_stats;
  unsigned int seed;
  char szref_file[200];
  char szFastqFileName[200];
};

// FASTQ
struct fastq_t {
  char *file;
  long num;
  long long len_total;
  long len_min, len_max;
  long num_filtered;
  long long len_total_filtered;
  long len_min_filtered, len_max_filtered;
  double len_mean_filtered, len_sd_filtered;
  double accuracy_mean_filtered, accuracy_sd_filtered;
};

// Reference
struct ref_t {
  char *file;
  char *seq;
  char id[REF_ID_LEN_MAX + 10];
  long len;
  long num_seq;
  long num;
};

// Mutation
struct mut_t {
  long sub_thre[100], ins_thre[100], del_thre;
  char *sub_nt_a, *sub_nt_t, *sub_nt_g, *sub_nt_c, *sub_nt_n, *ins_nt;
  char *qc, *new_qc, *tmp_qc, *seq, *new_seq, *maf_seq, *maf_ref_seq;
  long tmp_len_max;
  char seq_strand;
  long seq_left, seq_right;
};

// Quality code
struct qc_t {
  char character;
  double prob;
};

// Quality code of model
struct model_qc_t {
  int min;
  int max;
  double prob[100];
};
#pragma pack()




class CPBSim
{
	/////////////////////////////////////////
// Global variances                    //
/////////////////////////////////////////

FILE *fp_filtered, *fp_stats, *fp_ref, *fp_fq, *fp_maf;
struct sim_t sim;
struct fastq_t fastq;
struct ref_t ref;
struct mut_t mut;
struct qc_t qc[100];
struct model_qc_t model_qc[110]; // not too confident that pbsim code will not reference past 101 so allow additional

long freq_len[FASTQ_LEN_MAX+10]; // not too confident that pbsim code will not reference past FASTQ_LEN_MAX so allow additional
long freq_accuracy[100100];   // not too confident that pbsim code will not reference past 100000 so allow additional

/////////////////////////////////////////
// Prototypes of functions             //
/////////////////////////////////////////

int trim(char *line);
void init_sim_res(void);
int set_sim_param(void);
int get_ref_inf(void);
int get_ref_seq(void);
int get_fastq_inf(void);
int set_model_qc(void);
int set_mut(void);
int simulate_by_sampling(void);
int simulate_by_model(void);
int mutate(void);
int count_digit(long num);
void print_sim_param(void);
void print_fastq_stats(void);
void print_simulation_stats(void);
void print_helpvoid();
void revcomp(char* str);
long get_time();
public:
	CPBSim();
	~CPBSim();
	int Process(struct sim_t *psim);
};

