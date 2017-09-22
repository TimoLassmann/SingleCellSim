
#ifndef outdir_header

#define outdir_header

#define OUTDIR_PARAM "param"
#define OUTDIR_CHECKPOINTS  "checks" 
#define OUTDIR_LOG "logs"

#define LOGFILE "log.txt"

extern int create_output_directories(char* root);
extern int set_log_file(char* root, char* name);

extern int create_dir(char* name, int force);

extern int check_if_output_directories_exists(char* root);

extern char** get_files_names(char* directory, char* suffix, char** list, int* num_files,int depth);

extern char* concat_out_dirs(char* root, char* sub); 

#endif
