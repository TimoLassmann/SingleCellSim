
#define OPT_SAMPLE 1
#define OPT_TARGET_GENES 2
#define OPT_BOUNDARIES 3
#define OPT_DATFILE 4
#define OPT_NUMCELLS 5
#define OPT_DEPTH 6
#define OPT_SDDEPTH 7
#define OPT_CV 8
#define OPT_LOSS 9
#define OPT_NUM_THREAD 10
#define OPT_DROPOUT 11
#define OPT_FITONLY 12

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <ctype.h>
#include "tldevel.h"

#include "user_interface.h"

struct parameters* init_param(void);
int byg_count(char* pattern, char* text);
int print_model_help(int argc, char * argv[]);
int print_sim_help(int argc, char * argv[]);

int print_global_help(int argc, char * argv[])
{
        const char usage[] = " <command>";
        fprintf(stdout,"\nUsage: %s %s\n\n",basename(argv[0]) ,usage);
        fprintf(stdout,"Commands:\n\n");

        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"model","Models expression data." ,"[NA]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"sim","Simulats datasets." ,"[NA]"  );
        fprintf(stdout,"\n");

        return OK;
}

int print_model_help(int argc, char * argv[])
{
        const char usage[] = " model [-options] -i <input table> -o <outdir> ";
        fprintf(stdout,"\nUsage: %s %s\n\n",basename(argv[0]) ,usage);
        fprintf(stdout,"Options:\n\n");

        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--in","Input expression table" ,"[NA]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--out","Output directory." ,"[NA]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--nthread","Number of threads." ,"[8]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--fitonly","Only fit model." ,"[NA]"  );
        fprintf(stdout,"\n");
        return OK;
}


struct parameters* get_model_param(int argc, char * argv[])
{
        struct parameters* param = NULL;
        int c, help;

        RUN(print_program_header(argv,"models transcriptome complexity."));


        if(argc == 2){
                RUN(print_model_help(argc,argv));
                exit(EXIT_SUCCESS);
        }


        help = 0;
        c = 0;

        RUNP(param= init_param());

        while (1){
                static struct option long_options[] ={
                        {"sample",required_argument,0,OPT_SAMPLE},
                        {"in",required_argument,0, 'i'},
                        {"out",required_argument,0, 'o'},
                        {"dropout",required_argument,0,OPT_DROPOUT},
                        {"nthread",required_argument,0,OPT_NUM_THREAD},
                        {"ncell",required_argument,0,OPT_NUMCELLS},
                        {"depth",required_argument,0,OPT_DEPTH},
                        {"CV",required_argument,0,OPT_CV},
                        {"fitonly",0,0,OPT_FITONLY},
                        {"help",0,0,'h'},
                        {0, 0, 0, 0}
                };

                int option_index = 0;
                c = getopt_long_only (argc, argv,"o:i:h",long_options, &option_index);

                if (c == -1){
                        break;
                }

                switch(c) {
                case 0:
                        break;
                case OPT_FITONLY:
                        param->fitonly = 1;
                        break;
                case OPT_DROPOUT:
                        param->drop_Michaelisconstant = atof(optarg);
                        break;
                case OPT_CV:
                        param->CV = atof(optarg);
                        break;
                case OPT_NUMCELLS:
                        param->num_cells = atoi(optarg);
                        break;
                case OPT_DEPTH:
                        param->simulated_read_depth = atoi(optarg);
                        break;
                case OPT_NUM_THREAD:
                        param->num_threads = atoi(optarg);
                        break;
                case 'h':
                        help = 1;
                        break;
                case 'i':
                        param->infile = optarg;
                        break;
                case 'o':
                        param->outdir = optarg;
                        break;
                case '?':
                        exit(EXIT_FAILURE);
                        break;
                default:
                        WARNING_MSG("unrecognized option.");
                        break;
                }
        }

        if(help){
                RUN(print_model_help(argc,argv));
                exit(EXIT_SUCCESS);
        }
        ASSERT(param->infile != NULL,"No input file.");
        ASSERT(param->outdir != NULL,"No output directory.");


        RUN(log_command_line(argc, argv));

        LOG_MSG("Read param.");
        return param;
ERROR:
        LOG_MSG("Something went wrong. Try using the -h option.");
        free_param(param);
        return NULL;
}

int print_sim_help(int argc, char * argv[])
{
        const char usage[] = " sim [-options] -i <model root> -o <outdir> ";
        fprintf(stdout,"\nUsage: %s %s\n\n",basename(argv[0]) ,usage);
        fprintf(stdout,"Options:\n\n");

        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--in","Path to model directory." ,"[NA]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--ncell","Number of cells." ,"[96]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--depth","Read depth." ,"[80,000]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--depthSD","Read depth SD." ,"[10,000]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--nthread","Number of threads." ,"[8]"  );
        fprintf(stdout,"\n");
        return OK;
}


struct parameters* get_sim_param(int argc, char * argv[])
{
        struct parameters* param = NULL;
        int i,c,f,g, count, help;
        char *tmp;


        RUN(print_program_header(argv,"models trascriptome complexity."));

        if(argc == 2){
                RUN(print_sim_help(argc,argv));
                exit(EXIT_SUCCESS);
        }

        help = 0;
        c = 0;

        RUNP(param= init_param());

        while (1){
                static struct option long_options[] ={
                        {"gene",required_argument,0,OPT_TARGET_GENES},
                        {"sample",required_argument,0,OPT_SAMPLE},
                        {"in",required_argument,0, 'i'},
                        {"out",required_argument,0, 'o'},
                        {"dropout",required_argument,0,OPT_DROPOUT},
                        {"nthread",required_argument,0,OPT_NUM_THREAD},
                        {"ncell",required_argument,0,OPT_NUMCELLS},
                        {"depth",required_argument,0,OPT_DEPTH},
                        {"depthSD",required_argument,0,OPT_SDDEPTH},
                        {"CV",required_argument,0,OPT_CV},
                        {"help",0,0,'h'},
                        {0, 0, 0, 0}
                };

                int option_index = 0;
                c = getopt_long_only (argc, argv,"o:i:h",long_options, &option_index);

                if (c == -1){
                        break;
                }

                switch(c) {
                case 0:
                        break;
                case OPT_TARGET_GENES:
                        tmp = optarg;
                        count = byg_count(",", tmp);
                        //fprintf(stderr,"%d\n",count);
                        param->num_gene_targets = count+1;


                        MMALLOC(param->gene_names,sizeof(char*) * (count+1));

                        for(i = 0; i < count+1;i++){
                                param->gene_names[i] = 0;
                                MMALLOC(param->gene_names[i],sizeof(char)* 50);
                        }
                        f = 0;
                        g = 0;
                        for(i = 0;i < strlen(tmp);i++){
                                if(tmp[i] != ','){
                                        if(g != 49){
                                                param->gene_names[f][g] = tmp[i];
                                        }
                                        g++;

                                }else{
                                        param->gene_names[f][g] = 0;
                                        f++;
                                        g = 0;
                                }
                        }
                        param->gene_names[f][g] = 0;

                        break;
                case OPT_SAMPLE:
                        tmp = optarg;
                        count = byg_count(",", tmp);
                        //fprintf(stderr,"%d\n",count);
                        param->num_samples = count+1;


                        MMALLOC(param->sample_names,sizeof(char*) * (count+1));

                        for(i = 0; i < count+1;i++){
                                param->sample_names[i] = 0;
                                MMALLOC(param->sample_names[i],sizeof(char)* 50);
                        }
                        f = 0;
                        g = 0;
                        for(i = 0;i < strlen(tmp);i++){
                                if(tmp[i] != ','){
                                        if( g != 50-1){
                                                param->sample_names[f][g] = tmp[i];
                                        }
                                        g++;
                                }else{
                                        param->sample_names[f][g] = 0;
                                        f++;
                                        g = 0;
                                }
                        }
                        param->sample_names[f][g] = 0;

                        break;

                case OPT_DROPOUT:
                        param->drop_Michaelisconstant = atof(optarg);
                        break;
                case OPT_CV:
                        param->CV = atof(optarg);
                        break;
                case OPT_NUMCELLS:
                        param->num_cells = atoi(optarg);
                        break;
                case OPT_DEPTH:
                        param->simulated_read_depth = atoi(optarg);
                        break;
                case OPT_SDDEPTH:
                        param->simulated_read_depth_SD = atof(optarg);
                        break;
                case OPT_NUM_THREAD:
                        param->num_threads = atoi(optarg);
                        break;
                case 'h':
                        help = 1;
                        break;
                case 'i':
                        param->infile = optarg;
                        break;
                case 'o':
                        param->outdir = optarg;
                        break;
                case '?':
                        exit(EXIT_FAILURE);
                        break;
                default:
                        WARNING_MSG("unrecognized option.");
                        break;
                }
        }

        if(help){
                RUN(print_sim_help(argc,argv));
                exit(EXIT_SUCCESS);
        }
        ASSERT(param->infile != NULL,"No input file.");
        ASSERT(param->outdir != NULL,"No output directory.");
        RUN(log_command_line(argc, argv));

        return param;
ERROR:
        LOG_MSG("Something went wrong. Try using the -h option.");
        free_param(param);
        return NULL;
}


struct parameters* init_param(void)
{
        struct parameters* param = NULL;
        MMALLOC(param, sizeof(struct parameters));
        param->sample_names = NULL;
        param->gene_names = NULL;
        param->infile = NULL;
        param->outdir = NULL;
        param->num_threads = 8;
        param->num_genes = 0;
        param->num_gene_targets = 0;
        param->simulated_read_depth = 80000;
        param->simulated_read_depth_SD = 10000.0;
        param->CV = 0.75;
        param->drop_Michaelisconstant = 10.0;
        param->num_samples = 0;
        param->num_cells = 96;
        param->fitonly = 0;
        return param;
ERROR:
        return NULL;
}


void free_param(struct parameters* param)
{
        int i;
        if(param){
                if(param->num_samples){
                        for(i = 0; i < param->num_samples;i++){
                                MFREE(param->sample_names[i]);
                        }
                        MFREE(param->sample_names);
                }
                if(param->num_gene_targets){
                        for(i = 0; i < param->num_gene_targets;i++){
                                MFREE(param->gene_names[i]);
                        }
                        MFREE(param->gene_names);

                }
                MFREE(param);
        }
}



int byg_count(char* pattern,char*text)
{
        int Tc;
        int i  = 0;
        int s = 0;
        int T[256];
        for (i = 0;i < 256;i++){
                T[i] = 0;
        }

        int m = (int) strlen(pattern);
        int n = (int) strlen(text);
        int count = 0;

        if(m > n){
                return -1;
        }
        int mb = (1 << (m-1));

        for (i= 0;i < m;i++){
                T[(int)toupper(pattern[i])] |= (1 << i);
        }

        for (i = 0;i < n;i++){
                s <<= 1;
                s |= 1;
                Tc = T[(int)toupper(text[i])];
                s &= Tc;
                if(s & mb){
                        count++;
                }
        }
        return count;
}



