#include "scs.h"
#include "user_interface.h"
#include "outdir.h"
#include "gigp.h"

struct unique_transcript_count{
        double* x;
        int max;
        int alloc_size;
};

struct gene_ex{
        char* name;
        double expr;
};

struct shared_data{
        struct parameters* param;
        struct unique_transcript_count* utc;
        struct double_matrix* gigp_param_out_table;
        struct double_matrix* rel_abundances;
        struct gigp_param** working_gigp_param;
        struct gigp_param* gigp_param_best;
        struct gigp_param* gigp_param_initial_guess;
        struct thr_pool* pool;
        uint8_t* available_work_space;
        pthread_mutex_t avail_mtx;
       	int num_threads;
};

struct thread_data{
        struct shared_data* bsd;
        struct gigp_param* gigp_param;
        double* vector;
        double try_gamma;
        double try_b;
        double try_c;
        int target;
        int thread_id;
        int num_threads;
};

int run_model_scs(struct parameters* param);

int make_sim_expression_table(struct shared_data* bsd, struct double_matrix* m,int target);
int cmp_gene_ex_struct(const void *a, const void *b);

int run_sim_scs(struct parameters* param);

int filter_matrix_tpm(struct double_matrix* m, double tpm_threshold);

int init_gigp_param_out_table(struct shared_data* bsd, struct double_matrix* input);

int get_unique_transcript_vector(struct shared_data* bsd,struct double_matrix* m, int column);

int fit_models_to_table(struct shared_data* bsd, struct double_matrix* m);

int fit_model(struct shared_data* bsd);

int fit_gigp(struct gigp_param* gigp_param,double* fx_unique_transcript_count, int len);

int set_initial_guess(struct shared_data* bsd);

void* do_fit_model(void *threadarg);

int enter_gigp_param_into_table(struct double_matrix* table, struct gigp_param* best,int col);
int read_gigp_param_from_table(struct double_matrix* table, struct gigp_param* best,int col);

int write_gigp_param_table_to_file(struct shared_data* bsd);
int read_gigp_param_table_from_file(struct shared_data* bsd);

int fit_curves_make_growth_curve(struct shared_data* bsd, struct double_matrix* m);

int fill_rel_abundance_matrix(struct shared_data* bsd, struct double_matrix* m);
void* run_get_rel_abundances (void *threadarg);

int fingerprints(struct shared_data* bsd, struct double_matrix* m);

int print_fit(struct shared_data* bsd);

int calculate_rel_frequencies(struct shared_data* bsd);

int sample_GIGP_same_depth_as_input(struct shared_data* bsd);

int simulate_expression_table(struct shared_data* bsd);

int write_fitted_param_to_file(struct gigp_param* gigp_param, char* filename);

int read_fitted_param_from_file(struct gigp_param* gigp_param, char* filename);

struct shared_data* init_shared_data(struct parameters* param);

void free_shared_data(struct shared_data* bsd);

static int compare_double (const void * a, const void * b);

int fill_fitted_curve(double* fit,int num,struct gigp_param* param );

double* pick_abundances(struct gigp_param* gigp_param, double* random, double* abundances,int outer_low,int outer_high,double inner_low,double inner_high);



int main (int argc, char * argv[])
{
        struct parameters* param = NULL;


        print_program_header(argv, "Models transcriptomes.");

        if(argc == 1){
                /* print help */
                RUN(print_global_help(argc,argv));
                return OK;
        }

        if(strncmp(argv[1],"model",5) == 0){
                /* get parameters */
                RUNP(param = get_model_param(argc,argv));
                /* create output structure ...  */
                RUN(run_model_scs(param));
        }else if(strncmp(argv[1],"sim",3) == 0){
                RUNP(param = get_sim_param(argc,argv));
                RUN(run_sim_scs(param));
        }else{
                ERROR_MSG("Option %s not recognized.",argv[1]);
        }
        free_param(param);
        return EXIT_SUCCESS;
ERROR:
        free_param(param);
        return EXIT_FAILURE;
}

int run_model_scs(struct parameters* param)
{
        struct shared_data* bsd = NULL;
        struct double_matrix* m = NULL;

        char buffer[BUFFER_LEN];

        RUN(create_output_directories(param->outdir));

        RUN(set_log_file(param->outdir,"scs_net"));

        /* allocate shared data  */
        RUNP(bsd = init_shared_data(param));

        /* create checkpoint... */
        snprintf(buffer, BUFFER_LEN, "%s/%s/",param->outdir,OUTDIR_CHECKPOINTS);
        DECLARE_CHK(MAIN_CHECK, buffer);

        /* Read in counts..  */

        RUNP(m =  read_dm(param->infile,1,1));

        /* filter outlowly expressed stuff  */
        RUN(filter_matrix_tpm(m,1.0));
        /* init gigparam->out. */
        RUN(init_gigp_param_out_table(bsd,m));
        //print_double_matrix(bsd->gigp_param_out_table,stdout,1,1);

        /* do modelling for every sample in input table..  */
        LOG_MSG("Fit models.");
        snprintf(buffer, BUFFER_LEN, "infile: %s datadim:%d %d.",param->infile, m->ncol,m->nrow);
        RUN_CHECKPOINT(MAIN_CHECK,fit_models_to_table(bsd,m),buffer);

        LOG_MSG("Write preliminary model parameters.");
        snprintf(buffer, BUFFER_LEN, "Write preliminary parameters.");
        RUN_CHECKPOINT(MAIN_CHECK,write_gigp_param_table_to_file(bsd),buffer);

        LOG_MSG("Make growth curves.");
        /* fit curves, print and calculate S (estmated number of transcripts in samples.  */
        snprintf(buffer, BUFFER_LEN, "datadim:%d %d.",m->ncol,m->nrow);
        RUN_CHECKPOINT(MAIN_CHECK,fit_curves_make_growth_curve(bsd, m), buffer);

        /* critical! writes estimated total number of transcripts per cell / sample  */
        snprintf(buffer, BUFFER_LEN, "Write parameters.");
        RUN_CHECKPOINT(MAIN_CHECK,write_gigp_param_table_to_file(bsd),buffer);

        /* calculate relative abundances...  */
        //RUN(print_double_matrix(bsd->gigp_param_out_table,stdout ,1,1));
        snprintf(buffer, BUFFER_LEN, "Calculate relative abundances.");
        RUN_CHECKPOINT(MAIN_CHECK,fill_rel_abundance_matrix(bsd,m),buffer);

        /* calculare X genes in low 10%, low 20% etc... */
        snprintf(buffer, BUFFER_LEN, "Calculate fingerprint based on top genes.");
        RUN_CHECKPOINT(MAIN_CHECK,fingerprints(bsd,m), buffer);

        free_double_matrix(m);
        free_shared_data(bsd);
        DESTROY_CHK(MAIN_CHECK);
        return OK;
        exit(0);

        /* calcuate relative frequencies from distribution and simulate. */
        snprintf(buffer, BUFFER_LEN, "%lf %lf %lf",round(bsd->gigp_param_best->gamma*100000)/100000.0f, round(bsd->gigp_param_best->b*100000)/100000.0f,round(bsd->gigp_param_best->c*100000)/100000.0f);
        RUN_CHECKPOINT(MAIN_CHECK,calculate_rel_frequencies(bsd),buffer);

        snprintf(buffer, BUFFER_LEN, "%lf %lf %lf",round(bsd->gigp_param_best->gamma*100000)/100000.0f, round(bsd->gigp_param_best->b*100000)/100000.0f,round(bsd->gigp_param_best->c*100000)/100000.0f);
        RUN_CHECKPOINT(MAIN_CHECK,sample_GIGP_same_depth_as_input(bsd),buffer);

        snprintf(buffer, BUFFER_LEN, "gamma:%lf b:%lf c:%lf simdepth:%d dropout:%f",round(bsd->gigp_param_best->gamma*100000)/100000.0f, round(bsd->gigp_param_best->b*100000)/100000.0f,round(bsd->gigp_param_best->c*100000)/100000.0f, bsd->param->simulated_read_depth, bsd->param->drop_Michaelisconstant);

        RUN_CHECKPOINT(MAIN_CHECK,simulate_expression_table(bsd),buffer);


        free_double_matrix(m);
        free_shared_data(bsd);
        return OK;
ERROR:
        free_shared_data(bsd);
        free_double_matrix(m);
        return FAIL;
}

int run_sim_scs(struct parameters* param)
{
        char buffer[BUFFER_LEN];

        struct shared_data* bsd = NULL;
        struct double_matrix* m = NULL;

        int i,j;
        int hits = 0;

        ASSERT(param != NULL,"No parameters.");

        RUNP(bsd = init_shared_data(param));
        /* read parameters - why? dunno  */
        RUN(read_gigp_param_table_from_file(bsd));
        /* read rel_abundance  */
        if(!bsd->rel_abundances){
                LOG_MSG("Reading relative abundances.");
                snprintf(buffer, BUFFER_LEN, "%s/%s/Relabundances.csv",bsd->param->outdir,OUTDIR_MODEL);
                LOG_MSG("Read from file: %s",buffer);
                RUNP(bsd->rel_abundances = read_double_matrix(buffer,1,1));
                RUNP(bsd->rel_abundances = transpose_double_matrix(bsd->rel_abundances));
        }

        /* Read in counts..  */
        RUNP(m = read_double_matrix(param->infile,1,1));

        /* make sure expression table is the one used in the modelling */
        ASSERT(m->ncol == bsd->rel_abundances->nrow,"Number of samples in expression table differs from samples in model file!");

        for(i = 0;i < m->ncol;i++){
                ASSERT(strcmp(m->col_names[i],bsd->rel_abundances->row_names[i]) == 0,"Different samples: %s %s!",m->col_names[i],bsd->rel_abundances->row_names[i]);
        }

        //print_double_matrix(bsd->rel_abundances,stdout,1,1);
        /* loop through sample array in param and for each sample simulate...  */
        for(i = 0; i < bsd->param->num_samples;i++){
                hits = 0;
                for(j = 0; j < bsd->rel_abundances->nrow;j++){
                        if(strstr(bsd->rel_abundances->row_names[j],bsd->param->sample_names[i])!= NULL){
                                hits++;
                        }
                }
                switch (hits) {
                case 0: {
                        WARNING_MSG("Sample: \"%s\" not found.", bsd->param->sample_names[i]);
                        break;
                }
                case 1: {
                        for(j = 0; j < bsd->rel_abundances->nrow;j++){
                                if(strstr(bsd->rel_abundances->row_names[j],bsd->param->sample_names[i])!= NULL){

                                        RUN(make_sim_expression_table(bsd,m,j));
                                }
                        }
                        //Do stuff....
                        break;
                }
                default:
                        WARNING_MSG("Sample: \"%s\" mas multiple hits.", bsd->param->sample_names[i]);
                        break;
                }

        }
        //char *strstr(const char *s1, const char *s2);


        /* highlight genes...j */

        free_shared_data(bsd);
        return OK;
ERROR:
        free_shared_data(bsd);
        free_double_matrix(m);
        return FAIL;
}

int filter_matrix_tpm(struct double_matrix* m, double tpm_threshold)
{
        int i,j;
        double* sums = NULL;
        int* keep = NULL;
        int* dropped = NULL;

        MMALLOC(sums, sizeof(double)* m->ncol);
        MMALLOC(keep, sizeof(double)* m->ncol);
        MMALLOC(dropped, sizeof(double)* m->ncol);
        for(i = 0; i <  m->ncol;i++){
                sums[i] = 0.0;
                keep[i] = 0;
                dropped[i] = 0;
        }

        for(i = 0; i <  m->nrow;i++){
                for(j = 0; j < m->ncol;j++){
                        sums[j]+= m->matrix[i][j];
                }
        }

        for(i = 0; i <  m->nrow;i++){
                for(j = 0; j < m->ncol;j++){
                        if(m->matrix[i][j] / sums[j] * 1000000.0 < tpm_threshold){
                                m->matrix[i][j] = 0.0;
                                dropped[j]++;
                        }else{
                                keep[j]++;
                        }
                }
        }

        fprintf(stdout,"Sample\tkeep\tdropped\n");
        for(j = 0; j < m->ncol;j++){
                fprintf(stdout,"%s\t%d\t%d\n", m->col_names[j], keep[j],dropped[j]);
        }


        MFREE(sums);
        MFREE(keep);
        MFREE(dropped);



        return OK;
ERROR:
        MFREE(sums);
        MFREE(keep);
        MFREE(dropped);
        return FAIL;
}

int make_sim_expression_table(struct shared_data* bsd, struct double_matrix* m,int target)
{
        char buffer[BUFFER_LEN];
        struct double_matrix* sim_table = NULL;
        struct gene_ex** sort_genes = NULL;
        double* rel_abundances = NULL;
        double* cell_abundance = NULL;
        FILE* out_ptr = NULL;

        double CV;
        double libdepth;
        double drop_k;


        int i,j,c;
        int num_genes;
        int num_cells;
        int simulated_read_depth;
        ASSERT(bsd != NULL,"No shared data.");
        ASSERT(m != NULL,"No expression data.");

        rel_abundances = bsd->rel_abundances->matrix[target];

        num_genes = bsd->rel_abundances->ncol;
        for(i =0; i < bsd->rel_abundances->ncol;i++){
                //        fprintf(stdout,"%f\n", rel_abundances[i]);
                if(rel_abundances[i] == 0.0){
                        num_genes = i;
                        break;
                }
        }
        LOG_MSG("%d genes found in sample %d... ", num_genes,target);

        num_cells = bsd->param->num_cells;
        simulated_read_depth = bsd->param->simulated_read_depth;
        CV = bsd->param->CV;
        drop_k = bsd->param->drop_Michaelisconstant;




        /* red in relative abundancs */
        MMALLOC(cell_abundance, sizeof(double) * num_genes);

        /* set up random number generator  */
        gsl_rng_env_setup();
        const gsl_rng_type *T = gsl_rng_default;
        gsl_rng *rfff = gsl_rng_alloc (T);
        LOG_MSG("Allocating expression table: %d %d", bsd->param->num_cells, num_genes);
        /* alloc outout table and aux. vector  */
        RUNP(sim_table = alloc_double_matrix(num_cells,num_genes, BUFFER_LEN));
        /* fill fake names */
        for(i = 0; i < num_cells ;i++){
                snprintf(sim_table->col_names[i],BUFFER_LEN,"Cell%d",i+1);
        }
        for(j = 0; j < num_genes;j++){
                snprintf(sim_table->row_names[j],BUFFER_LEN,"Blankgene%d",j+1);
        }

        /* sort genes in target sample by expression */
        MMALLOC(sort_genes,sizeof(struct gene_ex*) * m->nrow);
        for(i = 0; i < m->nrow;i++){
                sort_genes[i] = NULL;
                MMALLOC(sort_genes[i],sizeof(struct gene_ex));
                sort_genes[i]->name = m->row_names[i];
                sort_genes[i]->expr = m->matrix[i][target];
        }
        qsort(sort_genes, m->nrow, sizeof(struct gene_ex*), cmp_gene_ex_struct);

        j = num_genes-1;
        for(i = 0; i < m->nrow;i++){
                if(j == 0){
                        break;
                }
                for(c = 0; c < bsd->param->num_gene_targets;c++){
                        if(strstr(sort_genes[i]->name, bsd->param->gene_names[c])!= NULL){
                                snprintf(sim_table->row_names[j],BUFFER_LEN,"%s",sort_genes[i]->name);
                        }
                }
                j--;
        }



        /* fill expression table  */

        //Calculate tpm..
        libdepth = 0.0;
        for(j = 0; j < num_genes;j++){
                libdepth += rel_abundances[j];
        }


        for(i = 0; i < num_cells ;i++){
                for(j = 0; j < num_genes;j++){
                        cell_abundance[j] = rel_abundances[j] + gsl_ran_gaussian (rfff,  rel_abundances[j] * CV);
                        if(cell_abundance[j] < 0){
                                cell_abundance[j] = 0.0;
                        }
                        double average = rel_abundances[j] / libdepth * 1000000.0;
                        double p_drop = 1.0 - (average) / (average + drop_k);
                        if(random_float_zero_to_x(1.0) <= p_drop){
                                cell_abundance[j] = 0.0;
                        }

                }
                gsl_ran_discrete_t *sample_abundances = gsl_ran_discrete_preproc (num_genes, cell_abundance);
                for(j = 0;j < simulated_read_depth;j++){
                        sim_table->matrix[(int)gsl_ran_discrete(rfff, sample_abundances)][i] += 1;
                }
                gsl_ran_discrete_free(sample_abundances);
        }

        snprintf(buffer, BUFFER_LEN, "%s/%s/extable_%s_%d_%d_CV%f.csv",bsd->param->outdir,OUTDIR_SIM,bsd->rel_abundances->row_names[target],num_cells, simulated_read_depth,CV);
        LOG_MSG("Write to %s\n",buffer);
        RUNP(out_ptr = fopen(buffer, "w"));
        RUN(print_double_matrix(sim_table,out_ptr, 1,1));
        fclose(out_ptr);

        gsl_rng_free(rfff);
        MFREE(cell_abundance);

        free_double_matrix(sim_table);
        return OK;
ERROR:
        return FAIL;
}

int fingerprints(struct shared_data* bsd, struct double_matrix* m)
{
        char buffer[BUFFER_LEN];
        FILE* file_ptr = NULL;
        struct double_matrix* f = NULL;
        double sum;
        int i,j,min_S;
        int row;
        int step = 10;
        int max_genes = 2000;

        int* top_intervals = NULL;

        //int top_intervals[] = {0,5,10,50,100,250,500,750,1000,1500,2000};
        ASSERT(bsd != NULL, "No shared data.");
        ASSERT(m   != NULL, "No matrix.");

        MMALLOC(top_intervals,sizeof(int) * (max_genes / step + 1));

        top_intervals[0] = 0;
        for(i = 1; i <= max_genes /step;i++){
                top_intervals[i] = top_intervals[i-1] + step;
        }

        RUN(read_gigp_param_table_from_file(bsd));

        if(!bsd->rel_abundances){
                LOG_MSG("Reading relative abundances.");
                snprintf(buffer, BUFFER_LEN, "%s/%s/Relabundances.csv",bsd->param->outdir,OUTDIR_MODEL);
                LOG_MSG("Read from file: %s",buffer);
                RUNP(bsd->rel_abundances = read_double_matrix(buffer,1,1));
        }


        /* convert distrete probabilities of 1 .. S genes into a 0 (low
         * expressed) - 1 ( high expressed vector) with X intervals */
        /* parsing top percentages .  */
        double resolution = 100;

        RUNP(f = alloc_double_matrix(m->ncol, resolution,BUFFER_LEN));
        for(i = 0; i < m->ncol;i++){
                snprintf(f->col_names[i],BUFFER_LEN, "%s",m->col_names[i]);
        }

        for(i = 0; i < resolution;i++){
                snprintf(f->row_names[i],BUFFER_LEN, "Prob%0.2f_%0.2f", (double) i / resolution,(double) (i+1) / resolution);
        }

        for(i = 0; i < m->ncol;i++){
                int S = bsd->gigp_param_out_table->matrix[GIGP_S_ROW][i];
                for(j = 0; j <= S ;j++){
                        f->matrix[ (int)floor ((double) (j) / (double) (S+0.0001) * resolution) ][i] += bsd->rel_abundances->matrix[j][i];
                }
        }
        snprintf(buffer, BUFFER_LEN, "%s/%s/fingerprint_dist.csv",bsd->param->outdir,OUTDIR_MODEL);
        LOG_MSG("Writing to file: %s",buffer);
        RUNP(file_ptr = fopen(buffer,"w"));

        print_double_matrix(f,file_ptr,1,1);

        fclose(file_ptr);

        free_double_matrix(f);
        f = NULL;

        RUNP(f = alloc_double_matrix(m->ncol,(max_genes / step + 1), BUFFER_LEN));

        for(i = 0; i < m->ncol;i++){
                snprintf(f->col_names[i],BUFFER_LEN, "%s",m->col_names[i]);
        }

        for(i = 1; i < (max_genes / step + 1);i++){
                snprintf(f->row_names[i-1], BUFFER_LEN, "TOP%d",  top_intervals[i]);
        }

        min_S = INT_MAX;
        for(i = 0; i < bsd->gigp_param_out_table->ncol;i++){
                if(bsd->gigp_param_out_table->matrix[GIGP_S_ROW][i] < min_S){
                        min_S = bsd->gigp_param_out_table->matrix[GIGP_S_ROW][i];
                }
        }

        for(i = 0; i < m->ncol;i++){
                sum = 0.0;
                for(row = 1;row < (max_genes / step + 1);row++){
                        if(top_intervals[row] > min_S){
                                f->matrix[row-1][i] = 1.0;
                        }else{
                                for(j = bsd->gigp_param_out_table->matrix[GIGP_S_ROW][i] - top_intervals[row-1] ; j > bsd->gigp_param_out_table->matrix[GIGP_S_ROW][i] - top_intervals[row];j--){
                                        sum += bsd->rel_abundances->matrix[j][i];
                                }
                                f->matrix[row-1][i] = sum;
                        }
                }
        }

        snprintf(buffer, BUFFER_LEN, "%s/%s/fingerprint.csv",bsd->param->outdir,OUTDIR_MODEL);
        LOG_MSG("Writing to file: %s",buffer);
        RUNP(file_ptr = fopen(buffer,"w"));

        print_double_matrix(f,file_ptr,1,1);

        fclose(file_ptr);

        free_double_matrix(f);
        return OK;
ERROR:
        if(f){
                free_double_matrix(f);
        }
        return FAIL;
}

int fit_models_to_table(struct shared_data* bsd, struct double_matrix* m)
{
        int i;
        ASSERT(bsd != NULL, "No shared data.");
        ASSERT(m   != NULL, "No matrix.");
        DECLARE_TIMER(t1);

        for(i = 0 ;i < m->ncol;i++){
                START_TIMER(t1);
                LOG_MSG("Working on sample:%s.",m->col_names[i]);
                /* Turn counts into frequency vector. */
                RUN(get_unique_transcript_vector(bsd,m,i));
                if(bsd->utc->max > 1000000){
                        WARNING_MSG("Super high abundant transcript detected: %d.", bsd->utc->max);
                        WARNING_MSG("All steps may take considerable amount of time from this point on.");
                }
                /* do modelling...  */
                RUN(fit_model(bsd));
                /* enter_estimated parameters in table..   */
                RUN(enter_gigp_param_into_table(bsd->gigp_param_out_table,bsd->gigp_param_best,i));

                /* reset best param.... */
                RUN(clear_gigp_param(bsd->gigp_param_best));
                STOP_TIMER(t1);

                /* Done  */
                LOG_MSG("Done in %f seconds.",GET_TIMING(t1));
        }
        return OK;
ERROR:
        return FAIL;
}

int fill_rel_abundance_matrix(struct shared_data* bsd, struct double_matrix* m)
{
        char buffer[BUFFER_LEN];
        FILE* out_ptr = NULL;
        struct thread_data** td = NULL;
        double sum;
        int status;
        int max_S;
        int i,j;
        ASSERT(bsd != NULL, "No shared data!");
        ASSERT(bsd->rel_abundances == NULL,"Rel abundance matrix not NULL.");

        RUN(read_gigp_param_table_from_file(bsd));

        max_S = 0;
        for(i = 0; i < bsd->gigp_param_out_table->ncol;i++){
                if(bsd->gigp_param_out_table->matrix[GIGP_S_ROW][i] > max_S){
                        max_S = bsd->gigp_param_out_table->matrix[GIGP_S_ROW][i];
                }
                fprintf(stdout,"%f %d\n", bsd->gigp_param_out_table->matrix[GIGP_S_ROW][i], max_S);
        }
        /* confusing! the matrix is flipped - samples in rows! */
        RUNP(bsd->rel_abundances = alloc_double_matrix(  max_S +1,m->ncol, BUFFER_LEN));

        for(i = 0; i < m->ncol;i++){
                snprintf(bsd->rel_abundances->row_names[i], BUFFER_LEN, "%s",m->col_names[i]);
        }

        for(i = 0; i < (max_S +1);i++){
                snprintf(bsd->rel_abundances->col_names[i], BUFFER_LEN, "%d",i);
        }
        LOG_MSG("Start creating relative abundances.");
        MMALLOC(td, sizeof(struct thread_data*) * m->ncol);

        for(i = 0; i < m->ncol;i++){
                td[i] = NULL;
                MMALLOC(td[i], sizeof(struct thread_data));
                td[i]->thread_id = i;
                td[i]->target = i;
                td[i]->bsd = bsd;
                td[i]->vector = NULL;
                MMALLOC(td[i]->vector,sizeof(double) * (max_S+1));
                td[i]->num_threads = bsd->param->num_threads;
                RUNP(td[i]->gigp_param = init_gigp_param());
                RUN(read_gigp_param_from_table(bsd->gigp_param_out_table,td[i]->gigp_param,i));
                td[i]->gigp_param->S = 100;

                if((status = thr_pool_queue(bsd->pool,run_get_rel_abundances,td[i])) == -1) fprintf(stderr,"Adding job to queue failed.");
        }
        thr_pool_wait(bsd->pool);
        RUNP(bsd->rel_abundances = transpose_double_matrix(bsd->rel_abundances));

        //make sure dist sums to one...
        for(i = 0; i < m->ncol;i++){
                sum = 0.0;
                for(j = 0;j <=  bsd->gigp_param_out_table->matrix[GIGP_S_ROW][i];j++){
                        sum += bsd->rel_abundances->matrix[j][i];
                }

                fprintf(stdout,"%s %f S:%f \n", m->col_names[i], sum,bsd->gigp_param_out_table->matrix[GIGP_S_ROW][i] );
                for(j = 0;j < bsd->gigp_param_out_table->matrix[GIGP_S_ROW][i];j++){
                        //bsd->rel_abundances->matrix[j][i] /= sum;
                }
        }

        snprintf(buffer, BUFFER_LEN, "%s/%s/Relabundances.csv",bsd->param->outdir,OUTDIR_MODEL);
        LOG_MSG("Write to file: %s",buffer);
        RUNP(out_ptr = fopen(buffer, "w" ));

        RUN(print_double_matrix(bsd->rel_abundances,out_ptr,1,1));
        fclose(out_ptr);

        for(i = 0; i < m->ncol;i++){
                for(j = 1;j < bsd->gigp_param_out_table->matrix[GIGP_S_ROW][i];j++){
                        bsd->rel_abundances->matrix[j][i] += bsd->rel_abundances->matrix[j-1][i];
                }
        }

        snprintf(buffer, BUFFER_LEN, "%s/%s/Relabundances_sum.csv",bsd->param->outdir,OUTDIR_MODEL);
        LOG_MSG("Write to file: %s",buffer);
        RUNP(out_ptr = fopen(buffer, "w" ));

        RUN(print_double_matrix(bsd->rel_abundances,out_ptr,1,1));
        fclose(out_ptr);

        free_double_matrix(bsd->rel_abundances);
        bsd->rel_abundances = NULL;

        for(i = 0; i < m->ncol;i++){
                MFREE(td[i]->vector);
                MFREE(td[i]);
        }
        MFREE(td);
        return OK;
ERROR:
        return FAIL;
}

void* run_get_rel_abundances (void *threadarg)
{
        struct thread_data *data = NULL;
        struct gigp_param* gigp = NULL;
        double* observed = NULL;
        double* rel_abundance = NULL;
        int target;
        int i;
        data = (struct thread_data *) threadarg;

        /* copy param to working  param struct..  */

        ASSERT(data->bsd != NULL,"No shared data.");

        target = data->target;
        gigp = data->gigp_param;

        observed = data->vector;
        rel_abundance = data->bsd->rel_abundances->matrix[target];


        i = 0;
        rel_abundance[i] = 0.0;
        observed[i] =  1.0 / (gigp->S+1.0); // random_float_zero_to_x(1.0);

        for(i = 1; i <(int) (gigp->S+1);i++){
                rel_abundance[i] = 0.0;
                observed[i] = observed[i-1] + 1.0 / (gigp->S + 1.0);//  random_float_zero_to_x(1.0);
        }
        qsort(observed, (int)gigp->S, sizeof(double), compare_double);

        /*for(i = 0; i < (int) (gigp->S+1) - 30000;i++){
                observed[i] = -1;
                }*/
        int lower_bound;

        lower_bound = MACRO_MAX(0, (int) (gigp->S+1) - 30000);
        lower_bound = 0;
        rel_abundance = pick_abundances(gigp, observed, rel_abundance, lower_bound ,(int)gigp->S,DBL_MIN,1.0);
        qsort(rel_abundance, (int)gigp->S, sizeof(double), compare_double);


        return NULL;
ERROR:
        return NULL;
}

int fit_curves_make_growth_curve(struct shared_data* bsd, struct double_matrix* m)
{
        char buffer[BUFFER_LEN];

        double* fitted = NULL;
        FILE* out_ptr = NULL;
        double sum = 0.0;
        int i,j;

        ASSERT(m != NULL,"No data.");
        ASSERT(bsd != NULL, "No shared data.");

        for(i = 0 ;i < m->ncol;i++){
                LOG_MSG("Fitting curve to sample:%s.",m->col_names[i]);
                /* read in parameters for sample i..  */
                RUN(read_gigp_param_from_table(bsd->gigp_param_out_table,bsd->gigp_param_best,i));

                /* Turn counts into frequency vector. */
                RUN(get_unique_transcript_vector(bsd,m,i));

                /* make fitted and after write estimates S back into tables...  */
                MMALLOC(fitted,sizeof(double) * bsd->utc->max);
                RUN(fill_fitted_curve(fitted ,bsd->utc->max, bsd->gigp_param_best));

                /* enter_estimated parameters in table..   */
                RUN(enter_gigp_param_into_table(bsd->gigp_param_out_table,bsd->gigp_param_best,i));

                /* print fitted file...  */
                snprintf(buffer, BUFFER_LEN, "%s/%s/GIGP_fit_%s.csv",bsd->param->outdir,OUTDIR_MODEL,m->col_names[i]);
                LOG_MSG("Write to file: %s",buffer);
                RUNP(out_ptr = fopen(buffer, "w" ));

                for(j = 1; j <  bsd->utc->max;j++){
                        if( bsd->utc->x[j]){
                                fprintf(out_ptr,"%d,%f,%f\n",  j,  bsd->utc->x[j],fitted[j]);
                        }
                }
                fclose(out_ptr);
                MFREE(fitted);

                /* print growth curve */
                snprintf(buffer, BUFFER_LEN, "%s/%s/GIGP_growth_curve_%s.csv",bsd->param->outdir,OUTDIR_MODEL,m->col_names[i]);
                LOG_MSG("Writing to: %s",buffer);
                RUNP(out_ptr = fopen(buffer, "w" ));
                sum = bsd->gigp_param_best->N;

                for(j = 1000000; j <= 50000000;j+= 1000000){
                        if(j < sum && (double)(j + 1000000.0) > sum){
                                bsd->gigp_param_best->N = sum;
                                fprintf(out_ptr,"%d,%d,%d\n",(int)sum, (int)(bsd->gigp_param_best->S * ( 1.0 - give_me_sichel_p0(bsd->gigp_param_best,bsd->utc->x))),(int) bsd->gigp_param_best->S);
                        }
                        bsd->gigp_param_best->N = j;
                        fprintf(out_ptr,"%d,%d,%d\n",j, (int)(bsd->gigp_param_best->S * ( 1.0 - give_me_sichel_p0(bsd->gigp_param_best,bsd->utc->x))),(int) bsd->gigp_param_best->S);
                }

                fclose(out_ptr);
                bsd->gigp_param_best->N = sum;



                /* Done  */
                LOG_MSG("Done.");
        }

        return OK;
ERROR:
        return FAIL;
}


int read_gigp_param_table_from_file(struct shared_data* bsd)
{
        char buffer[BUFFER_LEN];

        ASSERT(bsd != NULL,"No shared data");

        if(bsd->gigp_param_out_table){
                WARNING_MSG("GIGP param table not empty will delete and load from file...");
                free_double_matrix(bsd->gigp_param_out_table);
        }
        snprintf(buffer, BUFFER_LEN, "%s/%s/%s",bsd->param->outdir,OUTDIR_MODEL,"GIGP_parameters.csv");
        RUNP(bsd->gigp_param_out_table =read_double_matrix(buffer,1,1));

        return OK;
ERROR:
        return FAIL;
}

int write_gigp_param_table_to_file(struct shared_data* bsd)
{
        char buffer[BUFFER_LEN];
        FILE* file_ptr = NULL;
        ASSERT(bsd != NULL, "No shared parameters.");
        snprintf(buffer, BUFFER_LEN, "%s/%s/%s",bsd->param->outdir,OUTDIR_MODEL,"GIGP_parameters.csv");
        RUNP(file_ptr = fopen(buffer,"w"));
        RUN(print_double_matrix(bsd->gigp_param_out_table,file_ptr ,1,1));

        fclose(file_ptr);
        return OK;
ERROR:
        return FAIL;
}

int enter_gigp_param_into_table(struct double_matrix* table, struct gigp_param* best,int col)
{
        ASSERT(best != NULL,"No parameters to write.");
        ASSERT(table != NULL,"No table.");
        ASSERT(col < table->ncol,"Too few columns in parameter table.");

        table->matrix[GIGP_GAMMA_ROW][col] = best->gamma;
        table->matrix[GIGP_B_ROW][col] = best->b;
        table->matrix[GIGP_C_ROW][col] = best->c;
        table->matrix[GIGP_N_ROW][col] = best->N;
        table->matrix[GIGP_s_ROW][col] = best->s;
        table->matrix[GIGP_S_ROW][col] = best->S;
        table->matrix[GIGP_MAXCOUNT_ROW][col] = best->max_count;
        table->matrix[GIGP_FIT_ROW][col] = best->fit;

        return OK;
ERROR:
        return FAIL;
}

int read_gigp_param_from_table(struct double_matrix* table, struct gigp_param* best,int col)
{
        ASSERT(best != NULL,"No parameters to write.");
        ASSERT(table != NULL,"No table.");
        ASSERT(col < table->ncol,"Too few columns in parameter table.");

        best->gamma = table->matrix[GIGP_GAMMA_ROW][col];
        best->b = table->matrix[GIGP_B_ROW][col];
        best->c = table->matrix[GIGP_C_ROW][col];
        best->N = table->matrix[GIGP_N_ROW][col];
        best->s = table->matrix[GIGP_s_ROW][col];
        best->S = table->matrix[GIGP_S_ROW][col];
        best->max_count = table->matrix[GIGP_MAXCOUNT_ROW][col];
        best->fit = table->matrix[GIGP_FIT_ROW][col];
        return OK;
ERROR:
        return FAIL;
}


int init_gigp_param_out_table(struct shared_data* bsd, struct double_matrix* input)
{
        int i;
        ASSERT(bsd != NULL,"No shared data.");
        bsd->gigp_param_out_table = NULL;
        RUNP(bsd->gigp_param_out_table = alloc_double_matrix(input->ncol,8,BUFFER_LEN));

        for(i = 0; i< input->ncol;i++){
                snprintf(bsd->gigp_param_out_table->col_names[i],BUFFER_LEN, "%s",input->col_names[i]);
        }

        snprintf(bsd->gigp_param_out_table->row_names[GIGP_GAMMA_ROW],BUFFER_LEN,"Gamma");
        snprintf(bsd->gigp_param_out_table->row_names[GIGP_B_ROW],BUFFER_LEN,"B");
        snprintf(bsd->gigp_param_out_table->row_names[GIGP_C_ROW],BUFFER_LEN,"C");
        snprintf(bsd->gigp_param_out_table->row_names[GIGP_N_ROW],BUFFER_LEN,"N");
        snprintf(bsd->gigp_param_out_table->row_names[GIGP_s_ROW],BUFFER_LEN,"s");
        snprintf(bsd->gigp_param_out_table->row_names[GIGP_S_ROW],BUFFER_LEN,"S");
        snprintf(bsd->gigp_param_out_table->row_names[GIGP_MAXCOUNT_ROW],BUFFER_LEN,"max");
        snprintf(bsd->gigp_param_out_table->row_names[GIGP_FIT_ROW],BUFFER_LEN,"Fit");
        return OK;
ERROR:
        return FAIL;
}


int sample_GIGP_same_depth_as_input(struct shared_data* bsd)
{
        char buffer[BUFFER_LEN];
        FILE* file_ptr = NULL;

        struct gigp_param* gigp_param = NULL;
        double* observed = NULL;
        double* rel_abundance = NULL;
        double* fx_unique_transcript_count = NULL;
        int i,c,f,g;
        int max_transcript_frequency = 0;
        int t_count_size = 10000;

        /* red in relative abundancs */

        ASSERT(bsd != NULL,"no shared data.");
        gigp_param = bsd->gigp_param_best;
       	MMALLOC(observed, sizeof(double) * (int) (gigp_param->S+1));

        MMALLOC(rel_abundance, sizeof(double) * (int) (gigp_param->S+1));

        snprintf(buffer, BUFFER_LEN, "%s/%s/%s",bsd->param->outdir,OUTDIR_MODEL,"GIGP_rel_abundances.txt");

        RUNP(file_ptr = fopen(buffer, "r"));
        for(i = 0; i < (int)(gigp_param->S+1) ;i++){
                if(fscanf(file_ptr,"%lg\n",  &rel_abundance[i]) != 1){
                        ERROR_MSG("Could not read rel abundance %d.",i);
                }
        }
        fclose(file_ptr);

        gsl_rng_env_setup();
        const gsl_rng_type *T = gsl_rng_default;
        gsl_rng *rfff = gsl_rng_alloc (T);

        gsl_ran_discrete_t *sample_abundances = gsl_ran_discrete_preproc ((int)gigp_param->S, rel_abundance);

        for(i = 0; i <(int)gigp_param->S;i++){
                observed[i] = 0.0;
        }

        for(i = 1; i <= 14000000 ;i++){
                observed[(int)gsl_ran_discrete(rfff, sample_abundances)] += 1.0;
                if(i > 1 && (i == gigp_param->N  || (i % 1000000) == 0)){
                        max_transcript_frequency = 0;
                        t_count_size = 10000;
                        MMALLOC(fx_unique_transcript_count, sizeof(double)* t_count_size);
                        for(c= 0; c < t_count_size;c++){
                                fx_unique_transcript_count[c] = 0.0;
                        }
                        f = 0;
                        for(c = 0;c < (int)gigp_param->S;c++){
                                if(observed[c] != 0.0){
                                        while( (int)observed[c] >= t_count_size){
                                                fx_unique_transcript_count = realloc(fx_unique_transcript_count, sizeof(double)* t_count_size *2);
                                                for(g = t_count_size;g < t_count_size*2 ;  g++ ){
                                                        fx_unique_transcript_count[g] = 0.0;
                                                }
                                                t_count_size = t_count_size *2;
                                        }
                                        if(max_transcript_frequency <  (int)  observed[c]){
                                                max_transcript_frequency = (int) observed[c];
                                        }
                                        fx_unique_transcript_count[ (int)observed[c]] += 1.0;
                                }else{
                                        f++;
                                }
                        }
                        snprintf(buffer, BUFFER_LEN, "%s/%s/GIGP_sampled_%d.csv",bsd->param->outdir,OUTDIR_MODEL,i);
                        RUNP(file_ptr = fopen(buffer,"w"));
                        for(c = 1;c <= max_transcript_frequency;c++){
                                if(fx_unique_transcript_count[c]){
                                        fprintf(file_ptr,"%d,%f\n",  c,fx_unique_transcript_count[c]);
                                }
                        }
                        fclose(file_ptr);
                        MFREE(fx_unique_transcript_count);
                }
        }
        gsl_rng_free(rfff);
        gsl_ran_discrete_free(sample_abundances);
        return OK;
ERROR:
        return FAIL;
}

int simulate_expression_table(struct shared_data* bsd)
{
        struct double_matrix* m = NULL;
        char buffer[BUFFER_LEN];
        FILE* in_ptr = NULL;
        double* rel_abundance = NULL;
        double* cell_abundance = 0;
        double CV = 0.0f;
        double libdepth = 0.0f;
        double drop_k = 0.0f;
        int i,j;
        int num_genes;
        int num_cells;
        int simulated_read_depth = 0;

        ASSERT(bsd != NULL,"No shared data.");

        snprintf(buffer, BUFFER_LEN, "%s/%s/%s",bsd->param->outdir,OUTDIR_MODEL,"GIGP_parameters.txt");
        RUN(read_fitted_param_from_file( bsd->gigp_param_best,buffer));

        bsd->param->num_genes = (int)(bsd->gigp_param_best->S +1);

        drop_k = bsd->param->drop_Michaelisconstant;
        num_genes = bsd->param->num_genes;
        num_cells = bsd->param->num_cells;
        simulated_read_depth = bsd->param->simulated_read_depth;
        CV = bsd->param->CV;
        /* red in relative abundancs */

        MMALLOC(rel_abundance, sizeof(double) * (int) num_genes);
        MMALLOC(cell_abundance, sizeof(double) * num_genes);

        snprintf(buffer, BUFFER_LEN, "%s/%s/%s",bsd->param->outdir,OUTDIR_MODEL,"GIGP_rel_abundances.txt");

        RUNP(in_ptr = fopen(buffer, "r"));
        for(i = 0; i < num_genes ;i++){
                if(fscanf(in_ptr,"%lg\n",  &rel_abundance[i]) != 1){
                        ERROR_MSG("Could not read rel abundance %d.",i);
                }
        }
        fclose(in_ptr);

        /* set up random number generator  */
        gsl_rng_env_setup();
        const gsl_rng_type *T = gsl_rng_default;
        gsl_rng *rfff = gsl_rng_alloc (T);
        LOG_MSG("Allocating expression table: %d %d",bsd->param->num_cells,bsd->param->num_genes);
        /* alloc outout table and aux. vector  */
        RUNP(m = alloc_double_matrix(bsd->param->num_cells,bsd->param->num_genes,15));
        /* fill fake names */
        for(i = 0; i < num_cells ;i++){
                snprintf(m->col_names[i],15,"Cell%d",i+1);
        }
        for(j = 0; j < num_genes;j++){
                snprintf(m->row_names[j],15,"Gene%d",j+1);
        }


        /* fill expression table  */

        //Calculate tpm..
        libdepth = 0.0;
        for(j = 0; j < num_genes;j++){
                libdepth += rel_abundance[j];
        }

        for(i = 0; i < num_cells ;i++){
                for(j = 0; j < num_genes;j++){
                        cell_abundance[j] = rel_abundance[j] + gsl_ran_gaussian (rfff,  rel_abundance[j] * CV);
                        if(cell_abundance[j] < 0){
                                cell_abundance[j] = 0.0;
                        }
                        double average = rel_abundance[j] / libdepth * 1000000.0;
                        double p_drop = 1.0 - (average) / (average + drop_k);
                        if(random_float_zero_to_x(1.0) <= p_drop){
                                cell_abundance[j] = 0.0;
                        }

                }
                gsl_ran_discrete_t *sample_abundances = gsl_ran_discrete_preproc (num_genes, cell_abundance);


                for(j = 0;j < simulated_read_depth;j++){
                        m->matrix[(int)gsl_ran_discrete(rfff, sample_abundances)][i] += 1;
                }
                gsl_ran_discrete_free(sample_abundances);
        }
        //random deletion...
        /*for(i = 0; i < ncells;i++){
          for(j = 0; j < num_genes;j++){
          if(gsl_rng_uniform_pos(rfff) < loss){
          table[i][j] = 0;
          }
          }
          }*/


        snprintf(buffer, BUFFER_LEN, "%s/%s/extable_%d_%d_CV%f.csv",bsd->param->outdir,OUTDIR_MODEL,num_cells, simulated_read_depth,CV);
        RUNP(in_ptr = fopen(buffer, "w"));
        RUN(print_double_matrix(m,in_ptr, 1,1));
        fclose(in_ptr);

        gsl_rng_free(rfff);
        MFREE(cell_abundance);
        MFREE(rel_abundance);

        free_double_matrix(m);

        return OK;

ERROR:
        if(m){
                free_double_matrix(m);
        }
        if(rel_abundance){
                MFREE(rel_abundance);
        }
        return FAIL;
}

int calculate_rel_frequencies(struct shared_data* bsd)
{
        char buffer[BUFFER_LEN];
        FILE* out_ptr = NULL;
        double* rel_abundance = NULL;
        double* observed = NULL;
        int i;

        ASSERT(bsd != NULL,"No shared data.");


        snprintf(buffer, BUFFER_LEN, "%s/%s/%s",bsd->param->outdir,OUTDIR_MODEL,"GIGP_parameters.txt");
        RUN(read_fitted_param_from_file( bsd->gigp_param_best,buffer));

        MMALLOC(rel_abundance, sizeof(double) * (int) (bsd->gigp_param_best->S+1));
        MMALLOC(observed, sizeof(double) * (int)(bsd->gigp_param_best->S+1));

        for(i = 0; i <(int) (bsd->gigp_param_best->S+1);i++){
                rel_abundance[i] = 0.0;
                observed[i] =  random_float_zero_to_x(1.0);
        }

        qsort(observed, (int)bsd->gigp_param_best->S, sizeof(double), compare_double);

        rel_abundance = pick_abundances(bsd->gigp_param_best, observed, rel_abundance, 0 ,(int)bsd->gigp_param_best->S,1e-7,1);

        qsort(rel_abundance, (int)bsd->gigp_param_best->S, sizeof(double), compare_double);


        snprintf(buffer, BUFFER_LEN, "%s/%s/%s",bsd->param->outdir,OUTDIR_MODEL,"GIGP_rel_abundances.txt");

        RUNP(out_ptr = fopen(buffer, "w"));
        for(i = 0; i <(int)(bsd->gigp_param_best->S+1);i++){
                fprintf(out_ptr,"%10.10e\n",  rel_abundance[i]  );
        }
        fclose(out_ptr);

        //sampling...

        MFREE(rel_abundance);
        MFREE(observed);
        return OK;
ERROR:
        if(rel_abundance){
                MFREE(rel_abundance);
        }
        if(observed){
                MFREE(observed);
        }
        return FAIL;
}

int print_fit(struct shared_data* bsd)
{
        char buffer[BUFFER_LEN];
        FILE* out_ptr = NULL; 	/* is actually an input...  */
        double* fitted = NULL;
        double sum;

        int i;

        ASSERT(bsd != NULL,"No shared memory.");

        snprintf(buffer, BUFFER_LEN, "%s/%s/%s",bsd->param->outdir,OUTDIR_MODEL,"GIGP_parameters.txt");

        RUN(read_fitted_param_from_file(bsd->gigp_param_best,buffer));

        MMALLOC(fitted,sizeof(double) * bsd->utc->max);
        RUN(fill_fitted_curve(fitted ,bsd->utc->max, bsd->gigp_param_best));

        /* print fitted file...  */
        snprintf(buffer, BUFFER_LEN, "%s/%s/%s",bsd->param->outdir,OUTDIR_MODEL,"GIGP_fit.csv");
        LOG_MSG("REading from: %s",buffer);

        RUNP(out_ptr = fopen(buffer, "w" ));
        for(i = 1; i <  bsd->utc->max;i++){
                if( bsd->utc->x[i]){
                        fprintf(out_ptr,"%d,%f,%f\n",  i,  bsd->utc->x[i],fitted[i]);
                }
        }
        fclose(out_ptr);

        snprintf(buffer, BUFFER_LEN, "%s/%s/%s",bsd->param->outdir,OUTDIR_MODEL,"GIGP_growth_curve.csv");
        LOG_MSG("Writing to: %s",buffer);
        RUNP(out_ptr = fopen(buffer, "w" ));
        sum = bsd->gigp_param_best->N;

        for(i = 1000000; i <= 50000000;i+= 1000000){
                if(i < sum && (double)(i + 1000000.0) > sum){
                        bsd->gigp_param_best->N = sum;
                        fprintf(out_ptr,"%d,%d,%d\n",(int)sum, (int)(bsd->gigp_param_best->S * ( 1.0 - give_me_sichel_p0(bsd->gigp_param_best,bsd->utc->x))),(int) bsd->gigp_param_best->S);
                }
                bsd->gigp_param_best->N = i;
                fprintf(out_ptr,"%d,%d,%d\n",i, (int)(bsd->gigp_param_best->S * ( 1.0 - give_me_sichel_p0(bsd->gigp_param_best,bsd->utc->x))),(int) bsd->gigp_param_best->S);
        }

        fclose(out_ptr);
        bsd->gigp_param_best->N = sum;

        snprintf(buffer, BUFFER_LEN, "%s/%s/%s",bsd->param->outdir,OUTDIR_MODEL
                 ,"GIGP_parameters.txt");
        RUN(write_fitted_param_to_file( bsd->gigp_param_best,buffer));

        MFREE(fitted);
        return OK;
ERROR:
        if(fitted){
                MFREE(fitted);
        }
        return FAIL;
}

int set_initial_guess(struct shared_data* bsd)
{
        struct gigp_param* ptr = NULL;
        int i;

        ASSERT(bsd != NULL,"No shared data.");

        ptr = bsd->gigp_param_initial_guess;

        ptr->gamma = -2.0;
        ptr->b = 1.0;
        ptr->c = 0.1;
        ptr->N = 0.0;
        ptr->s = 0.0;
        ptr->S = 0.0;
        ptr->fit = DBL_MAX;
        ptr->max_count = bsd->utc->max;

        for(i = 0; i < bsd->utc->max;i++){
                ptr->N += i * bsd->utc->x[i];
                ptr->s +=  bsd->utc->x[i];
        }

        return OK;
ERROR:
        return FAIL;
}

int get_unique_transcript_vector(struct shared_data* bsd, struct double_matrix* m, int column)
{
        struct unique_transcript_count* utc = NULL;
        int i,c;
        ASSERT(bsd != NULL,"No shared data.");
        ASSERT(m != NULL,"Count matrix is NULL.");

        ASSERT(column < m->nrow,"Count matrix has %d columns but you ask for %d.", m->nrow,column );

        utc = bsd->utc;

        utc->max = 0;
        for(i = 0; i < utc->alloc_size;i++){
                utc->x[i] = 0.0;
        }
        for(i = 0; i < m->nrow;i++){
                if(m->matrix[i][column] != 0){
                        while( (int)m->matrix[i][column]  >= utc->alloc_size){
                                utc->x = realloc(utc->x, sizeof(double)*  (utc->alloc_size +64));
                                for(c =  utc->alloc_size;c <  (utc->alloc_size+64) ;  c++ ){
                                        utc->x[c] = 0.0;
                                }
                                utc->alloc_size = (utc->alloc_size +64);
                        }
                        if(utc->max < (int) m->matrix[i][column] ){
                                utc->max = (int) m->matrix[i][column];
                        }
                        utc->x[(int) m->matrix[i][column]] += 1.0;
                }
        }
        return OK;
ERROR:
        return FAIL;
}


int fit_model(struct  shared_data* bsd)
{

        struct gigp_param* gigp_param = NULL;

        struct thread_data** td = NULL;

        int alloc_td = 1000;

        int step = 3;
        double min_gamma = -3.0;
        double max_gamma = 3.0;

        double min_b = 0.1;
        double max_b = 0.4;

        double min_c = 0.01;
        double max_c = 0.04;


        int status;
        int i,j,c,g;

        g = 0;

        ASSERT(bsd != NULL,"No shared data.");

        /* set inital guesss */
        RUN(set_initial_guess(bsd));

        MMALLOC(td, sizeof(struct thread_data*) * alloc_td);

        for(i = 0 ; i < step; i++){
                for(j = 0; j < step; j++){
                        for(c = 0; c < step; c++){
                                td[g] = NULL;
                                MMALLOC(td[g], sizeof(struct thread_data));
                                td[g]->bsd = bsd;
                                td[g]->thread_id = 0;
                                td[g]->num_threads = bsd->param->num_threads;
                                td[g]->try_gamma = min_gamma + ((max_gamma - min_gamma  )/(double)step *(double)i);
                                td[g]->try_b  = min_b + ((max_b -min_b  )/(double)step *(double)j);
                                td[g]->try_c  = min_c + ((max_c -min_c  )/(double)step *(double)c);

                                RUNP(gigp_param = init_gigp_param());
                                td[g]->gigp_param = gigp_param;
                                gigp_param = NULL;

                                g++;
                                if(g == alloc_td){
                                        alloc_td = alloc_td << 1;
                                        MREALLOC(td,sizeof(struct thread_data*) * alloc_td);
                                }
                        }
                }
        }
        for(i = 0; i < g;i++){
                td[i]->thread_id = i;
                if((status = thr_pool_queue(bsd->pool,do_fit_model,td[i])) == -1) fprintf(stderr,"Adding job to queue failed.");
        }
        thr_pool_wait(bsd->pool);

        for(i = 0; i < g;i++){
                MFREE(td[i]->gigp_param);
                MFREE(td[i]);
        }
        MFREE(td);
        return OK;
ERROR:
        for(i = 0; i < g;i++){
                MFREE(td[i]->gigp_param);
                MFREE(td[i]);
        }
        MFREE(td);
        return FAIL;
}

int write_fitted_param_to_file(struct gigp_param* gigp_param, char* filename)
{
        FILE* out_ptr = NULL;

        RUNP(out_ptr = fopen(filename, "w" ));
        fprintf(out_ptr,"%10.10e\tgamma\n",  gigp_param->gamma);
        fprintf(out_ptr,"%10.10e\tb\n",  gigp_param->b );
        fprintf(out_ptr,"%10.10e\tc\n",  gigp_param->c );
        fprintf(out_ptr,"%f\tN\n", gigp_param->N );
        fprintf(out_ptr,"%f\ts\n", gigp_param->s);
        fprintf(out_ptr,"%f\tS\n", gigp_param->S );
        fprintf(out_ptr,"%f\tMax_transcript_count\n", gigp_param->max_count -1 );
        fprintf(out_ptr,"%f\tFit\n", gigp_param->fit);
        fclose(out_ptr);
        return OK;
ERROR:
        return FAIL;
}

int read_fitted_param_from_file(struct gigp_param* gigp_param, char* filename)
{
        FILE* out_ptr = NULL;

        LOG_MSG("Reading from: %s",filename);

        RUNP(out_ptr = fopen(filename, "r" ));

        if(fscanf(out_ptr,"%lg\tgamma\n",  &gigp_param->gamma )!= 1){
                ERROR_MSG("could not read gamma from %s.",filename);
        }
        if(fscanf(out_ptr,"%lg\tb\n",  &gigp_param->b )!= 1){
                ERROR_MSG("could not read b %s.",filename);
        }
        if(fscanf(out_ptr,"%lg\tc\n",  &gigp_param->c )!= 1){
                ERROR_MSG("could not read c from %s.",filename);
        }
        if(fscanf(out_ptr,"%lf\tN\n", &gigp_param->N )!= 1){
                ERROR_MSG("could not read N from %s.",filename);
        }
        if(fscanf(out_ptr,"%lf\ts\n",  &gigp_param->s)!= 1){
                ERROR_MSG("could not read s from %s.",filename);
        }
        if(fscanf(out_ptr,"%lf\tS\n", &gigp_param->S )!= 1){
                ERROR_MSG("could not read S from %s.",filename);
        }
        if(fscanf(out_ptr,"%lf\tMax_transcript_count\n", &gigp_param->max_count)!= 1){
                ERROR_MSG("could not read max count from %s.",filename);
        }
        if(fscanf(out_ptr,"%lf\tFit\n", &gigp_param->fit)!= 1){
                ERROR_MSG("could not read fir from %s.",filename);
        }
        fclose(out_ptr);

        gigp_param->max_count++;
        return OK;
ERROR:
        return FAIL;
}


void* do_fit_model(void *threadarg)
{
        struct thread_data *data = NULL;
        struct gigp_param* gigp = NULL;

        data = (struct thread_data *) threadarg;

        gigp = data->gigp_param;//  bsd->working_gigp_param[id];
        /* copy starting parameters from shares  */

        RUN(copy_gigp_param(data->bsd->gigp_param_initial_guess,gigp));

        /* copy run specific params...  */

        gigp->gamma = data->try_gamma;
        gigp->b = data->try_b;
        gigp->c = data->try_c;
        RUN(fit_gigp(gigp,data->bsd->utc->x, data->bsd->utc->max));

        pthread_mutex_lock(&data->bsd->avail_mtx);
        if(gigp->fit < data->bsd->gigp_param_best->fit){
                RUN(copy_gigp_param(gigp,data->bsd->gigp_param_best));
        }


        //data->bsd->available_work_space[id] = 1;
        pthread_mutex_unlock(&data->bsd->avail_mtx);

        return NULL;
ERROR:
        return NULL;
}

struct shared_data* init_shared_data(struct parameters* param)
{
        struct shared_data* bsd = NULL;
        struct gigp_param* tmp = NULL;
        int i;

        ASSERT(param != NULL,"No parameters");

        MMALLOC(bsd, sizeof(struct shared_data));
        bsd->param = param;
        bsd->utc = NULL;
        bsd->gigp_param_out_table = NULL;
        bsd->rel_abundances = NULL;
        bsd->gigp_param_best  = NULL;
        bsd->gigp_param_initial_guess  = NULL;

        bsd->working_gigp_param = NULL;
        bsd->num_threads = param->num_threads;
        bsd->available_work_space = NULL;
        MMALLOC(bsd->available_work_space, sizeof(uint8_t) *bsd->num_threads);
        for(i = 0; i < bsd->num_threads;i++){
                bsd->available_work_space[i] = 1;
        }

        pthread_mutex_init(&bsd->avail_mtx,NULL);

        bsd->pool = NULL;
        if((bsd->pool = thr_pool_create(param->num_threads+1, param->num_threads+1, 0, 0)) == NULL) ERROR_MSG("Creating pool thread failed.");

        RUNP(tmp = init_gigp_param());
        bsd->gigp_param_best = tmp;
        tmp = NULL;

        RUNP(tmp = init_gigp_param());
        bsd->gigp_param_initial_guess = tmp;
        tmp = NULL;

        MMALLOC(bsd->utc, sizeof(struct unique_transcript_count));
        bsd->utc->x = NULL;
        bsd->utc->alloc_size = 16384;
        bsd->utc->max = 0;

        MMALLOC(bsd->utc->x, sizeof(double)* bsd->utc->alloc_size);


        MMALLOC(bsd->working_gigp_param,sizeof(struct gigp_param*)* bsd->num_threads);
        for(i = 0; i < bsd->num_threads;i++){

                bsd->available_work_space[i] = 1;

                bsd->working_gigp_param[i] = NULL;

                RUNP(tmp = init_gigp_param());
                bsd->working_gigp_param[i]  = tmp;
                tmp = NULL;
        }

        return bsd;
ERROR:
        return NULL;
}

void free_shared_data(struct shared_data* bsd)
{
        int i;
        if(bsd){
                for(i = 0; i < bsd->num_threads;i++){
                        MFREE(bsd->working_gigp_param[i]);
                }
                MFREE(bsd->working_gigp_param);
                MFREE(bsd->gigp_param_initial_guess);
                MFREE(bsd->gigp_param_best);
                if(bsd->utc){
                        MFREE(bsd->utc->x);
                        MFREE(bsd->utc);
                }
                pthread_mutex_destroy(&bsd->avail_mtx);

                if(bsd->available_work_space){
                        MFREE(bsd->available_work_space);
                }
                if(bsd->rel_abundances){
                        free_double_matrix(bsd->rel_abundances);
                }

                if(bsd->gigp_param_out_table){
                        free_double_matrix(bsd->gigp_param_out_table);
                }

                if(bsd->pool){
                        thr_pool_destroy(bsd->pool);
                }
                MFREE(bsd);
        }
}


int fill_fitted_curve(double* fit,int num,struct gigp_param* param )
{

        double alpha = 0.0;
        double beta = 0.0;

        double gamma = 0.0;
        double res[3];


        double y = 0;
        gsl_sf_result result;

        double x = 0;
        double sum = 0.0;
        double bessel1 = 0.0;
        double bessel2 = 0.0;
        int status = 0;

        ASSERT(param != NULL,"No parameters.");

        alpha = param->b*sqrt(1.0 + param->c *param->N );
        beta = (param->b*param->c*param->N)/(2.0 * sqrt(1.0 + param->c *param->N ));

        gamma = param->gamma;

        x = 0;


        status = gsl_sf_bessel_Kn_e( gamma, sqrt( pow(alpha,2.0) - 2.0* alpha*beta )   ,  &result);

        if (status == GSL_SUCCESS) {
                //printf("%f	%f\n",i,result.val);
                bessel1 =  result.val;
        }else{
                ERROR_MSG("gsl_sf_bessel_Kn_e failed with param: %f %f.",gamma, sqrt( pow(alpha,2.0) - 2.0* alpha*beta));
        }

        status = gsl_sf_bessel_Kn_e( gamma + x,  alpha   ,  &result);

        if (status == GSL_SUCCESS) {
                //printf("%f	%f\n",i,result.val);
                bessel2 =  result.val;
        }else{
                ERROR_MSG("gsl_sf_bessel_Kn_e failed with param: %f %f.",gamma + x,  alpha );
        }






        y = pow((  sqrt( pow(alpha,2.0) - 2.0* alpha*beta  )   / alpha  ),gamma) * 1.0 / bessel1 * (pow(beta,x)/  1.0 )  * bessel2;
        fit[(int)x] =y;
        sum += y;
        //fprintf(stderr,"%f	%f	%f	%f	%f\n",x, fit[(int)x ],alpha,beta,gamma );
        //fprintf(stdout,"%f	%e	%f\n", x,y,y * 77152.0);
        res[0] = y;
        //fprintf(stderr,"S:%f	%f\n", (51470.0 / (1.0- y)),y) ;

        x = 1.0 ;
        status = gsl_sf_bessel_Kn_e( gamma, sqrt( pow(alpha,2.0) - 2.0* alpha*beta )   ,  &result);
        if (status == GSL_SUCCESS) {
                bessel1 =  result.val;
        }else{
                ERROR_MSG("gsl_sf_bessel_Kn_e failed with param: %f %f.",gamma, sqrt( pow(alpha,2.0) - 2.0* alpha*beta));
        }

        status = gsl_sf_bessel_Kn_e( gamma + x,  alpha   ,  &result);
        if (status == GSL_SUCCESS) {
                bessel2 =  result.val;
        }else{
                ERROR_MSG("gsl_sf_bessel_Kn_e failed with param: %f %f.",gamma + x,  alpha );
        }




        y = pow ((  sqrt( pow(alpha,2.0) - 2.0* alpha*beta  )   / alpha  ),gamma) * 1.0 / bessel1 * (pow(beta,x)/  1.0 )  * bessel2;
        fit[(int)x] =y;
        sum += y;

        //fprintf(stdout,"%f	%e	%f\n", x,y,y * 77152.0);
        res[1] = y;

        for(x = 2.0; x <   param->max_count ;x+=1.0){
                y = (2.0*beta)/ alpha * ((x + gamma -1.0)/ x) * res[1] + (pow(beta,2.0)/ ( x *(x-1.0) )  )* res[0];
                //fprintf(stdout,"%f	%e	%f\n", x,y,y * 77152.0);
                fit[(int)x] =y;
                sum += y;

                res[0] = res[1];
                res[1] = y;


        }
        //fprintf(stderr,"Fit0:%f\n", fit[0]);
        param->S = param->s * 1.0 / ( 1.0 - fit[0] / sum );


        for(x = 0; x < param->max_count;x++){
                fit[(int)x] = fit[(int)x] / sum * param->S;
        }

        return OK;
ERROR:
        return FAIL;
}


static int compare_double (const void * a, const void * b)
{
        if (*(double*)a > *(double*)b) return 1;
        else if (*(double*)a < *(double*)b) return -1;
        else return 0;
}


double* pick_abundances(struct gigp_param* gigp_param, double* random, double* abundances,int outer_low,int outer_high,double inner_low,double inner_high)
{
        if(outer_low > outer_high){
                return abundances;
        }
        gsl_function F;
        gsl_integration_workspace * w	= gsl_integration_workspace_alloc (1000);
        struct gigp_param param_4integration;
        param_4integration.gamma =  gigp_param->gamma;
        param_4integration.b = gigp_param->b;
        param_4integration.c = gigp_param->c;

        F.function = &gigp_dist;
        F.params = &param_4integration;
        double x,r;
        double integration_result, error,sum;
        //r =   (double)rand_r(&seed)/(double)RAND_MAX;
        x = 1.0;
        gsl_integration_qag (&F, 0.0, x,1e-12, 1e-16, 1000,GSL_INTEG_GAUSS61, w , &integration_result, &error);
        sum = integration_result;

        double inner_mid = 0.0;
        int outer_mid;

        outer_mid =  (outer_high + outer_low) / 2.0;

        double tmp[3];

        //if(outer_high    )

        tmp[0] = inner_low;
        tmp[1] = inner_mid;
        tmp[2] = inner_high;

        //fprintf(stdout,"L:%d H:%d m:%d workon\n",outer_low, outer_high, outer_mid);
        r = random[(int)outer_mid];
        if(r == -1.0){
                gsl_integration_workspace_free (w);
                return abundances;
        }
        while(1){
                inner_mid = (inner_high + inner_low )/ 2.0;
                gsl_integration_qag (&F, 0.0, inner_mid,1e-12, 1e-16, 1000,GSL_INTEG_GAUSS61,w , &integration_result, &error);
                integration_result = integration_result / sum;
                //fprintf(stderr,"%d L:%e\tH:%e\tmid:%e\t\t%f\t%f\n",outer_mid, inner_low,inner_high,inner_mid, integration_result, r  );

                if(fabs(integration_result - r) < DBL_EPSILON){
                        abundances[outer_mid] = inner_mid;
                        random[outer_mid] = -1.0;
                        break;
                }else if(integration_result < r){
                        inner_low = inner_mid;
                }else	if(integration_result > r ){
                        inner_high = inner_mid;
                }

                if(fabs(inner_high - inner_low)  < DBL_EPSILON){
                        abundances[outer_mid] = inner_mid;
                        random[outer_mid] = -1.0;
                        break;
                }
        }

        gsl_integration_workspace_free (w);

        /* if((outer_mid-1) -outer_low > 1000 ){ */
        /*         fprintf(stderr,"Working on : %d - %d\n",outer_low, outer_mid-1); */
        /* } */
        abundances = pick_abundances(gigp_param,random,abundances, outer_low, outer_mid-1,tmp[0],inner_mid);

        /* if(outer_high - (outer_mid+1)  > 1000){ */
        /*         fprintf(stderr,"Working on : %d - %d\n",outer_mid+1, outer_high); */
        /* } */
        abundances = pick_abundances(gigp_param,random,abundances, outer_mid+1, outer_high,  inner_mid,tmp[2]);

        return  abundances;
}




int cmp_gene_ex_struct(const void *a, const void *b)
{
        struct gene_ex* const *one = a;
        struct gene_ex* const *two = b;
        if ((*one)->expr >= (*two)->expr){
                return -1;
        }else{
                return 1;
        }
        return 1;
}
