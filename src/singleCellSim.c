#include "scs.h"
#include "outdir.h"		/*  */

#define OPT_SAMPLE 1
#define OPT_BOUNDARIES 2
#define OPT_DATFILE 3
#define OPT_NUMCELLS 4 
#define OPT_DEPTH 5
#define OPT_CV 6
#define OPT_LOSS 7

#define OPT_NUM_THREAD 8


struct gigp_param{
	double gamma;
	double b;
	double c;
	double N;
	double s;
	double S;
	double max_count;
	double fit;
};

struct unique_transcript_count{
	double* x;
	int max;
	int alloc_size;
};

struct parameters{
	char* infile;
	char* outdir;
	int num_threads;
};

struct shared_data{
	struct parameters* param;
	struct unique_transcript_count* utc;
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
	double try_gamma;
	double try_b;
	double try_c;
	int target;
	int thread_id;
	int num_threads;
};

int run_scs(struct parameters* param);

struct unique_transcript_count* get_unique_transcript_vector(struct float_matrix* m, int column);

int fit_model(struct unique_transcript_count* utc);

int fit_gigp_controller(struct unique_transcript_count* utc, struct gigp_param* gigp_param);

int fit_gigp(struct gigp_param* gigp_param,double* fx_unique_transcript_count, int len);

int set_initial_guess(struct shared_data* bsd);


struct shared_data* init_shared_data(struct parameters* param);
void free_shared_data(struct shared_data* bsd);

struct gigp_param* init_gigp_param(void);
int copy_gigp_param(struct gigp_param* source, struct gigp_param* target);


double give_me_sichel_p0(struct gigp_param* gigp_param ,double* p);
int make_fake_expression_table(double CV,double* rel_abundance,char* outfile, int ncells, int depth,int num_genes,double loss);

static int compare_double (const void * a, const void * b);


//struct gigp_param* fit_gigp(struct gigp_param* gigp_param,double* fx_unique_transcript_count, int len);
struct gigp_param* fit_exact_gigp(struct gigp_param* gigp_param,double* fx_unique_transcript_count, int len);
double gigp_dist (double x, void * p);


double* fill_fitted_curve(double* fit,int num,struct gigp_param* param );

double my_f (const gsl_vector *v, void *params);

double sichel_function(const gsl_vector *v, void *params);

double* read_count_distribution(int argc, char * argv[] ,int* max_transcript_frequency,int data_origin,char** sample_names,int num_samples);
double* pick_abundances(struct gigp_param* gigp_param, double* random, double* abundances,int outer_low,int outer_high,double inner_low,double inner_high);

int main (int argc, char * argv[])
{
	struct parameters* param = NULL;
	
	tlog.echo_build_config();
     	MMALLOC(param, sizeof(struct parameters));
	param->infile = NULL;
	param->outdir = NULL;
	param->num_threads = 8;
	
	FILE* file = 0;
	int i,c,f,g;
	char* tmp = 0;
	char* outfile_name = 0;
	int count;
	
	int num_samples = 0;
	char** sample_names = 0;
	int help = 0;
	int quiet = 0;
	char* outfile = 0;
	int number_of_cells = 0;
	int simulated_read_depth = 0;
	double CV = 0.75;
	double loss = 0.00;
	double* rel_abundance = 0;
	double* observed = 0;
	double sum;
	
	int data_origin = 0;
	
	unsigned int seed = (unsigned int) (time(NULL) * ( 42));
	
	while (1){
		static struct option long_options[] ={
			{"sample",required_argument,0,OPT_SAMPLE},
			{"in",required_argument,0, 'i'},
			{"out",required_argument,0, 'o'},
			{"nthread",required_argument,0,OPT_NUM_THREAD},			
			{"boundaries",0,0, OPT_BOUNDARIES},
			{"datafile",0,0, OPT_DATFILE},
			{"ncell",required_argument,0,OPT_NUMCELLS},
			{"depth",required_argument,0,OPT_DEPTH},
			{"loss",required_argument,0,OPT_LOSS },
			{"CV",required_argument,0,OPT_CV},
			{"quiet",0,0,'q'},
			{"help",0,0,'h'},
			{0, 0, 0, 0}
		};
		
		int option_index = 0;
		c = getopt_long_only (argc, argv,"o:qh:",long_options, &option_index);
		
		if (c == -1){
			break;
		}
		
		switch(c) {
		case 0:
			break;
		case OPT_LOSS:
			loss = atof (optarg);
			break;
		case OPT_CV:
			CV = atof(optarg);
			break;
		case OPT_NUMCELLS:
			number_of_cells = atoi(optarg);
			break;
		case OPT_DEPTH:
			simulated_read_depth = atoi(optarg);
			break;
		case OPT_BOUNDARIES:
			data_origin = 1;
			break;
		case OPT_DATFILE:
			data_origin = 3;
			break;
		case OPT_SAMPLE:
			tmp = optarg;
			//count = byg_count(",", tmp);
			//fprintf(stderr,"%d\n",count);
			num_samples = count+1;
				
				
			MMALLOC(sample_names,sizeof(char*) * (count+1));
				
			for(i = 0; i < count+1;i++){
				sample_names[i] = 0;
				MMALLOC(sample_names[i],sizeof(char)* 1000);
			}
			f = 0;
			g = 0;
			for(i = 0;i < strlen(tmp);i++){
				if(tmp[i] != ','){
					sample_names[f][g] = tmp[i];
					g++;
				}else{
					sample_names[f][g] = 0;
					f++;
					g = 0;
				}
			}
			sample_names[f][g] = 0;
				

			break;
		case OPT_NUM_THREAD:
			param->num_threads = atoi(optarg);
			break;					
		case 'q':
			quiet = 1;
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
			fprintf(stderr,"default\n\n\n\n");
			abort ();
		}
	}
	
	
	
	
	if(num_samples > 1){
		fprintf(stderr,"Warning: I will only use the first sample... \n");
	}
	
	if(!param->outdir){
		fprintf(stderr,"You need to specify an output file suffix with -o option...\n");
		exit(EXIT_FAILURE);
	}
	/* create output structure ...  */
	RUN(create_output_directories(param->outdir));
	
        RUN(set_log_file(param->outdir,"scs_net"));
	
	RUN(print_program_header(argv,"TESTUGN a"));
	
	RUN(log_command_line(argc, argv));

	
	RUN(run_scs(param));
	int t_count_size = 10000;
	
	int max_transcript_frequency = 0;
	double* fx_unique_transcript_count = 0;
	double* fitted = 0;
	
	
	struct gigp_param* gigp_param = 0;
	MMALLOC(gigp_param,sizeof(struct gigp_param));
	gigp_param->gamma = 0.0;
	gigp_param->b = 0.0;
	gigp_param->c = 0.0;
	gigp_param->N = 0.0;
	gigp_param->s = 0.0;
	
	gigp_param->S = 0.0;
	gigp_param->fit = 0.0;
	gigp_param->max_count = 0.0;

	MFREE(param);

	return EXIT_SUCCESS;
	MMALLOC(outfile_name,sizeof(char) *(strlen(outfile) + 50 ));
        
	sprintf (outfile_name, "%s_GIGP_param.csv",outfile);
	if(my_file_exists(outfile_name)) {
		if (!(file = fopen(outfile_name , "r" ))){
			fprintf(stderr,"Cannot open output file '%s'\n",outfile_name);
			exit(EXIT_FAILURE );
		}
		
		//sscanf (sentence,"%s %*s %d",str,&i);
		//printf ("%s -> %d\n",str,i);
		fscanf(file,"%lg\tgamma\n",  &gigp_param->gamma );
		fscanf(file,"%lg\tb\n",  &gigp_param->b );
		fscanf(file,"%lg\tc\n",  &gigp_param->c );
		fscanf(file,"%lf\tN\n",  &gigp_param->N );
		fscanf(file,"%lf\ts\n",  &gigp_param->s);
		fscanf(file,"%lf\tS\n",  &gigp_param->S );
		fscanf(file,"%lf\tMax_transcript_count\n",  &gigp_param->max_count);
		fscanf(file,"%lf\tFit\n",  &gigp_param->fit);
		
				
		fclose(file);
		gigp_param->max_count++;
		
		fprintf(stderr,"%10.10e	%10.10e	%10.10e	N:%f	s:%f	S:%f	%f\n",gigp_param->gamma,gigp_param->b,gigp_param->c,gigp_param->N,gigp_param->s,gigp_param->S ,gigp_param->fit);
		
	}else{
		
		
		fx_unique_transcript_count =  read_count_distribution(argc,argv ,&max_transcript_frequency,data_origin, sample_names, num_samples);
		
		fprintf(stderr,"READ:%d\n",max_transcript_frequency );
		//exit(0);
		
		
		
		
		
		
		
		gigp_param->gamma = -2.0;
		gigp_param->b = 1.0;
		gigp_param->c = 0.1;
		gigp_param->N = 0.0;
		gigp_param->s = 0.0;
		
		gigp_param->S = 0.0;
		gigp_param->fit = DBL_MAX;
		gigp_param->max_count = max_transcript_frequency;
		
		for(i = 0; i < max_transcript_frequency;i++){
			gigp_param->N += i * fx_unique_transcript_count[i];
			gigp_param->s += fx_unique_transcript_count[i];
		}
		
		//gigp_param = fit_gigp_controller(gigp_param,fx_unique_transcript_count,max_transcript_frequency  );
		
		
		
		MMALLOC(fitted,sizeof(double) *max_transcript_frequency );
		fitted =  fill_fitted_curve(fitted , max_transcript_frequency,gigp_param );
		fprintf(stderr,"%f	%f	%f	N:%f	s:%f	S:%f	%f\n",gigp_param->gamma,gigp_param->b,gigp_param->c,gigp_param->N,gigp_param->s,gigp_param->S ,gigp_param->fit);
		
		
		sprintf (outfile_name, "%s_GIGP_param.csv",outfile);
		
		if (!(file = fopen(outfile_name , "w" ))){
			fprintf(stderr,"Cannot open output file '%s'\n",outfile_name);
			exit(EXIT_FAILURE );
		}
		fprintf(file,"%10.10e\tgamma\n",  gigp_param->gamma );
		fprintf(file,"%10.10e\tb\n",  gigp_param->b );
		fprintf(file,"%10.10e\tc\n",  gigp_param->c );
		fprintf(file,"%f\tN\n",  gigp_param->N );
		fprintf(file,"%f\ts\n",  gigp_param->s);
		fprintf(file,"%f\tS\n",  gigp_param->S );
		fprintf(file,"%f\tMax_transcript_count\n",  gigp_param->max_count -1 );
		fprintf(file,"%f\tFit\n",  gigp_param->fit);
		
		fclose(file);
		
		
		sprintf (outfile_name, "%s_GIGP_fit.csv",outfile);
		if (!(file = fopen(outfile_name , "w" ))){
			fprintf(stderr,"Cannot open output file '%s'\n",outfile_name);
			exit(EXIT_FAILURE);
		}
		
		for(i = 1; i < max_transcript_frequency;i++){
			if(fx_unique_transcript_count[i]){
				fprintf(file,"%d\t%f\t%f\n",  i,fx_unique_transcript_count[i],fitted[i]);
			}
		}
		
		fclose(file);
		
		sprintf (outfile_name, "%s_GIGP_growth_curve.csv",outfile);
		if (!(file = fopen(outfile_name , "w" ))){
			fprintf(stderr,"Cannot open output file '%s'\n",outfile_name);
			exit(EXIT_FAILURE);
		}
		
		sum = gigp_param->N;
		
		for(i = 1000000; i <= 50000000;i+= 1000000){
			if(i < sum && (double)(i + 1000000.0) > sum){
				gigp_param->N = sum;
				fprintf(file,"%d\t%d\t%d\n",(int)sum, (int)(gigp_param->S * ( 1.0 - give_me_sichel_p0(gigp_param,fx_unique_transcript_count))),(int) gigp_param->S);
			}
			gigp_param->N = i;
			fprintf(file,"%d\t%d\t%d\n",i, (int)(gigp_param->S * ( 1.0 - give_me_sichel_p0(gigp_param,fx_unique_transcript_count))),(int) gigp_param->S);
		}
		
		gigp_param->N = sum;
		fclose(file);
		MFREE(fx_unique_transcript_count);
		
	}
	
	
	
	

	//exit(EXIT_SUCCESS);
	////
	//	Step two: select fake relative frequencies from distribution and simulate.
	/////
	
	//gigp_param->S = 20;
	
	
	/*
		
	MMALLOC(rel_abundance, sizeof(double) * (int) gigp_param->S);
	MMALLOC(observed, sizeof(double) * (int) gigp_param->S);
	
	
	
	sprintf (outfile_name, "%s_GIGP_rel_abundances.csv",outfile);
	if( access( outfile_name, F_OK ) != -1 ) {
		if (!(file = fopen(outfile_name , "r" ))){
			fprintf(stderr,"Cannot open output file '%s'\n",outfile_name);
			exit(EXIT_FAILURE );
		}
		for(i = 0; i <(int)gigp_param->S;i++){
			fscanf(file,"%lg\n",  &rel_abundance[i]  );
		}
		
		fclose(file);
	}else{
		for(i = 0; i <(int)gigp_param->S;i++){
			rel_abundance[i] = 0.0;
			observed[i] = (double)rand_r(&seed)/(double)RAND_MAX;
		}
		
		qsort(observed, (int)gigp_param->S, sizeof(double), compare_double);
		
		rel_abundance = pick_abundances(gigp_param, observed, rel_abundance, 0 ,(int)gigp_param->S,1e-7,1);
		
		qsort(rel_abundance, (int)gigp_param->S, sizeof(double), compare_double);
		
		if (!(file = fopen(outfile_name , "w" ))){
			fprintf(stderr,"Cannot open output file '%s'\n",outfile_name);
			exit(EXIT_FAILURE );
		}
		for(i = 0; i <(int)gigp_param->S;i++){
			fprintf(file,"%10.10e\n",  rel_abundance[i]  );
		}
		fclose(file);
	}
	//sampling... 
	
	
	make_fake_expression_table(CV,rel_abundance, outfile,number_of_cells,simulated_read_depth,(int)gigp_param->S,loss);
	*/
        
	return EXIT_SUCCESS;
ERROR:
	
//	free_float_matrix(m);
	MFREE(param);
	return EXIT_FAILURE;
}


int run_scs(struct parameters* param)
{
	struct shared_data* bsd = NULL;
	struct float_matrix* m = NULL;
	struct unique_transcript_count* utc = NULL;

	float* uniq_transcript_vector = NULL;	
	char buffer[BUFFER_LEN];	
	/* allocate shared data  */
	RUNP(bsd = init_shared_data(param));
	
	/* create checkpoint... */
	snprintf(buffer, BUFFER_LEN, "%s/%s/",param->outdir,OUTDIR_CHECKPOINTS);
	DECLARE_CHK(MAIN_CHECK, buffer);
        
	/* Read in counts..  */
	RUNP(m = read_float_matrix(param->infile,1,1));

	/* Turn counts into frequency vector. */
	RUNP(bsd->utc = get_unique_transcript_vector(m,1));

	/* set inital guesss */
	RUN(set_initial_guess(bsd));
	
	/* do modelling...  */
	snprintf(buffer, BUFFER_LEN, "GAGA");
	
	RUN_CHECKPOINT(MAIN_CHECK,fit_model(bsd->utc),buffer);
	DESTROY_CHK(MAIN_CHECK);
	
	free_float_matrix(m);
	free_shared_data(bsd);
	return OK;
ERROR:
	
	free_shared_data(bsd);
	free_float_matrix(m);
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


struct unique_transcript_count* get_unique_transcript_vector(struct float_matrix* m, int column)
{     
	struct unique_transcript_count* utc = NULL;
	int i,c; 
	ASSERT(m != NULL,"Count matrix is NULL.");

	ASSERT(column < m->nrow,"Count matrix has %d columns but you ask for %d.", m->nrow,column );

	MMALLOC(utc, sizeof(struct unique_transcript_count));
	utc->x = NULL;
	utc->alloc_size = 10000;
	utc->max = 0;
	

	MMALLOC(utc->x, sizeof(double)* utc->alloc_size);
	for(i = 0; i < utc->alloc_size;i++){
		utc->x[i] = 0.0;
	}
	for(i = 0; i < m->nrow;i++){
		if(m->matrix[i][column] != 0){
			while( (int)m->matrix[i][column]  >= utc->alloc_size){
				utc->x = realloc(utc->x, sizeof(double)*  utc->alloc_size *2);
				for(c =  utc->alloc_size;c <  utc->alloc_size*2 ;  c++ ){
					utc->x[c] = 0.0;
				}
				utc->alloc_size =  utc->alloc_size << 1;
			}
			if(utc->max < (int) m->matrix[i][column] ){
				utc->max = (int) m->matrix[i][column];
			}
			utc->x[(int) m->matrix[i][column]] += 1.0;
		}
	}
	return utc;
ERROR:
	return NULL;
}

int fit_model(struct unique_transcript_count* utc)
{

	struct gigp_param* gigp_param = NULL; 

	ASSERT(utc != NULL,"No counts....");

	MMALLOC(gigp_param,sizeof(struct gigp_param));
	gigp_param->gamma = 0.0;
	gigp_param->b = 0.0;
	gigp_param->c = 0.0;
	gigp_param->N = 0.0;
	gigp_param->s = 0.0;
	
	gigp_param->S = 0.0;
	gigp_param->fit = 0.0;
	gigp_param->max_count = 0.0;

	RUN(fit_gigp_controller(utc, gigp_param));//, utc->x, utc->max+1);
		
	return OK;
ERROR:
	return FAIL;
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
	bsd->gigp_param_best  = NULL;
	bsd->gigp_param_initial_guess  = NULL;
	
	bsd->working_gigp_param = NULL;
	bsd->num_threads = param->num_threads;
	bsd->available_work_space = NULL;
	MMALLOC(bsd->available_work_space, sizeof(uint8_t) *bsd->num_threads);
	

	pthread_mutex_init(&bsd->avail_mtx,NULL);

	bsd->pool = NULL;
	if((bsd->pool = thr_pool_create(param->num_threads+1, param->num_threads+1, 0, 0)) == NULL) ERROR_MSG("Creating pool thread failed.");
	
	RUNP(tmp = init_gigp_param());
	bsd->gigp_param_best = tmp;
	tmp = NULL;

	RUNP(tmp = init_gigp_param());
	bsd->gigp_param_initial_guess = tmp;
	tmp = NULL;
	
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
		if(bsd->utc){
			MFREE(bsd->utc->x);
			MFREE(bsd->utc);
		}
		pthread_mutex_destroy(&bsd->avail_mtx);

		if(bsd->available_work_space){
			MFREE(bsd->available_work_space);
		}

		if(bsd->pool){
			thr_pool_destroy(bsd->pool);
		}
		MFREE(bsd);
	}
}

struct gigp_param* init_gigp_param(void)
{
	struct gigp_param* tmp = NULL; 
	MMALLOC(tmp,sizeof(struct gigp_param));
	tmp->gamma = 0.0;
	tmp->b = 0.0;
	tmp->c = 0.0;
	tmp->N = 0.0;
	tmp->s = 0.0;
	tmp->S = 0.0;
	tmp->fit = 0.0;
	tmp->max_count = 0.0;
	return tmp;
ERROR:
	return NULL;
}

int copy_gigp_param(struct gigp_param* source, struct gigp_param* target)
{
	ASSERT(source != NULL," No source.");
	ASSERT(target != NULL," No target.");
	/* copy estimated parameters back  */
	target->gamma = source->gamma;
	target->b = source->b;
	target->c = source->c;
	target->N = source->N;
	target->s = source->s;	
	target->S = source->S;
	target->fit = source->fit;
	target->max_count = source->max_count;
	
	return OK;
ERROR:
	return FAIL;
}


int  make_fake_expression_table(double CV,double* rel_abundance,char* outfile, int ncells, int depth,int num_genes,double loss)
{
	FILE* file = 0;
	int i,j;
	gsl_rng_env_setup();
	const gsl_rng_type *T = gsl_rng_default;
	gsl_rng *rfff = gsl_rng_alloc (T);
	
	double* cell_abundance = 0;
	MMALLOC(cell_abundance, sizeof(double) * num_genes);
	
	char* outfile_name = 0;
	
	MMALLOC(outfile_name, sizeof(char) *(strlen(outfile) + 50 ));
	
	
	int** table = 0;
	MMALLOC(table, sizeof(int*) * ncells);
	for(i = 0; i < ncells;i++){
		table[i] = 0;
		MMALLOC(table[i], sizeof(int) * num_genes );
		for(j = 0; j < num_genes;j++){
			table[i][j] = 0;
		}
	}
	
	
	for(i = 0; i < ncells;i++){
		
		for(j = 0; j < num_genes;j++){
			cell_abundance[j] = rel_abundance[j] + gsl_ran_gaussian (rfff,  rel_abundance[j] * CV);
			if(cell_abundance[j] < 0){
				cell_abundance[j] = 0.0;
			//	fprintf(stderr,"%f	%f\n",rel_abundance[j],cell_abundance[j]);
			}
		}
		gsl_ran_discrete_t *sample_abundances = gsl_ran_discrete_preproc (num_genes, cell_abundance);

		
		for(j = 0;j < depth;j++){
			table[i][(int)gsl_ran_discrete(rfff, sample_abundances)] += 1;
		}
		
		gsl_ran_discrete_free(sample_abundances);
	}
	
	//random deletion...
	
	for(i = 0; i < ncells;i++){
		
		for(j = 0; j < num_genes;j++){
			if(gsl_rng_uniform_pos(rfff) < loss){
				table[i][j] = 0;
			}
		}
	}
	
	
	sprintf (outfile_name, "%s_extable_%d_%d_CV%f_loss%f.csv",outfile,ncells,depth,CV,loss);
	
	if (!(file = fopen(outfile_name , "w" ))){
		fprintf(stderr,"Cannot open output file '%s'\n",outfile_name);
		exit(EXIT_FAILURE );
	}
	
	///table[i] = malloc(sizeof(int) * num_genes );
	for(j = 0; j < num_genes;j++){
		for(i = 0; i < ncells;i++){
			fprintf(file,"%d\t",table[i][j] );
		}
		fprintf(file,"\n");
	}

	fclose(file);
	for(i = 0; i < ncells;i++){
		MFREE(table[i]);// = malloc(sizeof(int) * num_genes );
	}
	gsl_rng_free(rfff);
	MFREE(table);
	MFREE(cell_abundance);
	MFREE(outfile_name);
	return OK;
ERROR:
	return FAIL;
}


double* read_count_distribution(int argc, char * argv[] ,int* max_transcript_frequency,int data_origin,char** sample_names,int num_samples)
{
	struct tome_data* td = 0;
	FILE* file = 0;
	int numseq,i,c,gzipped = 0;
	char* column_spec = 0;
	
	
	char* db_file = 0;
	char* input_file = 0;
	
	/*if(num_samples > 1){
		fprintf(stderr,"Warning: I will only use the first sample... \n");
	}
	
	if(!outfile){
		fprintf(stderr,"You need to specify an output file suffix with -o option...\n");
		exit(EXIT_FAILURE);
	}*/
	
	if(data_origin == 1){
		
		db_file = argv[optind++];
		if(!db_file){
			fprintf(stderr,"No database file...\n");
			exit(EXIT_FAILURE);
		}
		
		
		input_file = argv[optind++];
		if(!input_file){
			fprintf(stderr,"No input file...\n");
			exit(EXIT_FAILURE);
		}
		
	}else if(data_origin == 3){
		db_file = argv[optind++];
		if(!db_file){
			fprintf(stderr,"No database file...\n");
			exit(EXIT_FAILURE);
		}
		input_file = argv[optind++];
		if(!input_file){
			fprintf(stderr,"No input file...\n");
			exit(EXIT_FAILURE);
		}
	}else{
		
		fprintf(stderr,"No input method specified.\n");
		exit(EXIT_SUCCESS);
	}
	
	
	
	
	double* fx_unique_transcript_count = 0;
	
	*max_transcript_frequency = 0;
	
	int t_count_size = 10000;
	MMALLOC(fx_unique_transcript_count, sizeof(double)* t_count_size);
	for(i= 0; i < t_count_size;i++){
		fx_unique_transcript_count[i] = 0.0;
	}
	//
	//sample = XXXX; which sample to select
	//
	if(input_file){
		/*	file = open_track(input_file, &gzipped);
		column_spec = get_column_specs_based_on_file_suffix(input_file);
		while ((numseq = read_data_chunk(db,td,column_spec,file))) {
			td->numseq = numseq;
			for(i = 0; i < td->numseq;i++){
				if(td->gc[i]->count != 0.0){
					while( (int) td->gc[i]->count >= t_count_size){
							
						fx_unique_transcript_count = realloc(fx_unique_transcript_count, sizeof(double)* t_count_size *2);
						for(c = t_count_size;c < t_count_size*2 ;  c++ ){
							fx_unique_transcript_count[c] = 0.0;
						}
						t_count_size = t_count_size *2;
					}
					if(*max_transcript_frequency < (int) td->gc[i]->count ){
						*max_transcript_frequency = (int) td->gc[i]->count ;
					}
					fx_unique_transcript_count[ (int) td->gc[i]->count ] += 1.0;
				}
			}
			}*/
		
	}
	//close_track(file, gzipped);
	//free_tome_data(td);
	
	////
	//	Step one: estimation on GIGP parameters using Nelder & Mead simplex method....
	/////
	
	
	*max_transcript_frequency = *max_transcript_frequency + 1;
	return fx_unique_transcript_count;
ERROR:
	return NULL;
}


int fit_gigp_controller(struct unique_transcript_count* utc, struct gigp_param* gigp_param)
{
	
	struct gigp_param* gigp_param_best = NULL;
	int i;	

	/* set initial guess...  */
	gigp_param->gamma = -2.0;
	gigp_param->b = 1.0;
	gigp_param->c = 0.1;
	gigp_param->N = 0.0;
	gigp_param->s = 0.0;
		
	gigp_param->S = 0.0;
	gigp_param->fit = DBL_MAX;
	gigp_param->max_count = utc->max;
		
	for(i = 0; i < utc->max;i++){
		gigp_param->N += i * utc->x[i];
		gigp_param->s +=  utc->x[i];
	}

	/* allov best fit  */
	MMALLOC(gigp_param_best, sizeof(struct gigp_param));
	gigp_param_best->gamma = -2.0;
	gigp_param_best->b = 1.0;
	gigp_param_best->c = 0.1;
	gigp_param_best->N = gigp_param->N;
	gigp_param_best->s = gigp_param->s;
	gigp_param_best->S = gigp_param->S;
	gigp_param_best->fit = DBL_MAX;
	gigp_param_best->max_count  = gigp_param->max_count ;
	/* fit model  */
	RUN(fit_gigp(gigp_param,utc->x, utc->max));
	
	if(gigp_param->fit < gigp_param_best->fit ){
		fprintf(stderr,"%f	%f	%f	N:%f	s:%f	S:%f	%f\n",gigp_param->gamma,gigp_param->b,gigp_param->c,gigp_param->N,gigp_param->s,gigp_param->S ,gigp_param->fit);
		gigp_param_best->fit = gigp_param->fit;
		gigp_param_best->gamma = gigp_param->gamma;
		gigp_param_best->b = gigp_param->b;
		gigp_param_best->c =  gigp_param->c;
	}
	/* try a couple of times more with differeny starting points...  */
	for(i = 1; i < 20;i++){
		gigp_param->gamma = 1 - (double) i * 0.2;
		gigp_param->b = 0.1 + (double) i * 0.01;
		gigp_param->c = 0.01 +  (double) i * 0.001;
		
		RUN(fit_gigp(gigp_param,utc->x, utc->max));
		fprintf(stderr,"i:%d	%f	%f	%f	N:%f	s:%f	S:%f	%f\n",i,gigp_param->gamma,gigp_param->b,gigp_param->c,gigp_param->N,gigp_param->s,gigp_param->S ,gigp_param->fit);

		if(gigp_param->fit < gigp_param_best->fit ){
			fprintf(stderr,"%f	%f	%f	N:%f	s:%f	S:%f	%f\n",gigp_param->gamma,gigp_param->b,gigp_param->c,gigp_param->N,gigp_param->s,gigp_param->S ,gigp_param->fit);
			gigp_param_best->fit = gigp_param->fit;
			gigp_param_best->gamma = gigp_param->gamma;
			gigp_param_best->b = gigp_param->b;
			gigp_param_best->c =  gigp_param->c;
		}
	}
	/* copy estimated parameters back  */
	gigp_param->gamma = gigp_param_best->gamma;
	gigp_param->b = gigp_param_best->b;
	gigp_param->c = gigp_param_best->c;
	gigp_param->N = gigp_param_best->N;
	gigp_param->s = gigp_param_best->s;
	
	gigp_param->S = gigp_param_best->S;
	gigp_param->fit = gigp_param_best->fit;
	gigp_param->max_count = gigp_param_best->max_count;

	/* free best struct..  */
	MFREE(gigp_param_best);

	return OK;
ERROR:
	if(gigp_param_best){
	  	MFREE(gigp_param_best);
	}
	return FAIL;
}




int fit_gigp(struct gigp_param* gigp_param,double* fx_unique_transcript_count, int len)
{
	//double par[5] = {1.0, 2.0, 10.0, 20.0, 30.0};
	
	const gsl_multimin_fminimizer_type *T =
	gsl_multimin_fminimizer_nmsimplex2rand;
	gsl_multimin_fminimizer *s = NULL;
	gsl_vector *ss, *x;
	gsl_multimin_function minex_func;
	
	size_t iter = 0;
	int status;
	double size;
	
	
	double* par = 0;
	MMALLOC(par,sizeof(double) * (len + 3));
	
	par[0] =gigp_param->N;
	par[1] =gigp_param->s;
	par[2] = gigp_param->max_count;
	for(status = 0;status < len;status++){
		par[3+status] = fx_unique_transcript_count[status];
	}
	/*for(status = 0; status < 10;status++){
		fprintf(stderr,"%d\t%f	%f\n",status,par[status],fx_unique_transcript_count[status]);
	}*/
	
	/* Starting point */
	//x = gsl_vector_alloc (2);
	//gsl_vector_set (x, 0, 5.0);
	//gsl_vector_set (x, 1, 7.0);
	

	x = gsl_vector_alloc (3);
	
	gsl_vector_set (x, 0,  gigp_param->gamma );
	gsl_vector_set (x, 1, gigp_param->b);
	gsl_vector_set (x, 2, gigp_param->c);
	//gsl_vector_set (x, 3, gigp_param->N);
	//gsl_vector_set (x, 4, gigp_param->s);
	//gsl_vector_set (x, 5, gigp_param->max_count);
	
	

	
	/* Set initial step sizes to 1 */
	ss = gsl_vector_alloc (3);
	gsl_vector_set (ss, 0,gigp_param->gamma /1.5);
	gsl_vector_set (ss, 1, gigp_param->b/1.5);
	gsl_vector_set (ss, 2, gigp_param->c/1.5);

	
	//gsl_vector_set_all (ss, 1e-4);
	
	/* Initialize method and iterate */
	minex_func.n = 3;
	minex_func.f = sichel_function;
	minex_func.params = par;
	
	s = gsl_multimin_fminimizer_alloc (T, 3);
	gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
	
	do
	{
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);
		
		if (status)
			break;
		
		size = gsl_multimin_fminimizer_size (s);
		status = gsl_multimin_test_size (size, 1e-9);
		
		/*if (status == GSL_SUCCESS)
		{
			fprintf (stderr,"converged to minimum at\n");
		}
		
		fprintf (stderr,"%d\t%10.3e\t%10.3e\t%10.3e\t\tf() = %7.3f size = %.3f\n",
		        iter,
		        gsl_vector_get (s->x, 0),
		        gsl_vector_get (s->x, 1),
		         gsl_vector_get (s->x, 2),
		        s->fval, size);*/
	}
	while (status == GSL_CONTINUE && iter < 100);
	
	gigp_param->gamma = gsl_vector_get (s->x, 0);
	gigp_param->b = gsl_vector_get (s->x, 1);
	gigp_param->c = gsl_vector_get (s->x, 2);
	gigp_param->fit = s->fval;
	
	

	
	gsl_vector_free(x);
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free (s);
        
	return OK;
ERROR:
	return FAIL;
}


double my_f (const gsl_vector *v, void *params)
{
	double x, y;
	double *p = (double *)params;
	
	x = gsl_vector_get(v, 0);
	y = gsl_vector_get(v, 1);
	
	return p[2] * (x - p[0]) * (x - p[0]) + p[3] * (y - p[1]) * (y - p[1]) + p[4];
}


double sichel_function_for_fly(void*data,int n, double* var, int num_var)
{
	double *p = (double *)data;
	double gamma,b,c,N,s,max_observed, log_likelihood ;
	gamma = var[0];
	b =var[1];
	c = var[2];
	N =var[3];
	s = var[4];// p[1]; //gsl_vector_get(v, 4);
	max_observed = var[5];// p[2];//sl_vector_get(v, 5);
	
	//p += 2;
	
	gsl_sf_result result;
	
	log_likelihood = 0.0;
	int status = 0;
	
	double res[3];
	
	double sum = 0.0;
	
	double alpha = b*sqrt(1.0 + c *N );
	double beta = (b*c*N)/(2.0 * sqrt(1.0 + c *N ));
	
	double bessel1 = 0.0;
	double bessel2 = 0.0;
	
	double x = 0;
	double y = 0;
	status = gsl_sf_bessel_Kn_e( gamma, sqrt( pow(alpha,2.0) - 2.0* alpha*beta )   ,  &result);
	
	if (status == GSL_SUCCESS) {
		//printf("%f	%f\n",i,result.val);
		bessel1 =  result.val;
	}
	
	status = gsl_sf_bessel_Kn_e( gamma + x,  alpha   ,  &result);
	
	if (status == GSL_SUCCESS) {
		//printf("%f	%f\n",i,result.val);
		bessel2 =  result.val;
	}
	
	
	
	
	y = pow((  sqrt( pow(alpha,2.0) - 2.0* alpha*beta  )   / alpha  ),gamma) * 1.0 / bessel1 * (pow(beta,x)/  1.0 )  * bessel2;
	sum+= y;
	//log_likelihood = log_likelihood - s * log(1.0 -y);
	//fprintf(stderr," , )
	//fprintf(stdout,"%f	%f	%f\n", x,y,y * 77152.0);
	res[0] = y;
	
	x = 1.0 ;
	status = gsl_sf_bessel_Kn_e( gamma, sqrt( pow(alpha,2.0) - 2.0* alpha*beta )   ,  &result);
	if (status == GSL_SUCCESS) {
		bessel1 =  result.val;
	}
	
	status = gsl_sf_bessel_Kn_e( gamma + x,  alpha   ,  &result);
	if (status == GSL_SUCCESS) {
		bessel2 =  result.val;
	}
	
	
	
	
	y = pow ((  sqrt( pow(alpha,2.0) - 2.0* alpha*beta  )   / alpha  ),gamma) * 1.0 / bessel1 * (pow(beta,x)/  1.0 )  * bessel2;
	if(p[(int)x]){
		sum+= y;
	//	log_likelihood += p[(int)x] *log(y);
	}
	
	//fprintf(stdout,"%f	%f	%f\n", x,y,y * 77152.0);
	res[1] = y;
	
	for(x = 2.0; x < max_observed;x+=1.0){
		y = (2.0*beta)/ alpha * ((x + gamma -1.0)/ x) * res[1] + (pow(beta,2.0)/ ( x *(x-1.0) )  )* res[0];
	//	fprintf(stdout,"%f	%f	%f\n", x,y,y * 77152.0);
		if(p[(int)x]){
			sum+= y;
			//log_likelihood += p[(int)x] *log(y);
		}
		res[0] = res[1];
		res[1] = y;
		
		
	}
	
	
	//fprintf(stderr,"SUM:%f\n",sum);
	x = 0;
	status = gsl_sf_bessel_Kn_e( gamma, sqrt( pow(alpha,2.0) - 2.0* alpha*beta )   ,  &result);
	
	if (status == GSL_SUCCESS) {
		//printf("%f	%f\n",i,result.val);
		bessel1 =  result.val;
	}
	
	status = gsl_sf_bessel_Kn_e( gamma + x,  alpha   ,  &result);
	
	if (status == GSL_SUCCESS) {
		//printf("%f	%f\n",i,result.val);
		bessel2 =  result.val;
	}
	
	
	
	
	y = pow((  sqrt( pow(alpha,2.0) - 2.0* alpha*beta  )   / alpha  ),gamma) * 1.0 / bessel1 * (pow(beta,x)/  1.0 )  * bessel2;
	//sum+= y;
	log_likelihood = log_likelihood - s * log(1.0 -(y/sum));
	//fprintf(stderr," , )
	//fprintf(stdout,"%f	%f	%f\n", x,y,y * 77152.0);
	res[0] = y;
	
	x = 1.0 ;
	status = gsl_sf_bessel_Kn_e( gamma, sqrt( pow(alpha,2.0) - 2.0* alpha*beta )   ,  &result);
	if (status == GSL_SUCCESS) {
		bessel1 =  result.val;
	}
	
	status = gsl_sf_bessel_Kn_e( gamma + x,  alpha   ,  &result);
	if (status == GSL_SUCCESS) {
		bessel2 =  result.val;
	}
	
	
	
	
	y = pow ((  sqrt( pow(alpha,2.0) - 2.0* alpha*beta  )   / alpha  ),gamma) * 1.0 / bessel1 * (pow(beta,x)/  1.0 )  * bessel2;
	if(p[(int)x]){
		//sum+= y;
		log_likelihood += p[(int)x] *log(y/sum);
	}
	
	//fprintf(stdout,"%f	%f	%f\n", x,y,y * 77152.0);
	res[1] = y;
	
	for(x = 2.0; x < max_observed;x+=1.0){
		y = (2.0*beta)/ alpha * ((x + gamma -1.0)/ x) * res[1] + (pow(beta,2.0)/ ( x *(x-1.0) )  )* res[0];
		//	fprintf(stdout,"%f	%f	%f\n", x,y,y * 77152.0);
		if(p[(int)x]){
			//sum+= y;
			log_likelihood += p[(int)x] *log(y/sum);
		}
		res[0] = res[1];
		res[1] = y;
		
		
	}

	
	//p-=2;
	if(isnan(log_likelihood)){
		return -100000000000000.0;
	}
	
	return   log_likelihood;
}

double give_me_sichel_p0(struct gigp_param* gigp_param,double* p )
{
	
	double gamma = gigp_param->gamma;
	double b = gigp_param->b;
	double c = gigp_param->c;
	double N = gigp_param->N;
	int max_observed = gigp_param->max_count;
	int status;
	
	gsl_sf_result result;
	
	double alpha = b*sqrt(1.0 + c *N );
	double beta = (b*c*N)/(2.0 * sqrt(1.0 + c *N ));
	
	double sum = 0.0;
		
	double res[3];
	double bessel1 = 0.0;
	double bessel2 = 0.0;
	
	double x = 0;
	double y = 0;
	
	status = gsl_sf_bessel_Kn_e( gamma, sqrt( pow(alpha,2.0) - 2.0* alpha*beta )   ,  &result);
	
	if (status == GSL_SUCCESS) {
		//printf("%f	%f\n",i,result.val);
		bessel1 =  result.val;
	}
	
	status = gsl_sf_bessel_Kn_e( gamma + x,  alpha   ,  &result);
	
	if (status == GSL_SUCCESS) {
		//printf("%f	%f\n",i,result.val);
		bessel2 =  result.val;
	}
	
	
	
	
	y = pow((  sqrt( pow(alpha,2.0) - 2.0* alpha*beta  )   / alpha  ),gamma) * 1.0 / bessel1 * (pow(beta,x)/  1.0 )  * bessel2;
	sum += y;
	
	res[0] = y;
	
	x = 1.0 ;
	status = gsl_sf_bessel_Kn_e( gamma, sqrt( pow(alpha,2.0) - 2.0* alpha*beta )   ,  &result);
	if (status == GSL_SUCCESS) {
		bessel1 =  result.val;
	}
	
	status = gsl_sf_bessel_Kn_e( gamma + x,  alpha   ,  &result);
	if (status == GSL_SUCCESS) {
		bessel2 =  result.val;
	}
	
	
	
	
	y = pow ((  sqrt( pow(alpha,2.0) - 2.0* alpha*beta  )   / alpha  ),gamma) * 1.0 / bessel1 * (pow(beta,x)/  1.0 )  * bessel2;
	//if(p[(int)x]){
		sum += y;
	//}
	
	res[1] = y;
	
	for(x = 2.0; x < max_observed;x+=1.0){
		y = (2.0*beta)/ alpha * ((x + gamma -1.0)/ x) * res[1] + (pow(beta,2.0)/ ( x *(x-1.0) )  )* res[0];
		//if(p[(int)x]){
			sum += y;
		//}
		
		res[0] = res[1];
		res[1] = y;
		
		
	}
	x = 0;
	
	status = gsl_sf_bessel_Kn_e( gamma, sqrt( pow(alpha,2.0) - 2.0* alpha*beta )   ,  &result);
	
	if (status == GSL_SUCCESS) {
		bessel1 =  result.val;
	}
	
	status = gsl_sf_bessel_Kn_e( gamma + x,  alpha   ,  &result);
	
	if (status == GSL_SUCCESS) {
		bessel2 =  result.val;
	}
	
	
	
	
	y = pow((  sqrt( pow(alpha,2.0) - 2.0* alpha*beta  )   / alpha  ),gamma) * 1.0 / bessel1 * (pow(beta,x)/  1.0 )  * bessel2;
	
	return y/sum;
}


double sichel_function(const gsl_vector *v, void *params)
{
	
	double *p = (double *)params;
	double gamma,b,c,N,s,max_observed, log_likelihood ;
	gamma = gsl_vector_get(v, 0);
	b = gsl_vector_get(v, 1);
	c = gsl_vector_get(v, 2);
	N = p[0];// gsl_vector_get(v, 3);
	s = p[1]; //gsl_vector_get(v, 4);
	max_observed = p[2];//sl_vector_get(v, 5);
	
	
	double sum = 0.0;
	
	gsl_sf_result result;

	log_likelihood = 0.0;
	int status = 0;
	
	double res[3];
	if(b < 0.0 ){
		return DBL_MAX;
		b = DBL_MIN;
	}
	
	if(c < 0.0){
		return DBL_MAX;
		c = DBL_MIN;
	}
	
	if(c > b){
		return DBL_MAX;
	}
	
	
	double alpha = b*sqrt(1.0 + c *N );
	double beta = (b*c*N)/(2.0 * sqrt(1.0 + c *N ));
	
		
	double bessel1 = 0;
	double bessel2 = 0;
	
	double x = 0;
	double y = 0;
	
	status = gsl_sf_bessel_Kn_e( gamma, sqrt( pow(alpha,2.0) - 2.0* alpha*beta )   ,  &result);
	
	if (status == GSL_SUCCESS) {
		//printf("%f	%f\n",i,result.val);
		bessel1 =  result.val;
	}else{
		ERROR_MSG("Blah");
	}
	
	status = gsl_sf_bessel_Kn_e( gamma + x,  alpha   ,  &result);
	
	if (status == GSL_SUCCESS) {
		//printf("%f	%f\n",i,result.val);
		bessel2 =  result.val;
	}else{
		ERROR_MSG("Blah");
	}
	
	
	
	
	y = pow((  sqrt( pow(alpha,2.0) - 2.0* alpha*beta  )   / alpha  ),gamma) * 1.0 / bessel1 * (pow(beta,x)/  1.0 )  * bessel2;
	sum += y;
	
	res[0] = y;
	
	x = 1.0 ;
	status = gsl_sf_bessel_Kn_e( gamma, sqrt( pow(alpha,2.0) - 2.0* alpha*beta )   ,  &result);
	if (status == GSL_SUCCESS) {
		bessel1 =  result.val;
	}else{
		ERROR_MSG("Blah");
	}
	
	status = gsl_sf_bessel_Kn_e( gamma + x,  alpha   ,  &result);
	if (status == GSL_SUCCESS) {
		bessel2 =  result.val;
	}else{
		ERROR_MSG("Blah");
	}
	
	
	
	
	y = pow ((  sqrt( pow(alpha,2.0) - 2.0* alpha*beta  )   / alpha  ),gamma) * 1.0 / bessel1 * (pow(beta,x)/  1.0 )  * bessel2;
	//if(p[(int)x+3]){
		sum += y;
	//}
	
	res[1] = y;
	
	for(x = 2.0; x < max_observed;x+=1.0){
		y = (2.0*beta)/ alpha * ((x + gamma -1.0)/ x) * res[1] + (pow(beta,2.0)/ ( x *(x-1.0) )  )* res[0];
		//if(p[(int)x+3]){
			sum += y;
		//}
		
		res[0] = res[1];
		res[1] = y;
		
		
	}
	x = 0;
	
	status = gsl_sf_bessel_Kn_e( gamma, sqrt( pow(alpha,2.0) - 2.0* alpha*beta )   ,  &result);
	
	if (status == GSL_SUCCESS) {
		bessel1 =  result.val;
	}else{
		ERROR_MSG("Blah");
	}
	status = gsl_sf_bessel_Kn_e( gamma + x,  alpha   ,  &result);
	
	if (status == GSL_SUCCESS) {
		bessel2 =  result.val;
	}else{
		ERROR_MSG("Blah");
	}
	
	
	
	
	y = pow((  sqrt( pow(alpha,2.0) - 2.0* alpha*beta  )   / alpha  ),gamma) * 1.0 / bessel1 * (pow(beta,x)/  1.0 )  * bessel2;
	log_likelihood = log_likelihood - s * log(1.0 -y/sum);
	res[0] = y;
	
	x = 1.0 ;
	status = gsl_sf_bessel_Kn_e( gamma, sqrt( pow(alpha,2.0) - 2.0* alpha*beta )   ,  &result);
	if (status == GSL_SUCCESS) {
		bessel1 =  result.val;
	}else{
		ERROR_MSG("Blah");
	}
	
	status = gsl_sf_bessel_Kn_e( gamma + x,  alpha   ,  &result);
	if (status == GSL_SUCCESS) {
		bessel2 =  result.val;
	}else{
		ERROR_MSG("Blah");
	}
	
	
	
	y = pow ((  sqrt( pow(alpha,2.0) - 2.0* alpha*beta  )   / alpha  ),gamma) * 1.0 / bessel1 * (pow(beta,x)/  1.0 )  * bessel2;
	if(p[(int)x+3]){
		log_likelihood += p[(int)x+3] *log(y/sum);
	}
	res[1] = y;
	
	for(x = 2.0; x < max_observed;x+=1.0){
		y = (2.0*beta)/ alpha * ((x + gamma -1.0)/ x) * res[1] + (pow(beta,2.0)/ ( x *(x-1.0) )  )* res[0];
		if(p[(int)x+3]){
			log_likelihood += p[(int)x+3] *log(y/sum);
		}
		
		res[0] = res[1];
		res[1] = y;
	}
	
	
	//fprintf (stderr,"%10.3e\t%10.3e\t%10.3e	%f\n",gamma,alpha,beta,log_likelihood);

	if(isnan(log_likelihood)){
		return DBL_MAX;
	}
	
	return -1.0 * log_likelihood;
ERROR:
	return -1;
}



double* fill_fitted_curve(double* fit,int num,struct gigp_param* param )
{
	double alpha = param->b*sqrt(1.0 + param->c *param->N );
	double beta = (param->b*param->c*param->N)/(2.0 * sqrt(1.0 + param->c *param->N ));

	double gamma = param->gamma;
	double res[3];
	
	
	double y = 0;
	gsl_sf_result result;
	
	double x = 0;
	double sum = 0.0;
	double bessel1 = 0.0;
	double bessel2 = 0.0;
	int status = 0;

	x = 0;
	
	
	status = gsl_sf_bessel_Kn_e( gamma, sqrt( pow(alpha,2.0) - 2.0* alpha*beta )   ,  &result);
	
	if (status == GSL_SUCCESS) {
		//printf("%f	%f\n",i,result.val);
		bessel1 =  result.val;
	}
	
	status = gsl_sf_bessel_Kn_e( gamma + x,  alpha   ,  &result);
	
	if (status == GSL_SUCCESS) {
		//printf("%f	%f\n",i,result.val);
		bessel2 =  result.val;
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
	}
	
	status = gsl_sf_bessel_Kn_e( gamma + x,  alpha   ,  &result);
	if (status == GSL_SUCCESS) {
		bessel2 =  result.val;
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
	
	return fit;
}


double gigp_dist (double x, void * p)
{
	struct gigp_param * params = (struct gigp_param *)p;
	double gamma = (params->gamma);
	double b = (params->b);
	double c = (params->c);
	
	int status = 0;
	double bessel1 = 0.0;
	double y = 0.0;
	gsl_sf_result result;
	
	status = gsl_sf_bessel_Kn_e( gamma, b    ,  &result);
	
	if (status == GSL_SUCCESS) {
		//printf("%f	%f\n",i,result.val);
		bessel1 =  result.val;
	}
	y = (  pow(2.0/(b*c),gamma) /(2.0 * bessel1)  ) * pow(x,gamma-1.0) * exp(  -1.0*(x/c) -((b*b*c)/ (4.0*x)) );
	
	
	
	return y;
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
	gsl_integration_workspace * w	= gsl_integration_workspace_alloc (100);
	struct gigp_param param_4integration;
	param_4integration.gamma =  gigp_param->gamma;
	param_4integration.b = gigp_param->b;
	param_4integration.c = gigp_param->c;
	
	F.function = &gigp_dist;
	F.params = &param_4integration;
	double x,r;
	double integration_result, error,sum;
	//r =   (double)rand_r(&seed)/(double)RAND_MAX;
	x = 1;
	gsl_integration_qag (&F, 1e-7, x, 0, 1e-7, 100,1,    w, &integration_result, &error);
	sum = integration_result;
	
	
	
	double inner_mid = 0.0;
	int outer_mid;
	
	
	
	outer_mid =  (outer_high + outer_low) / 2;
	
	double tmp[3];
	
	tmp[0] = inner_low;
	tmp[1] = inner_mid;
	tmp[2] = inner_high;
	
	
	//fprintf(stderr,"Working on:%d	%d	%d	%10.10e	%10.10e\n", (int) outer_mid,   outer_low,outer_high,inner_low,inner_high  );
	
	
	r = random[(int)outer_mid];
	if(r == -1.0){
		gsl_integration_workspace_free (w);
		return abundances;
	}
	while(1){
		inner_mid = (inner_high + inner_low )/ 2.0;
		gsl_integration_qag (&F, 1e-7, inner_mid, 0, 1e-7, 100,1,    w, &integration_result, &error);
		integration_result = integration_result / sum;
		//fprintf(stderr,"L:%e\tH:%e\tmid:%e\t\t%f\t%f\n", inner_low,inner_high,inner_mid, integration_result, r  );
		
		if(fabs(integration_result - r) < DBL_EPSILON){
			abundances[outer_mid] = inner_mid;
			random[outer_mid] = -1.0;
			break;
		}else if(integration_result < r){
			inner_low = inner_mid;
		}else	 if(integration_result > r ){
			inner_high = inner_mid;
		}
		
		if(fabs(inner_high - inner_low)  < DBL_EPSILON){
			abundances[outer_mid] = inner_mid;
			random[outer_mid] = -1.0;
			break;
		}
	}
	
	gsl_integration_workspace_free (w);
	
	if((outer_mid-1) -outer_low > 1000 ){
		fprintf(stderr,"Working on : %d - %d\n",outer_low, outer_mid-1);
	}
	abundances = pick_abundances(gigp_param,random,abundances, outer_low, outer_mid-1,tmp[0],inner_mid);
		
	if(outer_high - (outer_mid+1)  > 1000){
		fprintf(stderr,"Working on : %d - %d\n",outer_mid+1, outer_high);
	}
	abundances = pick_abundances(gigp_param,random,abundances, outer_mid+1, outer_high,  inner_mid,tmp[2]);
	
	return  abundances;
}






