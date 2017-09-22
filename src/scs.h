
#ifndef SCS_HEADER

#define SCS_HEADER 

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h> 
#include <stdio.h>
#include <getopt.h>
#include <stdint.h>
#include <float.h>

#include <stdio.h>
#include <math.h>
#include <string.h>



#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
 
#include "tldevel.h"

#include "thr_pool.h"

struct float_matrix{
	float** matrix; 
	char** col_names;
	char** row_names;	
        void* matrix_mem;
	uint8_t* label;
	int nrow;
	int ncol;
	int real_sample; 
};


struct float_matrix* read_float_matrix(char* filename, int has_col_names,int has_row_names);
struct float_matrix* alloc_float_matrix(int ncol,int nrow, int name_len);
int add_rows_float_matrix(struct float_matrix* m, int n_extra);
int remove_rows_float_matrix(struct float_matrix* m, int n_extra);
int fill_random_matrix(struct float_matrix* m);
int shuffle_float_matrix(struct float_matrix* m);
int print_float_matrix(struct float_matrix* m,FILE* file, int has_col_names,int has_row_names);
void free_float_matrix(struct float_matrix* m);

#endif
