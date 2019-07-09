#ifndef MATRIX_IO_H
#define MATRIX_IO_H


#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "tldevel.h"

#include <string.h>
#include "ctype.h"



struct double_matrix{
        double** matrix;
        char** col_names;
        char** row_names;
        void* matrix_mem;
        uint8_t* label;
        int nrow;
        int ncol;
        int real_sample;

};

extern struct double_matrix* read_dm(char* filename, int has_col_names,int has_row_names);
extern struct double_matrix* read_double_matrix(char* filename, int has_col_names,int has_row_names);
extern struct double_matrix* alloc_double_matrix(int ncol,int nrow, int name_len);
extern void free_double_matrix(struct double_matrix* m);

extern int print_double_matrix(struct double_matrix* m,FILE* file, int has_col_names,int has_row_names);

extern int fill_random_matrix(struct double_matrix* m);
extern int shuffle_double_matrix(struct double_matrix* m);

//extern int analyze_data_file(char* filename);

extern struct double_matrix* transpose_double_matrix(struct double_matrix* m);

#endif
