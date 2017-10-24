#include "scs.h"

struct double_matrix* read_double_matrix(char* filename, int has_col_names,int has_row_names)
{
        struct double_matrix* m = NULL;
        FILE* file = NULL;
        char* tmp_storage = NULL;
        int c,n,lastchar,i,j;
        int i_mem,j_mem;
        int32_t max_columns = INT32_MIN;
        int32_t min_columns = INT32_MAX;
        int cur_columns = 0;
        int nrows = 0;
        int longest_entry= 0;
        n = 0;

        RUNP(file = fopen(filename, "r" ));

        lastchar = 1;
        while((c = fgetc(file))){
                switch (c) {
                case EOF:
                        if(n > longest_entry){
                                longest_entry = n;
                        }
                        n = -1000;
                        break;
                case '\t':
                case ',':
                        ASSERT(lastchar != 0, "Multiple");
                        if(n > longest_entry){
                                longest_entry = n;
                        }
                        n = 0;
                        lastchar = 0;
                        cur_columns++;
                        break;
                case '\n':
                        ASSERT(lastchar != 0, "Multiple");
                        if(n > longest_entry){
                                longest_entry = n;
                        }
                        n = 0;
                        lastchar = 0;
                        cur_columns++;
                        if(cur_columns >max_columns){
                                max_columns =cur_columns;
                        }
                        if(cur_columns < min_columns){
                                min_columns = cur_columns;
                        }
                        nrows++;
                        cur_columns = 0;
                        break;
                default:
                        n++;
                        lastchar = 1;
                        break;
                }
                if(n == -1000){
                        break;
                }
        }

        ASSERT(max_columns == min_columns, "Rows seem to have different number of columns.%d %d",min_columns,max_columns);

        rewind(file);

        if(longest_entry < 50){
                longest_entry = 50;
        }
        
        MMALLOC(tmp_storage,sizeof(char)*  (longest_entry+1));

        RUNP(m = alloc_double_matrix(max_columns,nrows,longest_entry+1));

        if(has_col_names){
                MFREE(m->col_names[m->ncol-1]);
        }
        if(has_row_names){
                MFREE(m->row_names[m->nrow-1]);
        }
        m->ncol = -1;
        m->nrow = -1;
        i_mem = 0;
        j_mem = 0;
        i = 0;
        j = 0;
        n = 0;
        while((c = fgetc(file))){
                switch (c) {
                case EOF:
                        n = -1000;

                        break;
                case '\t':
                case ',':
                        tmp_storage[n] = 0;
                        if(!i){
                                if(has_col_names){
                                        if(has_row_names){
                                                if(j){
                                                        for(lastchar = 0; lastchar <= n ;lastchar++){
                                                                m->col_names[j_mem][lastchar] = tmp_storage[lastchar];

                                                        }
                                                        j_mem++;
                                                }
                                        }else{
                                                //fprintf(stderr,"First roow:%s %d	%d\n",tmp_storage,j,has_row_names);
                                                for(lastchar = 0; lastchar <= n ;lastchar++){
                                                        m->col_names[j_mem][lastchar] = tmp_storage[lastchar];
                                                }
                                                j_mem++;
                                        }
                                }else{
                                        if(has_row_names){
                                                if(j){

                                                        if(sscanf(tmp_storage,"%lf",&m->matrix[i_mem][j_mem]) != 1){
                                                                ERROR_MSG("could not read");
                                                        }
                                                        j_mem++;
                                                }else{
                                                        for(lastchar = 0; lastchar <= n ;lastchar++){
                                                                m->row_names[i_mem][lastchar] = tmp_storage[lastchar];
                                                        }
                                                        fprintf(stdout,"%d: %s\n",i_mem,m->row_names[i_mem]);
                                                }
                                        }else{
                                                if(sscanf(tmp_storage,"%lf",&m->matrix[i_mem][j_mem]) != 1){
                                                        ERROR_MSG("could not read");
                                                }
                                                j_mem++;
                                        }
                                }
                        }else{
                                if(has_row_names){
                                        if(j){
                                                if(sscanf(tmp_storage,"%lf",&m->matrix[i_mem][j_mem]) != 1){
                                                        ERROR_MSG("could not read");
                                                }
                                        
                                                j_mem++;
                                        }else{
                                                for(lastchar = 0; lastchar <= n ;lastchar++){
                                                        m->row_names[i_mem][lastchar] = tmp_storage[lastchar];
                                                }
                                        }
                                }else{
                                        if(sscanf(tmp_storage,"%lf",&m->matrix[i_mem][j_mem]) != 1){
                                                ERROR_MSG("could not read");
                                        }
                                        j_mem++;
                                }
                        }
                        j++;
                        n = 0;
                        break;
                case '\n':
                        tmp_storage[n] = 0;
                        if(i == 0){ // first row;
                                if(has_col_names){
                                        for(lastchar = 0; lastchar <= n ;lastchar++){
                                                m->col_names[j_mem][lastchar] = tmp_storage[lastchar];
                                        }
                                }else{
                                        if(sscanf(tmp_storage,"%lf",&m->matrix[i_mem][j_mem]) != 1){
                                                ERROR_MSG("could not read");
                                        }
                                        
                                        j_mem++;
                                        i_mem++;
                                }
                        }else{
                                if(sscanf(tmp_storage,"%lf",&m->matrix[i_mem][j_mem]) != 1){
                                        ERROR_MSG("could not read");
                                }
                                j_mem++;
                                i_mem++;
                        }
                        i++;
                        j = 0;
                        if(j_mem > m->ncol){
                                m->ncol = j_mem;
                        }
                        j_mem = 0;
                        n = 0;
                        break;
                default:
                        tmp_storage[n] = c;
                        n++;
                        break;
                }
                if(n == -1000){
                        break;
                }
        }

        //fprintf(stderr,"alloc %d real %d\n",nrows,i_mem);

        if(nrows > i_mem){
                for(i = i_mem;i < nrows;i++){
                        MFREE(m->matrix[i]);
                }
        }
        if(i_mem > m->nrow){
                m->nrow = i_mem;
        }
        n = 0;
        for(i = 0; i < m->ncol;i++){
                for(j = i+1;j < m->ncol;j++){
                        if(strcmp(m->col_names[i],m->col_names[j]) == 0){
                                //fprintf(stdout,"same:%s %s\n",m->col_names[i],m->col_names[j]);
                                n = 1;
                                i = m->real_sample;
                                j = m->real_sample;
                                break;
                        }
                }
        }

        if(n){
                for(i = 0; i < m->ncol;i++){
                        snprintf(tmp_storage,longest_entry,"%s_%d",m->col_names[i],i+1);
                        snprintf(m->col_names[i],longest_entry,"%s",tmp_storage);
                }
        }

        fclose(file);
        MFREE(tmp_storage);
        return m;
ERROR:
        return NULL;
}

struct double_matrix* alloc_double_matrix(int ncol,int nrow, int name_len)
{
        struct double_matrix* m = NULL;
        int i,j;
        //unsigned long int p;
        MMALLOC(m,sizeof(struct double_matrix));
        m->real_sample = ncol;
        m->ncol = ncol;
        m->nrow = nrow;
        m->name_len = name_len;
        m->matrix_mem = NULL;
        m->col_names = NULL;
        m->row_names = NULL;
        m->matrix = NULL;
        m->label = NULL;
        //MMALLOC(m->matrix_mem ,sizeof(double) *(int)(m->ncol + m->ncol%4)   *m->nrow  + 15);
        MMALLOC(m->col_names,sizeof(char*) * m->ncol);
        MMALLOC(m->label,sizeof(uint8_t) * m->ncol);
        MMALLOC(m->row_names,sizeof(char*) * m->nrow);
        MMALLOC(m->matrix ,sizeof(double*) * m->nrow);
        //p = (((unsigned long int) m->matrix_mem + 15) & (~0xf));
        for(i = 0; i <  m->nrow;i++){
                m->matrix[i] = NULL;
                MMALLOC(m->matrix[i], sizeof(double) * m->ncol);
                //m->matrix[i] = (double*) p;
                //p +=(unsigned long int)  (sizeof(double) * (int)(m->ncol + m->ncol%4));
                for(j = 0; j < m->ncol;j++){
                        m->matrix[i][j] = 0.0;
                }
        }

        for(i = 0; i < m->ncol;i++){
                if(i < m->real_sample){
                        m->label[i] = 1;
                }else{
                        m->label[i] = 0;
                }
                m->col_names[i]= NULL;
                MMALLOC(m->col_names[i],sizeof(char) * name_len);
                m->col_names[i][0] = 0;

        }

        for(i = 0; i < m->nrow ;i++){
                m->row_names[i]= NULL;
                MMALLOC(m->row_names[i], sizeof(char) * name_len);
                m->row_names[i][0]= 0;
        }
        return m;
ERROR:
        free_double_matrix(m);
        return NULL;
}

struct double_matrix* transpose_double_matrix(struct double_matrix* m)
{
        struct double_matrix* t = NULL;
        int i,j;
        ASSERT(m != NULL,"No input matrix.");
        
        /* get dim of new matrix */      
        RUNP(t = alloc_double_matrix(m->nrow, m->ncol,m->name_len));
        
        /* copy names  */
        for(i = 0; i < m->ncol;i++){
                snprintf(t->row_names[i],m->name_len,"%s", m->col_names[i]);
        }
        
        for(i = 0; i < m->nrow;i++){
                snprintf(t->col_names[i],m->name_len,"%s", m->row_names[i]);
        }
        /* copy values  */
        for(j = 0; j < m->nrow;j++){
                for(i = 0; i < m->ncol;i++){
                        t->matrix[i][j] = m->matrix[j][i];
                }
        }
        free_double_matrix(m);      
        return t;
ERROR:
        return NULL;
}

int add_rows_double_matrix(struct double_matrix* m, int n_extra)
{
        int i,j;
        int new_size;
        //unsigned long int p;

        new_size = m->nrow + n_extra;

        //MREALLOC(m->matrix_mem ,sizeof(double) *(int)(m->ncol + m->ncol%4)   * new_size  + 15);
        MREALLOC(m->row_names, sizeof(char* ) * new_size);
        MREALLOC(m->matrix   , sizeof(double*) * new_size);
        //p = (((unsigned long int) m->matrix_mem + 15) & (~0xf));
        //for(i = 0; i <  m->nrow;i++){
        //	p +=(unsigned long int)  (sizeof(double) * (int)(m->ncol + m->ncol%4));
        //}
        for(i = m->nrow; i <  new_size;i++){
                m->matrix[i] = NULL;
                MMALLOC(m->matrix[i], sizeof(double) * m->ncol);
                //m->matrix[i] = (double*) p;
                //p +=(unsigned long int)  (sizeof(double) * (int)(m->ncol + m->ncol%4));
                for(j = 0; j < m->ncol;j++){
                        m->matrix[i][j] = 0.0;
                }
        }

        j = 1;
        for(i = m->nrow; i < new_size ;i++){
                m->row_names[i]= NULL;
                MMALLOC(m->row_names[i], sizeof(char) * 16);
                snprintf(m->row_names[i],16,"RAND%d",j);
                j++;
        }
        m->nrow = new_size;
        return OK;
ERROR:
        return FAIL;

}

int remove_rows_double_matrix(struct double_matrix* m, int n_extra)
{
        int i;
        int new_size;
        new_size = m->nrow - n_extra;

        for(i = new_size; i < m->nrow;i++){
                MFREE(m->matrix[i]);//, sizeof(double) * m->ncol);
                MFREE(m->row_names[i]);
        }
        MREALLOC(m->row_names, sizeof(char* ) * new_size);
        MREALLOC(m->matrix   , sizeof(double*) * new_size);
        m->nrow = new_size;
        return OK;
ERROR:
        return FAIL;

}



int fill_random_matrix(struct double_matrix* m)
{
        int i,j,a;
        double tmp;
        int real = m->real_sample;
        for(i = 0; i < m->nrow;i++){
                /* COPY real values */
                for(j = 0;j < real;j++){
                        m->matrix[i][j+real] = m->matrix[i][j];

                }
                /* Shuffle values */
                for(j = m->real_sample-1; j > 0;j--){
                        a = random_int_zero_to_x(j) + real;
                        tmp = m->matrix[i][a];
                        m->matrix[i][a] = m->matrix[i][j+real];
                        m->matrix[i][j+real] = tmp;
                }
        }
        return OK;
}


int shuffle_double_matrix(struct double_matrix* m)
{
        int i,j,a;
        double tmp = 0.0;
        for(i = 0; i < m->nrow;i++){
                for(j = m->ncol-1;j > 0; j--){
                        a = random_int_zero_to_x(j);
                        tmp = m->matrix[i][a];
                        m->matrix[i][a] =  m->matrix[i][j];
                        m->matrix[i][j] = tmp;
                }
        }
        return OK;
}

int print_double_matrix(struct double_matrix* m,FILE* file, int has_col_names,int has_row_names)
{
        int i,j;
        if(has_col_names){
                if(has_row_names){
                        fprintf(file,",");
                }
                for(j = 0; j <  m->ncol;j++){
                        fprintf(file,"%s", m->col_names[j]);
                        if(j != m->ncol-1){
                                fprintf(file,",");
                        }

                }
                fprintf(file,"\n");
        }
        for(i = 0; i < m->nrow;i++){
                if(has_row_names){
                        fprintf(file,"%s,",m->row_names[i]);
                }
                for(j = 0; j < m->ncol;j++){
                        fprintf(file,"%10.10e", m->matrix[i][j]);
                        if(j != m->ncol-1){
                                fprintf(file,",");
                        }

                }
                fprintf(file,"\n");
        }
        return OK;
}

void free_double_matrix(struct double_matrix* m)
{
        int i;
        if(m){
                for(i = 0; i < m->ncol;i++){
                        MFREE(m->col_names[i]);// = malloc(sizeof(char) * (longest_entry+1));
                }
                for(i = 0; i < m->nrow ;i++){
                        MFREE(m->matrix[i]);

                        MFREE(m->row_names[i]);// = malloc(sizeof(char) * (longest_entry+1));
                }
                MFREE(m->label);
                if(m->matrix_mem){
                        MFREE(m->matrix_mem);/// = malloc( sizeof(double) *m->ncol  *m->nrow  + 15);
                }
                MFREE(m->col_names);/// = malloc(sizeof(char*) * m->ncol );
                MFREE(m->row_names);/// = malloc(sizeof(char*) * m->nrow);
                MFREE(m->matrix);/// = malloc(sizeof(double*) * m->nrow);
                MFREE(m);// = malloc(sizeof(struct sserdt_matrix));
        }
}
