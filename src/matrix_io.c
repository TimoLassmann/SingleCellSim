#include "matrix_io.h"

/* Goal: I want to have a program that reads in a matrix of numbers and recognises the presence of column / rownames */

#include "hash_table.h"

#define DELIMCOMMA ','
#define DELIMSPACE ' '
#define DELIMTAB '\t'

#define DELIMCOMMAINDEX 0
#define DELIMSPACEINDEX 1
#define DELIMTABINDEX 2

struct matrix_file_specs{
        char delim[2];
        int has_row_names;
        int has_col_names;
        int ncol;
        int nrow;
        int max_name_len;
        int max_line_len;
};

HT_GLOBAL_INIT(HTI, int);
struct hash_table_node_HTI_t;



int analyze_token(char*t,int len, int* res);

struct row_col_poperties{
        double entry_len;
        double is_alpha;
        double is_cat;
        int num;
};
struct row_col_poperties* alloc_rc_prop(int num);

int detect_outlier_rc(struct row_col_poperties* p, int len, double threshold, int* ret);


static int fill_double_matrix(struct double_matrix* m,char* filename, struct matrix_file_specs* spec);
static struct matrix_file_specs*  analyze_data_file(char* filename);

struct double_matrix* transpose_double_matrix(struct double_matrix* m)
{
        struct double_matrix* t = NULL;
        int i,j;
        ASSERT(m != NULL,"No input matrix.");

        /* get dim of new matrix */
        RUNP(t = alloc_double_matrix(m->nrow, m->ncol, BUFFER_LEN));

        /* copy names */
        for(i = 0; i < m->ncol;i++){
                snprintf(t->row_names[i],BUFFER_LEN,"%s", m->col_names[i]);
        }

        for(i = 0; i < m->nrow;i++){
                snprintf(t->col_names[i],BUFFER_LEN,"%s", m->row_names[i]);
        }
        /* copy values */
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

struct double_matrix* read_dm(char* filename, int has_col_names,int has_row_names)
{
        struct double_matrix* m = NULL;

        struct matrix_file_specs* spec = NULL;


        RUNP(spec = analyze_data_file(filename));


        LOG_MSG("File: %s", filename);
        LOG_MSG("%s\tdelimiter",spec->delim);
        LOG_MSG("%d\trownames",spec->has_row_names);
        LOG_MSG("%d\tcolnames",spec->has_col_names);
        LOG_MSG("%d\tncol",spec->ncol);
        LOG_MSG("%d\tnrow",spec->nrow);


        if(has_col_names != -1){

                if(spec->has_col_names == 1 && has_col_names == 0){
                        WARNING_MSG("Detected column names in file %s contradicting command line argument!",filename);
                        spec->nrow++;
                }else if(spec->has_col_names == 0 && has_col_names == 1){
                        WARNING_MSG("Detected no column names in file %s contradicting command line argument!",filename);
                        spec->nrow--;
                }


                spec->has_col_names = has_col_names;
        }

        if(has_row_names != -1){

                if(spec->has_row_names == 1 && has_row_names == 0){
                        WARNING_MSG("Detected row names in file %s contradicting command line argument!",filename);
                        spec->ncol++;
                }else if(spec->has_row_names == 0 && has_row_names == 1){
                        WARNING_MSG("Detected no row names in file %s contradicting command line argument!",filename);
                        spec->ncol--;
                }

                spec->has_row_names = has_row_names;
        }
        m = alloc_double_matrix(spec->ncol,spec->nrow,spec->max_name_len+1);


        RUN(fill_double_matrix( m,filename, spec));
        MFREE(spec);
        return m;
ERROR:
        if(spec){
                MFREE(spec);
        }
        return NULL;
}

int fill_double_matrix(struct double_matrix* m,char* filename, struct matrix_file_specs* spec)
{
        FILE* f_ptr = NULL;

        char* line_buffer = NULL;
        char *token;
        int c;
        int col_index;
        int row_index;


        line_buffer = galloc(line_buffer,spec->max_line_len);

        RUNP(f_ptr = fopen(filename, "r"));
        /* First row */
        if(spec->has_col_names){
                fgets(line_buffer,spec->max_line_len,f_ptr);
                line_buffer[strcspn(line_buffer, "\r\n")] = 0;
                token = strtok(line_buffer, spec->delim);
                col_index = 0;
                while(token != NULL){
                        c = strnlen(token, spec->max_line_len);
                        if(col_index == 0 && spec->has_row_names){
                        }else{

                                strncpy(m->col_names[col_index - spec->has_row_names],token,c+1);
                                //fprintf(stdout,"%d %s  %s\n",col_index- spec->has_row_names,token,m->col_names[col_index - spec->has_row_names]);
                        }
                        token = strtok(NULL, spec->delim);
                        col_index++;
                }
        }
        row_index = 0;
        while(fgets(line_buffer,spec->max_line_len,f_ptr)){
                line_buffer[strcspn(line_buffer, "\r\n")] = 0;
                token = strtok(line_buffer, spec->delim);
                col_index = 0;
                while(token != NULL){
                        c = strnlen(token, spec->max_line_len);
                        if(col_index == 0 && spec->has_row_names){
                                //fprintf(stdout,"%d %s\n",row_index,token);
                                strncpy(m->row_names[row_index],token, c+1);
                        }else if(col_index !=0  && spec->has_row_names){
                                m->matrix[row_index][col_index-1] = atof(token);
                        }else{
                                m->matrix[row_index][col_index] = atof(token);
                        }
                        token = strtok(NULL, spec->delim);
                        col_index++;
                }
                row_index++;
        }
        //exit(0);
        fclose(f_ptr);
        gfree(line_buffer);
        return OK;
ERROR:
        gfree(line_buffer);
        return FAIL;
}

struct matrix_file_specs*  analyze_data_file(char* filename)
{
        HT_TYPE(HTI)* ht = NULL;

        FILE* f_ptr = NULL;
        char* line_buffer = NULL;
        char delim[2];
        int** delimiter_count = NULL;
        int sum_delim[3];

        struct matrix_file_specs* spec = NULL;




        int i,j;
        int l;
        int c;
        int key;
        int line_number = 1;
        int line_length = 0;
        int max_line_length = -1;
        int has_col_names = -1;
        int has_row_names = -1;
        int n_columns;
        int max_entry_len;
        int col_number;
        ASSERT(filename != NULL, "No file name.");

        ASSERT(my_file_exists(filename),"Could not find file %s.",filename);


        ht = HT_INIT(HTI,65536);
        LOG_MSG("reading from %s",filename);
        RUNP(f_ptr = fopen(filename, "r" ));
        while((c = fgetc(f_ptr))){

                if(c == EOF){
                        break;
                }
                if(c == '\n'){
                        line_length++;
                        if(max_line_length < line_length){
                                max_line_length = line_length;

                        }
                        line_length = 0;
                        line_number++;
                        //if(line_number ==256){
                        //        break;
                        //}
                }else{
                        line_length++;
                        //fprintf(stdout,"%c %d %d\n", (char)c, line_number,(c << 8) | line_number);

                        if(line_number < 256){
                                key = (c << 8) | line_number;

                                RUN(HT_INSERT(HTI,ht,key,NULL));
                                key = (c << 8) | 0;
                                RUN(HT_INSERT(HTI,ht,key,NULL));
                        }
                }

        }

        //LOG_MSG("Found %d lines.",line_number);
        LOG_MSG("Max line length is %d characters.",max_line_length);
        HT_FLATTEN(HTI,ht);

        /*for(i = 0; i < ht->num_items;i++){
                fprintf(stdout,"HASH:%d %d  '%c' %d %d\n",i,ht->flat[i]->key, (char) ((ht->flat[i]->key >> 8) &0xFF), ht->flat[i]->key&0xFF, ht->flat[i]->count);
                }*/

        delimiter_count = galloc(delimiter_count,3,line_number-1,0);
        /* test for delimiter  */
        for(i = 0; i < ht->num_items;i++){
                c = ((ht->flat[i]->key >> 8) &0xFF);

                l = ht->flat[i]->key & 0xFF;
                if(l){
                        l = l-1;
                        if(c == DELIMCOMMA){
                                delimiter_count[DELIMCOMMAINDEX][l] = ht->flat[i]->count;
                        }else if(c == DELIMSPACE){
                                delimiter_count[DELIMSPACEINDEX][l] = ht->flat[i]->count;
                        }else if(c == DELIMTAB){
                                delimiter_count[DELIMTABINDEX][l] = ht->flat[i]->count;
                        }
                }
        }
        /* Free hash table */
        HT_FREE(HTI, ht);
        /* We are only looking at the top 256 lines to work out delimiters */
        line_number = MACRO_MIN(line_number-1,255);
        for(i = 0; i < 3;i++){
                sum_delim[i] = 0;
        }
        for(i = 0; i < 3;i++){
                for(l = 0; l < line_number-1;l++){
                        sum_delim[i] += delimiter_count[i][l];
                }
        }
        for(i = 0; i < 3;i++){
                if(sum_delim[i]){
                        for(j = 0; j < line_number-1;j++){
                                for(c = j+1;c <line_number-1;c++){
                                        if( delimiter_count[i][j] != delimiter_count[i][c]){
                                                sum_delim[i] = 0;
                                        }
                                }
                        }
                }
        }
        n_columns = -1;

        c = -1;
        for(i = 0; i < 3;i++){
                if(sum_delim[i]){
                        n_columns =  delimiter_count[i][0] +1;
                        switch (i) {
                        case DELIMCOMMAINDEX: {
                                c = DELIMCOMMA;
                                break;
                        }
                        case DELIMSPACEINDEX: {
                                c = DELIMSPACE;
                                break;
                        }
                        case DELIMTABINDEX: {
                                c = DELIMTAB;
                                break;
                        }
                        }

                }
        }
        if(c == -1){
                ERROR_MSG("No delimiter was found.");
        }
        if(n_columns == -1){
                ERROR_MSG("Could not determine number of columns.");
        }

        LOG_MSG("MaxLineLength = %d", max_line_length);
        max_line_length++;



        struct row_col_poperties* col_properties = NULL;
        struct row_col_poperties* row_properties = NULL;

        col_properties  = alloc_rc_prop(n_columns);
        row_properties  = alloc_rc_prop(256);



        delim[0] = c;
        delim[1] = 0;
        line_buffer = galloc(line_buffer,max_line_length);
        line_number = 0;
        max_entry_len = -1;
        rewind(f_ptr);
        while(fgets(line_buffer,max_line_length,f_ptr)){
                if(line_number < 5){
                fprintf(stdout,"LINE:%d %.10s\n",line_number,line_buffer);
                }
                line_buffer[strcspn(line_buffer, "\r\n")] = 0;

                if(line_number < 256){
                        col_number = 0;

                        char *token = strtok(line_buffer, delim);
                        while(token != NULL) {
                                c = strnlen(token, max_line_length);
                                if(c > max_entry_len ){
                                        max_entry_len = c;
                                }
                                analyze_token(token,c,&l);
                                //fprintf(stdout,"%s(%d,%d) ",token,c,l);
                                token = strtok(NULL, delim);
                                col_properties[col_number].entry_len += (double)c;
                                col_properties[col_number].is_alpha += l;
                                col_properties[col_number].is_cat += 0.0;
                                col_properties[col_number].num += 1;

                                row_properties[line_number].entry_len += (double)c;
                                row_properties[line_number].is_alpha += (double)l;
                                row_properties[line_number].is_cat += 0.0;
                                row_properties[line_number].num  +=1;
                                col_number++;
                        }

                        //fprintf(stdout,"newline...\n");
                }
                line_number++;

        }


        /* Now I should have column and row properties */
        RUN(detect_outlier_rc(col_properties, n_columns, 1.0f, &has_row_names));
        RUN(detect_outlier_rc(row_properties, MACRO_MIN(line_number,256), 1.0f, &has_col_names));

        fclose(f_ptr);

        MMALLOC(spec, sizeof(struct matrix_file_specs));
        spec->delim[0] = delim[0];
        spec->delim[1] = delim[1];
        spec->has_col_names = has_col_names;
        spec->has_row_names = has_row_names;
        spec->max_name_len = max_entry_len;
        spec->ncol = n_columns - has_row_names;
        spec->nrow = line_number-has_col_names;
        spec->max_line_len = max_line_length;



        MFREE(col_properties);
        MFREE(row_properties);
        gfree(delimiter_count);
        gfree(line_buffer);

        return spec;
ERROR:
        gfree(delimiter_count);
        gfree(line_buffer);
        if(f_ptr){
                fclose(f_ptr);
        }
        return NULL;
}

int detect_outlier_rc(struct row_col_poperties* p, int len, double threshold, int* ret)
{
        double s1[3];
        double s2[3];
        int i;
        for(i = 0; i < 3;i++){
                s1[i] = 0.0;
                s2[i] = 0.0;
        }
        /* Now I should have column and row properties */
        for(i = 0; i < len;i++){
                p[i].entry_len /= (double) p[i].num;
                p[i].is_alpha /= (double) p[i].num;
                p[i].is_cat /= (double) p[i].num;
                p[i].num = 0;
                if(i){          /* Look at first row.. */
                        s1[0] += p[i].entry_len ;

                        s2[0] += p[i].entry_len * p[i].entry_len;

                        s1[1] += p[i].is_alpha ;
                        s2[1] += p[i].is_alpha * p[i].is_alpha;

                        s1[2] += p[i].is_cat ;
                        s2[2] += p[i].is_cat * p[i].is_cat;
                }

        }

        for(i = 0; i < 3;i++){
                s2[i] = sqrt(((double) len * s2[i] - s1[i] * s1[i])/ ((double)len * ((double)len -1.0)));
                s1[i] = s1[i] / (double) len;
        }
        //for(i = 0; i < len;i++){
        *ret  =0;
        i = 0;
                if(s2[0]){
                        if(fabs(p[i].entry_len-s1[0]) /  s2[0] >= threshold){
                                //fprintf(stdout,"%d %f %f %f  TEST:%f %f \n",i,p[i].entry_len,s1[0],s2[0],fabs(p[i].entry_len-s1[0]) /  s2[0],fabs(p[i].entry_len-s1[0]));
                                p[i].num++;
                        }
                }
                if(s2[1]){
                        if(fabs(p[i].is_alpha -s1[1]) /  s2[1] >= threshold){

                                p[i].num++;
                        }
                }
                if(s2[2]){
                        if(fabs(p[i].is_cat -s1[2]) /  s2[2] >= threshold){

                                p[i].num++;
                        }
                }
                if(p[i].num){
                        LOG_MSG("%d might be an outlier",i);
                        *ret = 1;
                }
                //}
        return OK;
}

struct row_col_poperties* alloc_rc_prop(int num)
{
        struct row_col_poperties* p = NULL;
        int i;

        MMALLOC(p, sizeof(struct row_col_poperties) * num);

        for(i = 0; i < num;i++){
                p[i].entry_len = 0.0;
                p[i].is_alpha = 0.0;
                p[i].is_cat = 0.0;
                p[i].num = 0;
        }
        return p;
ERROR:
        return NULL;
}

int analyze_token(char*t,int len, int* res)
{
        int i;
        int alpha;
        int num;
        alpha = 0;
        num = 0;
        for(i = 0; i < len;i++){
                if(isalpha((int)t[i])){
                        alpha++;
                }
                if(isdigit((int)t[i])){
                        num++;
                }
        }
        if(num > alpha){
                *res = 0;
        }else{
                *res = 1;
        }

        return OK;
}

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

        longest_entry += 10; //in case I need to add stuff
        MMALLOC(tmp_storage,sizeof(char)*  (longest_entry+1));

        RUNP(m = alloc_double_matrix(max_columns,nrows,longest_entry+1));

        /*if(has_col_names){
          MFREE(m->col_names[m->ncol-1]);
          }
          if(has_row_names){
          MFREE(m->row_names[m->nrow-1]);
          }*/
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

                                                        m->matrix[i_mem][j_mem] = atof(tmp_storage);
                                                        j_mem++;
                                                }else{
                                                        for(lastchar = 0; lastchar <= n ;lastchar++){
                                                                m->row_names[i_mem][lastchar] = tmp_storage[lastchar];
                                                        }
                                                        //fprintf(stdout,"%d: %s\n",i_mem,m->row_names[i_mem]);
                                                }
                                        }else{
                                                m->matrix[i_mem][j_mem]  = atof(tmp_storage);
                                                j_mem++;
                                        }
                                }
                        }else{
                                if(has_row_names){
                                        if(j){
                                                m->matrix[i_mem][j_mem]  = atof(tmp_storage);
                                                j_mem++;
                                        }else{
                                                for(lastchar = 0; lastchar <= n ;lastchar++){
                                                        m->row_names[i_mem][lastchar] = tmp_storage[lastchar];
                                                }
                                        }
                                }else{
                                        m->matrix[i_mem][j_mem] = atof(tmp_storage);
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
                                                //fprintf(stdout,"%d %d\n", j_mem,lastchar);
                                                m->col_names[j_mem][lastchar] = tmp_storage[lastchar];
                                        }
                                }else{
                                        m->matrix[i_mem][j_mem] = atof(tmp_storage);
                                        j_mem++;
                                        i_mem++;
                                }
                        }else{
                                m->matrix[i_mem][j_mem]  = atof(tmp_storage);
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


        if(nrows > i_mem){
                for(i = i_mem;i < nrows;i++){
                        MFREE(m->matrix[i]);
                }
                MFREE(m->row_names[nrows-1]);
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
                fprintf(file,"%s", m->col_names[0]);
                for(j = 1; j < m->ncol;j++){
                        fprintf(file,",%s", m->col_names[j]);
                }
                fprintf(file,"\n");
        }
        /*if(has_row_names){
          fprintf(file,"  ");
          }
          for(j = 0; j < MACRO_MIN( m->ncol,1000);j++){
          fprintf(file,"%4u ", m->label[j]);
          }
          fprintf(file,"\n");*/
        for(i = 0; i < m->nrow;i++){
                if(has_row_names){
                        fprintf(file,"%s",m->row_names[i]);
                }
                for(j = 0; j < m->ncol;j++){
                        fprintf(file,",%10.10e", m->matrix[i][j]);
                        //fprintf(file,",%2.4f", m->matrix[i][j]);

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
