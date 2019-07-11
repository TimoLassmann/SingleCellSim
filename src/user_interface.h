#ifndef USER_INTERFACE

#define USER_INTERFACE


struct parameters{
        char** sample_names;
        char** gene_names;
        
        char* infile;
        char* outdir;

        double CV;
        double drop_Michaelisconstant;

        int simulated_read_depth;
        double simulated_read_depth_SD;
        int num_samples;
        int num_gene_targets; 

        int fitonly;
        int num_threads;
        int num_genes;
        int num_cells;
};

extern int print_global_help(int argc, char * argv[]);

extern struct parameters* get_model_param(int argc, char * argv[]);
extern struct parameters* get_sim_param(int argc, char * argv[]);

extern void free_param(struct parameters* param);


#endif // USER_INTERFACE

