#ifndef gigp_header

#define gigp_header


#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define GIGP_GAMMA_ROW 0
#define GIGP_B_ROW 1
#define GIGP_C_ROW 2
#define GIGP_N_ROW 3
#define GIGP_s_ROW 4
#define GIGP_S_ROW 5
#define GIGP_MAXCOUNT_ROW 6
#define GIGP_FIT_ROW 7




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


extern struct gigp_param* fit_exact_gigp(struct gigp_param* gigp_param,double* fx_unique_transcript_count, int len);
extern double gigp_dist (double x, void * p);
extern double sichel_function(const gsl_vector *v, void *params);
extern struct gigp_param* init_gigp_param(void);
extern int copy_gigp_param(struct gigp_param* source, struct gigp_param* target);
extern double give_me_sichel_p0(struct gigp_param* gigp_param ,double* p);


/* housekeeping stuff */
extern struct gigp_param* init_gigp_param(void);
extern int clear_gigp_param(struct gigp_param* tmp);
extern int copy_gigp_param(struct gigp_param* source, struct gigp_param* target);

#endif
