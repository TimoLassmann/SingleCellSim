#ifndef gigp_header

#define gigp_header


#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>



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

#endif
