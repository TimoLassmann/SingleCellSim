#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "tldevel.h"
#include "gigp.h"




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

	par[0] = gigp_param->N;
	par[1] = gigp_param->s;
	par[2] = gigp_param->max_count;
	for(status = 0;status < len;status++){ /*  */
		par[3+status] = fx_unique_transcript_count[status];
	}
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
	MFREE(par);

	return OK;
ERROR:
	return FAIL;
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



struct gigp_param* init_gigp_param(void)
{
        struct gigp_param* tmp = NULL;
        MMALLOC(tmp,sizeof(struct gigp_param));

        RUN(clear_gigp_param(tmp));
        return tmp;
ERROR:
        return NULL;
}


int clear_gigp_param(struct gigp_param* tmp)
{
        ASSERT(tmp != NULL,"No parameters..");
        tmp->gamma = 0.0;
        tmp->b = 0.0;
        tmp->c = 0.0;
        tmp->N = 0.0;
        tmp->s = 0.0;
        tmp->S = 0.0;
        tmp->fit = DBL_MAX;
        tmp->max_count = 0.0;
        return OK;
ERROR:
        return FAIL;
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
