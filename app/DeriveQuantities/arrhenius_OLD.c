#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_multifit_nlin.h>
#include"arrhenius.h"

typedef struct {
   size_t  n;
   double *t;
   double *y;
} fit_data;

#undef  DEBUG_LS    /*!< Debug mode turned off*/
#define NB_PAR 3    /*!< Number of modified Arrhenius law fitting parameters*/

/*!
 *  This function performs a three-parameter least square fit based on the modified Arrhenius law
 */
void fit_(size_t *ntmps, double *t, double *y, double *a, double *b, double *c, double *err) {

  size_t i;
  const gsl_multifit_fdfsolver_type *T;
  const gsl_rng_type * type;
  gsl_multifit_fdfsolver    *s;
  gsl_multifit_function_fdf  f;
  gsl_vector_view x;
  gsl_rng * r;
  int nb_temp;
  int status;
  unsigned int iter = 0;
  double chi, dof;
  double x_init[NB_PAR];
  fit_data d;

  gsl_matrix *covar = gsl_matrix_alloc (NB_PAR, NB_PAR);

  /* Data to be fitted are given in input */
  nb_temp = (*ntmps);
  d.n = nb_temp;
  d.t = t;
  d.y = y;

  /* Evaluate the logarithm of the rate coefficients */
  for (i = 0;i < nb_temp;i++) 
    d.y[i] = log(y[i]);

  /* Parameter guess values */
  x_init[0] = 1.0;
  x_init[1] = 10.0;
  x_init[2] = 1000.;
  x = gsl_vector_view_array(x_init, NB_PAR);

  gsl_rng_env_setup();

  type = gsl_rng_default;
  r = gsl_rng_alloc(type);

  f.f   = &arrhenius_f;
  f.df  = &arrhenius_df;
  f.fdf = &arrhenius_fdf;

  f.n = nb_temp;
  f.p = NB_PAR;
  f.params = &d;

  T = gsl_multifit_fdfsolver_lmsder;
  s = gsl_multifit_fdfsolver_alloc(T, nb_temp, NB_PAR);
  gsl_multifit_fdfsolver_set(s, &f, &x.vector);

#ifdef DEBUG_LS
  print_state (iter, s);
#endif

  do {
   iter++;
   status = gsl_multifit_fdfsolver_iterate(s);
#ifdef DEBUG_LS
   print_state (iter, s);
#endif
   if (status)
   break;
   status = gsl_multifit_test_delta (s->dx, s->x, 1e-12, 1e-12);
  } while (status == GSL_CONTINUE && iter < 1000);

  gsl_multifit_covar(s->J, 0.0, covar);

  /* Interolation coefficients */  
  #define FIT(i) gsl_vector_get(s->x, i)
  #define ERR(i) sqrt(gsl_matrix_get(covar,i,i))
  *a = exp(FIT(0));
  *b = FIT(1);
  *c = FIT(2);

  /* chi squared*/
  chi = gsl_blas_dnrm2(s->f);
  
  dof    = (double)(nb_temp) - (double)(NB_PAR);
  (*err) = GSL_MAX_DBL(1.,chi/sqrt(dof));

#ifdef DEBUG_LS
  {
    printf("chisq/dof = %g\n", pow(chi,2.)/dof);
    printf ("a = %.5f +/- %.5f\n", FIT(0), (*err)*ERR(0));
    printf ("b = %.5f +/- %.5f\n", FIT(1), (*err)*ERR(1));
    printf ("c = %.5f +/- %.5f\n", FIT(2), (*err)*ERR(2));
  }
#endif

  gsl_multifit_fdfsolver_free(s);

  gsl_matrix_free(covar);
  gsl_rng_free(r);
}

/*!
 *  This function provides the interpolation model (modified Arrhenius law)
 */
int arrhenius_f(const gsl_vector * x, void *data, gsl_vector * f) {

  size_t  i, nb_temp;
  double  a, b, c, Yi;
  double *y;
  double *t;

  /*Number of temperature*/
  nb_temp = ((fit_data *)data)->n;

  /* Rate data*/
  y = ((fit_data *)data)->y;
   
  /*Temperature vector*/
  t = ((fit_data *)data)->t;

  /*Fit parameters*/
  a = gsl_vector_get(x, 0);
  b = gsl_vector_get(x, 1);
  c = gsl_vector_get(x, 2);

  for (i = 0; i < nb_temp; i++) {
    Yi = a + b*log(t[i]) - c/t[i];
    gsl_vector_set(f, i, (Yi - y[i]));
  }

  return GSL_SUCCESS;
}

/*!
 *  This function provides the interpolation model Jacobian (modified Arrhenius law)
 */
int arrhenius_df(const gsl_vector * x, void *data, gsl_matrix * J) {

  size_t i, nb_temp;

  double *t;

  /*Number of temperature*/
  nb_temp = ((fit_data *)data)->n;

  /*Temperature vector*/
  t = ((fit_data *)data)->t;

  for (i = 0; i < nb_temp; i++) {
    gsl_matrix_set(J, i, 0, 1.);
    gsl_matrix_set(J, i, 1, log(t[i]));
    gsl_matrix_set(J, i, 2, - 1./t[i]);
  }

  return GSL_SUCCESS;
}

/*!
 *  This function provides the interpolation model function and the related Jacobian
 */
int arrhenius_fdf(const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J) {
  arrhenius_f(x, data, f);
  arrhenius_df(x, data, J);
  return GSL_SUCCESS;
}

/*!
 *  This function prints the state
 */
void print_state(size_t iter, gsl_multifit_fdfsolver * s) {

  printf ("iter: %lu x = % 15.8f % 15.8f % 15.8f " "|f(x)| = %g\n", iter, 
          gsl_vector_get (s->x, 0),
          gsl_vector_get (s->x, 1),
          gsl_vector_get (s->x, 2),
          gsl_blas_dnrm2 (s->f));
}

