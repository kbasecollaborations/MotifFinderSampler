#include "macros.h"
#include "log_erfc.h"
#include "io.h"

/*
	regress

	Least squares regression on points (x,y) to give
		y = mx + b

	Returns the root mean squared error of the fit.
*/
double regress(
  int n,			/* number of points */
  double *x,			/* x values */
  double *y,			/* y values */
  double *m,			/* slope */
  double *b 			/* y intercept */
)
{
  int i;
  double sx=0, sy=0, sxx=0, sxy=0;
  double mse=0;

  for (i=0; i<n; i++) {
    sx += x[i];
    sy += y[i];
    sxx += x[i]*x[i];
    sxy += x[i]*y[i];
  }

  double denom = n*sxx - sx*sx;
  if (denom != 0) {
    *m = (n*sxy - sy*sx) / denom;
    *b = (sy - *m*sx)/n;
  } else {
    *m = 0;
    *b = 0;
  }

  for (i=0; i<n; i++) {
    double err = y[i] - (*m*x[i] + *b);
    mse += err * err;
  }
  mse = sqrt(mse);
  mse /= n;

  return mse;
}
/*
	w_regress

	Weighted least squares regression on points (x,y,w) to give
		y = mx + b

	Returns the weighted root mean squared error of the fit.
*/
double w_regress(
  int n,			/* number of points */
  double *x,			/* x values */
  double *y,			/* y values */
  double *w,			/* weights */
  double *m,			/* slope */
  double *b 			/* y intercept */
)
{
  int i;
  double s=0, sx=0, sy=0, sxx=0, sxy=0;
  double mse=0;

  for (i=0; i<n; i++) {
    double ww = w[i];
    s += ww;
    sx += x[i]*ww;
    sy += y[i]*ww;
    sxx += x[i]*x[i]*ww;
    sxy += x[i]*y[i]*ww;
  }

  double denom = s*sxx - sx*sx;
  if (denom != 0) {
    *m = (s*sxy - sx*sy) / denom;
    *b = (sxx*sy - sx*sxy)/ denom;
  } else {
    *m = 0;
    *b = 0;
  }

  for (i=0; i<n; i++) {
    double err = (y[i] - (*m*x[i] + *b)) * w[i];
    mse += err * err;
  }
  mse = sqrt(mse);
  mse /= n;

  return mse;
}

/******************************************************************************
*	pearson_correlation
*	
*	Returns the sample Pearson correlation coefficient of two sets of points	
*	and computes an estimate of its significance (that it is large, one-tailed) 
*	using the Fisher transform.  Also computes the regression line
*	and its mean-squared error unless those variables are null.
*
******************************************************************************/
double pearson_correlation(
  int n,			/* number of points */
  double *x,			/* x values */
  double *y,			/* y values */
  double *m,			/* slope */
  double *b,			/* intercept */
  double *mse,			/* mean-squared error */
  double *log_pv 		/* log of Fisher transform p-value */
)
{
  int i;
  double sx=0, sy=0, sxy=0, sxx=0, syy=0; 
  double r=0;		// correlation
  double z=0;		// Fisher transform of r


  for (i=0; i<n; i++) {
    sx += x[i];
    sy += y[i];
    sxy += x[i]*y[i];
    sxx += x[i]*x[i];
    syy += y[i]*y[i];
  }

  // Get regression line.
  if (m != NULL && b != NULL) {
    double denom = n*sxx - sx*sx;
    if (denom != 0) {
      *m = (n*sxy - sx*sy) / denom;
      *b = (sy - *m*sx) / n;
    } else {
      *m = 0;
      *b = 0;
    }

    // Get mean-squared error.
    if (mse != NULL) {
      for (i=0; i<n; i++) {
	double err = y[i] - (*m*x[i] + *b);
	*mse += err * err;
      }
      *mse = sqrt(*mse);
      *mse /= n;
    }
  }

  // Pearson sample correlation coefficient
  double denom = (n*sxx - sx*sx) * (n*syy - sy*sy);
  r = denom == 0 ? 0 : (n*sxy - sx*sy) / sqrt(denom);

  // Estimate the signficance.
  if (n < 3) {
    *log_pv = 0;		// Can't estimate on fewer than 3 points.
  } else if (r < 1) {
    z = 0.5*sqrt(n-3)*log((1+r)/(1-r));
    *log_pv = log(0.5) + log_erfc(z/sqrt(2));	// Right tail of normal.
  } else {
    *log_pv = -BIG;		// p-value is 0
  }

  return r;
}

#ifdef PEARSON_MAIN
/************************************************************************/
/*
	Reads [<x> <y>]+ from named file or '-' (standard input)
	and prints the Pearson correlation coefficient and its p-value.
*/
/************************************************************************/
#define BUFSIZE 100
int main(int argc, char **argv) {
  FILE *fp = NULL;
  char *line = NULL;
  size_t len = 0;
  ssize_t read;
  int i = 0; 
  int n = 0;
  int asize = 0;
  double *x = NULL, *y = NULL;

  if (argc != 2) {
    fprintf(stderr, "Usage: pearson <filename>\n"
	"\tReads [<x> <y>]+ from named file or '-' (standard input)\n"
	"\tand prints the Pearson correlation coefficient and its p-value.\n"
    );
    return(1);
  }
  char *filename = argv[1];

  if (strcmp(filename, "-") == 0) {
    fp = stdin;
  } else {
    fp = fopen(filename, "r");
    if (fp == NULL) {
      fprintf(stderr, "Error: Cannot open file '%s' for reading.\n", filename);
      return(1);
    }
  }

  while ((read = getline(&line, &len, fp)) != -1) {
    i++;
    if (line[0] == '#') continue;               /* skip comments */
    if (n >= asize) {
      asize += BUFSIZE;
      Resize(x, asize, double);
      Resize(y, asize, double);
    }
    if (sscanf(line, "%lf %lf", &(x[n]), &(y[n])) != 2) {
      fprintf(stderr, "Error: line %d does not have the format '<x> <y>'\n%s\n", i, line);
      return(1);
    }
    if (0) fprintf(stderr, "n %d X %f Y %f asize %d\n", n, x[n], y[n], asize);
    n++;
  }
  if (line) free(line);
  double m, b, mse, log_pv;
  double cc = pearson_correlation(
    n,			/* number of points */
    x,			/* x values */
    y,			/* y values */
    &m,			/* slope */
    &b,			/* intercept */
    &mse,		/* mean-squared error */
    &log_pv 		/* log of Fisher transform p-value */
  );
  fprintf(stdout, "Pearson CC %.3f p-value ", cc);
  print_log_value(stdout, log_pv, 2);
  fprintf(stdout, "\n");

  return(0);
} // main
#endif
