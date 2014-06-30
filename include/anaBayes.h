#ifndef anaBayes_H
#define anaBayes_H

#include "print.h"
#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <cstdio>
#include <gmp.h>
#include <gmpxx.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include "alglibinternal.h"
#include "ap.h"
#include "specialfunctions.h"
#include <gsl/gsl_sf_erf.h>
#include <algorithm>
#include <numeric>
#include <vector>
#include <iterator>
#include <sys/stat.h>
#include <sys/types.h>
#include <sstream>

using namespace std;

class anaBayes
{
  struct indata {
    int NOn, NOff;
    double ExpOn, ExpOff; } *idata;
  struct outdata {
    double excess, alpha, significance; } *odata;
  
  gsl_matrix *vpi;
  int iSize, sSize;
  int posteriorAccuracy;
  gsl_matrix *integrated_posterior_prob;
  bool verbose;
  gsl_vector *vs, *vps;
  string saveFilename, saveDirname;
  string readFilename, readDirname;

  double Pi( int i, double sExpOn );  
  mpz_class Ci( int i );
  double Integrate( gsl_vector* vs, gsl_vector* vps, int iMin, int iMax );
  void saveVector( string filename, gsl_vector* x, gsl_vector* y, string dir );
  int readVector( );

 public:
  /* verbose sets information to be pronted on screen */
  anaBayes( bool _verbose );
  ~anaBayes();

  /* set data to be analysed */
  void set( int _NOn, int _NOff, double _ExpOn, double _ExpOff );

  /* print data */
  void print( );

  /* run analysis; only appropriate vectors are filled and posterior probability is calculated; savePosterior saves calculated data; if succeded returns 0, if not returns 1 */
  int analyse( bool _savePosterior=false );

  /* sets directories to read and save posterior distribution data; both are set to 'test' by default*/
  void setSaveDir( string _dir );
  void setReadDir( string _dir );

  /* reads previously calculated posterior probability; returns 1 if there is no file available and 0 if file has been found and read */
  int readAnalysis( );

  /* returns mode, median, mean and significance value */
  double getMode( );
  double getMedian( );
  double getMean( );
  double getSignificance( );

  /* methods to calculte credible bayesian intervals; before calculation calculator has to be initialized with first method */
  void initialize_credible_interval_calculator( );
  double* getInterval( double sigma );

  /* returns vs and vps vectors with posterior probability */
  gsl_vector* get_vs( );
  gsl_vector* get_vps( );
};

#endif
