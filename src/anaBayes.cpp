#include "anaBayes.h"

anaBayes::anaBayes( bool _verbose ) {
  verbose = _verbose;
  if( verbose ) { print_info("[BAYES] Creating Bayesian analysis object"); }
  idata = new struct indata;
  odata = new struct outdata;
  posteriorAccuracy = 10;
  saveDirname = "test";
  saveFilename = "test.dat";
  readDirname = "test";
  readFilename = "test.dat"; }

anaBayes::~anaBayes() {
  if( verbose ) { print_info("[BAYES] Deleting Bayesian analysis object"); }
  delete idata;
  idata = NULL;
  delete odata;
  odata = NULL;
  gsl_vector_free( vs );
  gsl_vector_free( vps );
  gsl_matrix_free( vpi );
}

void anaBayes::setSaveDir( string _dir ){
  saveDirname = _dir; }

void anaBayes::setReadDir( string _dir ){
  readDirname = _dir; }

void anaBayes::set( int _NOn, int _NOff, double _ExpOn, double _ExpOff ) {
  if( verbose ) { print_info("[BAYES] Setting Data"); }
  idata->NOn=_NOn;
  idata->NOff=_NOff;
  idata->ExpOn=_ExpOn;
  idata->ExpOff=_ExpOff;
  
  iSize = idata->NOn+1; /* size of i vectors from 0,1,..., NOn */
  sSize = (idata->NOn-((idata->ExpOn/idata->ExpOff)*idata->NOff));
  if( sSize>10 ) { sSize = 10; }
  sSize *= 20*posteriorAccuracy;
  if( sSize <= 0 ) { sSize = 20*posteriorAccuracy; }
  if( verbose ) { 
    print_info("[BAYES] iSize",iSize);
    print_info("[BAYES] sSize",sSize);
  }
  
  stringstream s;
  s << "anaBayes" << "_ON" << _NOn <<"_OFF" << _NOff << "_ALPHA" << std::fixed << setprecision(3) << _ExpOn/_ExpOff << ".dat";
  saveFilename = s.str();

  /* initialize vectors */
  vpi = gsl_matrix_alloc(sSize, iSize);
  vs = gsl_vector_alloc(sSize);
  vps = gsl_vector_alloc(sSize);
  gsl_vector_set_zero( vs );
  gsl_vector_set_zero( vps );
  gsl_matrix_set_zero( vpi );
}

void anaBayes::print( ) {
  if( verbose ) { 
    print_info("[BAYES] Print data");
    print_info("[BAYES] NOn",idata->NOn);
    print_info("[BAYES] NOff",idata->NOff);
    print_info("[BAYES] ExpOn",idata->ExpOn);
    print_info("[BAYES] ExpOff",idata->ExpOff); }
}

void anaBayes::saveVector( string filename, gsl_vector* x, gsl_vector* y, string dir ) {
  mkdir( dir.c_str(), 0777 );
  FILE *output;
  stringstream s;
  s << dir << "/" << filename;
  output = fopen(s.str().c_str(),"w");
  
  for( int i=0;i<y->size;i++ ) { fprintf( output,"%le %le\n", gsl_vector_get(x,i), gsl_vector_get(y,i) ); }
  
  fclose(output); }

int anaBayes::readVector( ) {
  FILE *input;
  input = fopen( readFilename.c_str(), "r" );
  if ( input == NULL ) { return 1; }
  else
    {
      double x, y;
      for( int i=0;i<sSize;i++ )
	{
	  if( !fscanf( input,"%le %le\n",&x,&y ) ) 
	    if( verbose ) { print_info("[BAYES] Error reading file."); }
	  gsl_vector_set( vs, i, x );
	  gsl_vector_set( vps, i, y );
	}
      fclose(input);
      return 0;
    }
}

int anaBayes::readAnalysis( ) {
  stringstream s;
  s << readDirname << "/" << "anaBayes" << "_ON" << idata->NOn <<"_OFF" << idata->NOff << "_ALPHA" << std::fixed << setprecision(3) << idata->ExpOn/idata->ExpOff << ".dat";
  readFilename = s.str();
  if( verbose ) { print_info("[BAYES] Attempting to find file",readFilename); }
  if( readVector( ) ) { 
    if( verbose ) { print_info("[BAYES] File does not exist."); }
    return 1; }
  else { 
    if( verbose ) { print_info("[BAYES] File found."); }
    return 0; }
}

int anaBayes::analyse( bool _savePosterior ) {
  if( verbose ) { print_info("[BAYES] Start analysis"); }
  if( idata->NOn == 0 )
    { 
      if( verbose ) { print_info("[BAYES] NOn=0. Quit."); }
      return 1;
    }
  else {
    for( int i=0;i<sSize;i++ ) { gsl_vector_set( vs, i, (1.0/posteriorAccuracy)*(i+1)); }
    for( int i=0;i<sSize;i++ ) {
      gsl_matrix_set( vpi, i, 0, gsl_expm1(-gsl_vector_get(vs, i))+1.0 );
      for( int j=1;j<iSize;j++ ) { gsl_matrix_set( vpi, i, j, Pi( j, gsl_vector_get(vs, i) ) ); } }
    
    mpf_class CiSum = 0.0; /* CiSum  should be = 1.0 */
    /* alocate vectors for Ci */
    mpz_class vC[idata->NOn+1];
    mpf_class vCi[idata->NOn+1];
    mpz_class Csum = 0;
    for( int i=0; i<iSize; i++ ) {
      vC[i] = Ci(i);
      Csum += vC[i]; }
    for( int k=0; k<iSize; k++ ) {
      vCi[k] = vC[k]/(mpf_class)Csum;
      CiSum += vCi[k]; }
    if( verbose ) { cout << " [BAYES] CiSum " << CiSum << endl; }
    
    mpf_class vpsSum;
    for( int i=0;i<sSize;i++ ) {
      vpsSum = 0.0;
      for( int j=0;j<iSize;j++ )	{ vpsSum += idata->ExpOn*vCi[j]*gsl_matrix_get( vpi, i, j ); }
      gsl_vector_set( vps, i, mpf_get_d(vpsSum.get_mpf_t()) ); }
    
    /* posterior Norm check and normalize */
    double posteriorNorm = Integrate(vs,vps,0,vs->size);
    if( verbose ) { print_info("[BAYES] Posterior norm check",posteriorNorm); }
    if( verbose ) { print_info("[BAYES] Normalize posterior"); }
    gsl_vector_scale (vps, 1.0/posteriorNorm );
    if( verbose ) { print_info("[BAYES] Posterior norm check",Integrate(vs,vps,0,vs->size)); }
    
    /* save posterior distribution */
    if( _savePosterior ) {
      if( verbose ) { print_info("[BAYES] Save posterior distribution",saveDirname+"/"+saveFilename); }
    saveVector( saveFilename, vs, vps, saveDirname );
    }
    return 0;
  }     
}
//
//   /* Show results */
//   if( verbose ) { 
//     print_info("[BAYES] Analysis results");   
//     print_info("[BAYES] Mean",MeanValue(vs,vps,sSize));
//     print_info("[BAYES] Mode",ModeValue(vs,vps,sSize));
//   }
//
//   /* Calculate credible intervals */
//   if( verbose )
//     {
//       initialize_credible_interval_calculator();
//       calculate_credible_interval( 1.0 );
//       calculate_credible_interval( 2.0 );
//       calculate_credible_interval( 3.0 );
//     }
//   
//   if( verbose ) { print_info("[BAYES] Median",MedianValue(vs,vps,sSize)); }
//
//   /* Calculate significance */
//   double bayesianSignificance;
//   //bayesianSignificance = calculate_significance();
//
//   if( verbose ) { 
//     print_info("[BAYES] Significance",bayesianSignificance);
//     print_info("[BAYES] End analysis"); }
//
//   return ModeValue(vs,vps,sSize);
//}

void anaBayes::initialize_credible_interval_calculator( ) {
  integrated_posterior_prob = gsl_matrix_alloc(vs->size,vs->size);
  gsl_matrix_set_zero( integrated_posterior_prob );
  
  for( int i=0;i<vs->size;i++ )
    for( int j=1;j<vs->size;j++ )
      if( i < j ) {
	gsl_matrix_set( integrated_posterior_prob, i, j, gsl_matrix_get( integrated_posterior_prob, i, j-1)+gsl_vector_get( vps, j-1 ) ); }
  gsl_matrix_scale( integrated_posterior_prob, gsl_vector_get( vs, 1 )-gsl_vector_get( vs, 0 ) );
}

double anaBayes::Pi( int i, double sExpOn ) { return alglib::poissondistribution(i, sExpOn)-alglib::poissondistribution(i-1, sExpOn); }

mpz_class anaBayes::Ci( int i ) {
  mpz_class _temp = 1;
  mpz_class _i = (unsigned int)i;
  mpz_class _NOn = (unsigned int)idata->NOn;
  mpz_class _NOff = (unsigned int)idata->NOff;
  
  for ( mpz_class k = 1;k <=_NOff;k++ ) { _temp *= _NOn-_i+k; }

  mpz_class __temp;
  mpz_class __1ExpOff = 1 + idata->ExpOff;
  mpz_pow_ui(__temp.get_mpz_t(), __1ExpOff.get_mpz_t(), (unsigned long int)i);
  
  double __tempd;
  double Td;
  Td = 1.0+(idata->ExpOff/idata->ExpOn);
  __tempd = pow(Td,i);
  return __tempd*_temp; }

double anaBayes::getMode( ) {
  double temp = 0.0;
  int iMax;
  double mode;

  for( int i=0; i<sSize; i++ ) {
    if( gsl_vector_get( vps, i ) >= temp ) {
      temp = gsl_vector_get( vps, i );
      iMax = i; } }

  if( iMax == 0 ) { mode = 0.0; }
  else { mode = gsl_vector_get( vs, iMax ); }

  if( verbose ) { print_info("[BAYES] Mode",mode); }
  return mode; }

double anaBayes::getMean( ) {
  double temp = 0.0;
  for( int i=0; i<sSize; i++ ) { temp += gsl_vector_get( vs, i )*gsl_vector_get( vps, i ); }
  if( verbose ) { print_info("[BAYES] Mean",temp*(gsl_vector_get( vs, 1 )-gsl_vector_get( vs, 0 ))); }
  return temp*(gsl_vector_get( vs, 1 )-gsl_vector_get( vs, 0 )); }

double anaBayes::getMedian( ) {
  double Intmin = 0.0, Intmax = 0.0;
  int med_min = 1;
  int med_max = sSize-1; 
  double median;

  do {
    Intmin = Integrate(vs,vps,0,med_min);
    if( Intmin < 0.5 ) { med_min += 1; }
    else { median = gsl_vector_get( vs, med_min ); }
    Intmax = Integrate(vs,vps,med_max,sSize-1);
    if( Intmax < 0.5 ) { med_max -= 1; }
    else { median = gsl_vector_get( vs, med_min ); } }
  while( Intmin < 0.5 && Intmax < 0.5 );
  
  if( verbose ) { print_info("[BAYES] Median",median); }
  return median; }

double anaBayes::Integrate( gsl_vector* vs, gsl_vector* vps, int iMin, int iMax ) {
  double temp = 0.0;
  for( int i=iMin; i<iMax; i++ ) { temp += gsl_vector_get( vps, i ); }
  return (gsl_vector_get( vs, 1 )-gsl_vector_get( vs, 0 ))*temp; }

double* anaBayes::getInterval( double sigma ) {
  double K = gsl_sf_erf(sigma/sqrt(2));
  gsl_vector *posterior_path_prob = gsl_vector_alloc( vs->size );
  gsl_vector *posterior_path_length = gsl_vector_alloc( vs->size );
  gsl_vector *posterior_path_j = gsl_vector_alloc( vs->size );
  gsl_vector_set_zero( posterior_path_prob );
  gsl_vector_set_all( posterior_path_length, 9999 );
  gsl_vector_set_zero( posterior_path_j );

  for( int i=0;i<vs->size;i++ ) {
    for( int j=0;j<vs->size-1;j++ ) {
      if( gsl_matrix_get(integrated_posterior_prob, i, j) <= K && gsl_matrix_get(integrated_posterior_prob, i, j+1) >= K && i<j ) {
	gsl_vector_set( posterior_path_prob, i, 0.5*(gsl_matrix_get(integrated_posterior_prob, i, j)+gsl_matrix_get(integrated_posterior_prob, i, j+1)));
	gsl_vector_set( posterior_path_length, i, gsl_vector_get(vs,j)-gsl_vector_get(vs,i));
	gsl_vector_set( posterior_path_j, i, gsl_vector_get(vs,j));
	break; }
    } }

  double xmin;
  if( gsl_vector_get( vs, gsl_vector_min_index( posterior_path_length )) == gsl_vector_get( vs, 0 ) ) { xmin = 0.0; }
  else { xmin = gsl_vector_get( vs, gsl_vector_min_index( posterior_path_length )); }
  if( verbose ) { 
    cout << " [BAYES] " << sigma << " sigma credible interval (" << xmin << ";" << gsl_vector_get( posterior_path_j, gsl_vector_min_index( posterior_path_length )) << ")" << endl; }
  
  static double returnArray[2];
  returnArray[0] = xmin;
  returnArray[1] = gsl_vector_get( posterior_path_j, gsl_vector_min_index( posterior_path_length ));

  gsl_vector_free( posterior_path_prob );
  gsl_vector_free( posterior_path_length );
  gsl_vector_free( posterior_path_j );

  return returnArray;
 }
  
double anaBayes::getSignificance( ) {
//  std::cout << gsl_vector_get(vps,0) << std::endl;
//  std::cout << 1.0-gsl_vector_get(vps,0) << std::endl;
  
  if( gsl_vector_get(vps,sSize-1)>=gsl_vector_get(vps,0) )
    { print_info("[BAYES] Error: too low accuracy to compute significance.");
      return -1.0; }
  else {
    int imaxSig;
    for( int i = gsl_vector_max_index( vps );i<sSize-1;i++ ) {
      if( (gsl_vector_get(vps,0) <= gsl_vector_get(vps,i)) && (gsl_vector_get(vps,0) >= gsl_vector_get(vps,i+1)) ) { imaxSig = i;
	break; } }

    double maxSingificance = 5.0;
    long sigPrecision = 1.0e2;
    double bayesIntervalProbability;
    /* the case when posterior has a maximum in 0 */
    if( imaxSig == 0 )
      { bayesIntervalProbability = 0.0; }
    else
      {
	bayesIntervalProbability = 0.5*(Integrate(vs,vps,0,imaxSig)+Integrate(vs,vps,0,imaxSig+1));
      }
    
    if( bayesIntervalProbability > 1.0 ) bayesIntervalProbability = 1.0;
    if( verbose ) { print_info("[BAYES] Sig Credible interval",bayesIntervalProbability); }
    
    if( bayesIntervalProbability <= 1.0e-3 ) { return 0; }
    else {
 
  gsl_vector *significance_vector = gsl_vector_alloc(maxSingificance*sigPrecision);
  gsl_vector *erf_vector = gsl_vector_alloc(maxSingificance*sigPrecision);
      
  double bayesianSignificance = 0.0;
  for( int i=0;i<maxSingificance*sigPrecision;i++ ) { 
    gsl_vector_set( significance_vector, i, (1.0/((double)sigPrecision))*i );
    gsl_vector_set( erf_vector, i, gsl_sf_erf(gsl_vector_get( significance_vector, i )/sqrt(2)) ); } 
  
  for( int i=0;i<maxSingificance*sigPrecision-1;i++ ) { 
    if( gsl_vector_get( erf_vector, i ) <= bayesIntervalProbability && gsl_vector_get( erf_vector, i+1 ) >= bayesIntervalProbability ) {
      bayesianSignificance = 0.5*(gsl_vector_get( significance_vector, i )+gsl_vector_get( significance_vector, i+1 ));
	break; } } 
    
  if( bayesianSignificance >= 5.0 || bayesianSignificance == 0.0 ) { bayesianSignificance = -1.0; }
  
  gsl_vector_free( significance_vector );
  gsl_vector_free( erf_vector );
    
  if( verbose ) { print_info("[BAYES] Significance",bayesianSignificance); }
  
  return bayesianSignificance;
    }
  }
}
  
gsl_vector* anaBayes::get_vs( ) { return vs; }
gsl_vector* anaBayes::get_vps( ) { return vps; }
