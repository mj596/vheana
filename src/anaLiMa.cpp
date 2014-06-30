#include "anaLiMa.h"

anaLiMa::anaLiMa( bool _verbose ) {
  verbose = _verbose;
  if( verbose ) { print_info("[LIMA] Creating LiMa analysis object"); }
  idata = new struct indata;
  odata = new struct outdata; }

anaLiMa::~anaLiMa() {
  if( verbose ) { print_info("[LIMA] Deleting LiMa analysis object"); }
  delete idata;
  delete odata;
  idata = NULL;
  odata = NULL; }

void anaLiMa::set( int _NOn, int _NOff, double _ExpOn, double _ExpOff ) {
  if( verbose ) { print_info("[LIMA] Setting Data"); }
  idata->NOn=_NOn;
  idata->NOff=_NOff;
  idata->ExpOn=_ExpOn;
  idata->ExpOff=_ExpOff; }

void anaLiMa::print() {
  print_info("[LIMA] Print data");
  print_info("[LIMA] NOn",idata->NOn);
  print_info("[LIMA] NOff",idata->NOff);
  print_info("[LIMA] ExpOn",idata->ExpOn);
  print_info("[LIMA] ExpOff",idata->ExpOff); }

double anaLiMa::analyse( ) {
  if( verbose ) { print_info("[LIMA] Start analysis"); }
  odata->alpha = idata->ExpOn/idata->ExpOff;
  odata->stddev = sqrt( idata->NOn+pow( odata->alpha, 2 )*idata->NOff );
  odata->excess = getExcess( );
  odata->significance = getSignificance( );
  if( verbose ) { 
    print_info("[LIMA] Analysis results");
    print_info("[LIMA] Alpha",odata->alpha);
    print_info("[LIMA] Excess",odata->excess);
    print_info("[LIMA] StdDev",odata->stddev);
    for( int i=1;i<=3;i++ ) {
      cout << " [LIMA] " << i << " sigma confidence interval (" << setprecision(3) << odata->excess-i*odata->stddev << ";" << odata->excess+i*odata->stddev << ")" << endl; }
    print_info("[LIMA] Significance",odata->significance);
    print_info("[LIMA] End analysis"); }

  return odata->excess;
}
  
double anaLiMa::getExcess( ) {
  return idata->NOn-(odata->alpha*idata->NOff); }

double anaLiMa::getSignificance( ) {
  return sqrt(2)*sqrt( idata->NOn*log(((1.0+odata->alpha)*idata->NOn)/(odata->alpha*(idata->NOn+idata->NOff)))+
		idata->NOff*log(((1.0+odata->alpha)*idata->NOff)/(idata->NOn+idata->NOff))
		       ); }

double* anaLiMa::getInterval( double sigma ) {
  static double interval[2];
  interval[0] = odata->excess-sigma*odata->stddev;
  interval[1] = odata->excess+sigma*odata->stddev;
  return interval; }

  
