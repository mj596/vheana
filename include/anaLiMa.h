#ifndef anaLiMa_H
#define anaLiMa_H

#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <cstdio>
#include "print.h"

using namespace std;

class anaLiMa
{
  struct indata {
    int NOn, NOff;
    double ExpOn, ExpOff; };
  struct outdata {
    double excess, alpha, significance, stddev; };
  indata *idata;
  outdata *odata;
  bool verbose;

 public:
  anaLiMa( bool _verbose );
  ~anaLiMa();
  void set( int _NOn, int _NOff, double _ExpOn, double _ExpOff );
  void print( );
  
  /* analyse does all the simple analysis including writting results on screen if verbose option is true */
  double analyse( );
  double getExcess( );
  double getSignificance( );
  double* getInterval( double sigma );
};

#endif
