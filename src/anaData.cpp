#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <cstdio>
#include "anaData.h"
#include "print.h"

using namespace std;

anaData::anaData( bool _verbose ) {
  verbose = _verbose;
  if( verbose ) { print_info("[DATA] Creating Data object"); }
  idata = new struct indata;
}

anaData::~anaData() {
  if( verbose ) { print_info("[DATA] Deleting Data object"); }
  delete idata;
  idata = NULL;
}

int anaData::set( int _NOn, int _NOff, double _ExpOn, double _ExpOff ) {
  if( verbose ){ print_info("[DATA] Setting Data"); }
  idata->NOn=_NOn;
  idata->NOff=_NOff;
  idata->ExpOn=_ExpOn;
  idata->ExpOff=_ExpOff;
  if( idata->ExpOn >= idata->ExpOff ) { return 1;}
  else { return 0; }
 }

void anaData::print() {
  print_info("[DATA] Print data");
  print_info("[DATA] NOn",idata->NOn);
  print_info("[DATA] NOff",idata->NOff);
  print_info("[DATA] ExpOn",idata->ExpOn);
  print_info("[DATA] ExpOff",idata->ExpOff); }
  
