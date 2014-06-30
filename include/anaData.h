#ifndef anaData_H
#define anaData_H

using namespace std;

class anaData
{
  bool verbose;

  struct indata {
    int NOn, NOff;
    double ExpOn, ExpOff; };
  indata *idata;

public:
  anaData( bool _verbose = false );
  ~anaData( );
  /* sets Non, Noff, ExpOn and ExpOff */
  int set( int _NOn, int _NOff, double _ExpOn, double _ExpOff );
  void print( );
  /* methods to get particular data */
  int get_NOn() { return idata->NOn; }
  int get_NOff() { return idata->NOff; }
  double get_ExpOn() { return idata->ExpOn; }
  double get_ExpOff() { return idata->ExpOff; }
};

#endif
