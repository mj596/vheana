#ifndef _PRINT_H
#define _PRINT_H 1

#include "math.h"
#include <stdio.h>
#include <iostream>
#include <iomanip>

using namespace std;

void print_info( string name, double val, string add );
void print_info( string name, int val, string add );
void print_info( string name, double val );
void print_info( string name, int val );
void print_info( string name, string val );
void print_info( string name );
void print_header( );
void print_section( );

#endif /* _PRINT_H */