#include "print.h"

void print_info( string name, string val ) {
  cout << " " << name << "  " << setfill(' ') << setw(28-name.length()) << "  " << setprecision(6) << val << setfill(' ') << endl; }

void print_info( string name, int val )	{
  cout << " " << name << "  " << setfill(' ') << setw(28-name.length()) << "  " << setprecision(6) << val << setfill(' ') << endl; }
void print_info( string name, int val, string add ) {
  cout << " " << name << "  " << setfill(' ') << setw(28-name.length()) << "  " << setprecision(6) << val << setfill(' ') << " " << add << endl; }

void print_info( string name, double val ) {
  cout << " " << name << "  " << setfill(' ') << setw(28-name.length()) << "  " << setprecision(6) << val << endl; }

void print_info( string name, double val, string add ) {
  cout << " " << name << "  " << setfill(' ') << setw(28-name.length()) << "  " << setprecision(6) << val << " " << add << endl; }

void print_info( string name ) {
  cout << " " << name << endl; }

void print_header( ) {
  cout << " " << setfill('-') << setw(30) << "  " << endl; }

void print_section( ) {
  cout << " " << setfill('#') << setw(30) << "  " << endl; }
