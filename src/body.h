#ifndef BODY_H
#define BODY_H

#include "vector"
#include <string>
// #include <cstring>

//Class to populate fields (i.e. velocity etc) at HALF grid spacing
class Body {
public:
  double * mass;
  double * a;
  double * e;
  double * ecosE;
  double * esinE;
  double * M0;
  double * M;
  double * E0;
  double * E;
  double * r;
  double * r0;
  double * v0;
  double * v;
  double * n;
  char * name;
  std::string name_str;

  // std::vector<std::vector<double> > r_list;
  // std::vector<std::vector<double> > v_list;

  Body(double, double, double, double, double, double*, double*, std::string);

  void update_position();

  void update_velocity();

  void output_position();


};

#endif
