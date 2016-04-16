#include "vector"
#include "body.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cstring>

Body::Body(double mass_, double a_, double e_, double M_, double E_, double * r_, double * v_, std::string name_)
{
    name_str = name_;//new char[32];

    // for (int i; i<32; i++) {
    //   name[i] = name_[i];
    // }
    //
    // strcpy(name_str, name_.c_str());

    mass = new double;
    a = new long double;
    e = new long double;
    ecosE = new double;
    esinE = new double;
    M = new double;
    M0 = new double;
    E = new double;
    E0 = new double;
    n = new double;

    r0 = new double[3];
    v0 = new double[3];

    for (int i=0; i<3; i++) {
      r0[i] = r_[i];
      v0[i] = v_[i];
    }

    *mass = mass_;
    *a = a_;
    *e = e_;
    *M0 = M_;
    *M = M_;
    *E0 = E_;
    *E = E_;
    r = r_;
    v = v_;



    // std::vector<double> r_add (3);
    // std::vector<double> v_add (3);
    //
    // for (int i=0; i<3; i++) {
    //   r_add[i] = r[i];
    //   v_add[i] = v[i];
    // }
    //
    // r_list.push_back(r_add);
    // v_list.push_back(v_add);
};


void Body::update_position(void) {
  // std::vector<double> r_add (3);

  for (int i=0; i<3; i++) {
    // r_add[i] = r[i];
    r0[i] = r[i];
  }

  // r_list.push_back(r_add);
}

void Body::update_velocity(void) {
  std::vector<double> v_add (3);

  // for (int i=0; i<3; i++) {
  //   v_add[i] = v[i];
  // }
  //
  // v_list.push_back(v_add);
}

void Body::output_position(void) {
  FILE * file_out;

  file_out = fopen(&("/home/hamish/Dropbox/LPL/classes/2016/dynamics/hw9/"+name_str)[0], "a+");
  fprintf(file_out, "%.30f \t %.30f \n", r[0],r[1]);
  fclose(file_out);

}
