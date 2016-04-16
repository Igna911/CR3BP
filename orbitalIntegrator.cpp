#ifndef ORBITALINTEGRATOR_H
#define ORBITALINTEGRATOR_H

#include "body.h"
#include "orbitalIntegrator.h"
#include <math.h>
#include <quadmath.h>
#include <iostream>
#include <iomanip>

void cross(double * axb, double * a, double * b) {
  axb[0] = a[1]*b[2] - a[2]*b[1];
  axb[1] = -(a[0]*b[2] - a[2]*b[0]);
  axb[2] = a[0]*b[1] - a[1]*b[0];
};

void dot(double * adotb, double * a, double * b) {
  *adotb = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
  // axb[1] = -(a[0]*b[2] - a[2]*b[0]);
  // axb[2] = a[0]*b[1] - a[1]*b[0];
};

double mag(double * array) {
  return sqrt(pow(array[0],2.0) + pow(array[1],2.0) + pow(array[2],2.0));
};

void step_M(Body * body, double dt) {
  double pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303819644288109756659334461284756;
  *body->M0 = *body->M;
  *body->M = *body->M + *body->n * dt;
  //
  if (*body->M0 >= 2.0*pi) {
    *body->M0 -= 2.0*pi;
    *body->M -= 2.0*pi;
  }

}

void calc_orb_elements(Body * body, double mu) {
  double r0 = mag(body->r);
  double v0 = mag(body->v);
  double * h = new double[3];
  double * vxh = new double[3];
  double * rdotv = new double[3];
  double * e_vec = new double[3];
  double p = 0.0;

  *body->a = 1.0/(2.0/r0 - pow(v0, 2.0)/mu);


  *body->n = sqrt(mu/pow(*body->a,3.0));

  dot(rdotv, body->r, body->v);

  *body->ecosE = 1.0 - r0/(*body->a);
  *body->esinE = *rdotv/(*body->n * pow(*body->a,2.0));

  // std::cout<<*body->a<<"\t"<<*body->n<<"\t"<<*body->esinE<<"\t"<<*body->ecosE<<std::endl;
  //
  cross(h,body->r,body->v);
  cross(vxh,body->v,h);

  for (int i=0; i<3; i++) {
    e_vec[i] = vxh[i] - body->r[i]/r0;
  }
  //
  *body->e = mag(e_vec);

  // std::cout<<*body->e<<std::endl;




  // if (*body->e < 1e-5) *body->e = 0.0;



  //
  // p = pow(mag(h),2.0);
  //
  // // std::cout<<p/(*body->a)<<"\t"<<sqrt(1.0 - p/(*body->a))<<std::endl;
  //
  // if (p/(*body->a) > 1.0) *body->e;
  // else *body->e = sqrt(1.0 - p/(*body->a));

  // *body->e = sqrt(1.0 - pow(mag(h),2.0)/(*body->a));
  // std::cout<<*body->e<<std::endl;
  // std::cout<<fabs(*body->e - sqrt(1.0 - p/(*body->a)))<<std::endl;

  delete[] rdotv;
  delete[] h;
  delete[] vxh;
  delete[] e_vec;
}

double kepler(double dE, double ecosE, double esinE) {
  // return E - e*sin(E);
  return dE - ecosE*sin(dE) + esinE*(1.0 - cos(dE));
}


double kepler_deriv(double dE, double ecosE, double esinE) {
  // return 1.0 - e*cos(E);
  return 1.0 - ecosE*cos(dE) + esinE * sin(dE);
}

void newton_raphson(double  * dE, double * ecosE, double * esinE, double dM) {
  // return E - (kepler(M,e,E)-M)/kepler_deriv(e,E);
  *dE = *dE - (kepler(*dE, *ecosE, *esinE) - dM)/kepler_deriv(*dE, *ecosE, *esinE);
}

void calc_E(Body * body, double dt) {
  double pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303819644288109756659334461284756;
  double e = 0.0;
  double E0 = *body->E;//*body->M + sin(*body->M)*0.85*(*body->e);
  double E = 0.0;
  double dM = *body->n * dt;//*body->M - *body->M0;
  // double dE = *body->E - *body->E0;//(*body->M + sin(*body->M)*0.85*(*body->e)) - E0;
  // double dE0 = 0.0;
  double * dE = new double;
  double * esinE = new double;
  double * ecosE = new double;
  int count = 0;

  *esinE = *body->esinE;
  *ecosE = *body->ecosE;

  *dE = (dM+ sin(dM)*0.85*(*body->e));// - E0;

  double test = fabs(kepler(*dE, *ecosE, *esinE) - dM);
  // std::cout<<std::setprecision(25)<<"Entering: "<<std::endl;
  while (test > 1e-18 && !isnan(test) && count < 10) {


    newton_raphson(dE, body->ecosE, body->esinE, dM);
    // std::cout<<std::setprecision(10)<<"HELP"<<std::endl;
    // if (fabs(kepler(*dE, *ecosE, *esinE) - dM) > 1.0) {
    //   *dE = (dM+ sin(dM)*0.851*(*body->e));
    //   // std::cout<<"HALP";
    //   newton_raphson(dE, body->ecosE, body->esinE, dM);
    test = fabs(kepler(*dE, *ecosE, *esinE) - dM);
  // }
    // std::cout<<"HALP"<<std::endl;
    // dE0 = dE;
    count ++;
  }
  // std::cout<<std::setprecision(25)<<count<<std::endl;


  *body->E0 = *body->E;
  *body->E += *dE;

  *body->M0 = *body->M;
  *body->M += dM;

  // std::cout<<std::setprecision(10)<<*body->M<<"\t"<<*body->E<<std::endl;

  if (*body->E0 >= 2*pi) {
    *body->E -= 2*pi;
    *body->E0 -= 2*pi;
    *body->M -= 2*pi;
    *body->M0 -= 2*pi;
  }
  // std::cout<<*body->E<<std::endl;
}


// double kepler(double E, double e) {
//   return E - e*sin(E);
//   // return dE - e*cos(E0)*sin(dE) + e*sin(E0)*(1.0 - cos(dE));
// }
//
//
// double kepler_deriv(double E, double e) {
//   return 1.0 - e*cos(E);
//   // return 1.0 - e*cos(E0)*cos(dE) + e*sin(E0) * sin(dE);
// }
//
// double newton_raphson(double M, double e, double E) {
//   return E - (kepler(E,e)-M)/kepler_deriv(E,e);
//   // return dE - (kepler(dE, e, E0) - dM)/kepler_deriv(dE, e, E0);
// }
//
// void calc_E(Body * body) {
//   double E = *body->M + sin(*body->M)*0.85*(*body->e);
//   // double E = 0.0;
//   double M = *body->M;// - *body->M0;
//   // double dE = *body->E - *body->E0;//(*body->M + sin(*body->M)*0.85*(*body->e)) - E0;
//   // double dE0 = 0.0;
//   // double dE = (dM+ cos(M)*0.85*(*body->e));// - E0;
//
//   // E = *body->E0;
//   // std::cout<<std::setprecision(25)<<*body->e<<std::endl;
//   // std::cout<<std::setprecision(25)<<E<<"\t"<<M<<std::endl;
//   // std::cout<<std::setprecision(25)<<kepler(E,*body->e)<<std::endl;
//   while (fabs(M - kepler(E,*body->e)) > 1e-12) {
//     // std::cout<<std::setprecision(10)<<"HELP"<<std::endl;
//     E = newton_raphson(M, *body->e, E);
//     // dE0 = dE;
//   }
//   *body->E0 = *body->E;
//   *body->E = E;
// }








void calc_fg(Body * body, double *f, double * g, double dt, double dE) {
  *f = 1.0 + *body->a/mag(body->r) * (cos(dE) - 1.0);
  *g = dt + (sin(dE) - (dE))/(*body->n);
}

void calc_fdot_gdot(Body * body, double *fdot, double * gdot, double dE, double fp) {
  double * r0 = new double[3];

  for (int i=0; i<3; i++) {
    // r0[i] = body->r_list[body->r_list.size()-2][i];
    r0[i] = body->r0[i];
  }

  // *fdot = -(pow(*body->a,2.0) / (mag(body->r) * mag(r0))) * (*body->n)*sin(*body->E - *body->E0);

  *fdot = -(*body->a / (fp * mag(r0))) * (*body->n)*sin(dE);


  // *gdot = 1.0 + (*body->a)/mag(body->r) * (cos(*body->E - *body->E0) - 1.0);

  *gdot = 1.0 + (cos(dE) - 1.0)/fp;


  delete[] r0;
}

void update_position(Body * body, double dt) {
  double * f = new double;
  double * g = new double;
  double dE = (*body->E - *body->E0);

  // std::cout<<"Updating position: "<<std::endl<<std::endl;

  calc_fg(body, f, g, dt, dE);

  for (int i=0; i<3; i++) {
    body->r0[i] = body->r[i];
    body->r[i] = (*f) * body->r[i] + (*g) * body->v[i];
    // std::cout<<std::setprecision(10)<<"New position: "<<body->r[i]<<std::endl;
  }

  // body->update_position();

  delete[] f;
  delete[] g;

}

void update_velocity(Body * body, double force, double dt) {
  double * fdot = new double;
  double * gdot = new double;

  double dE = (*body->E - *body->E0);
  double fp = 1.0 - *body->ecosE *cos(dE) + *body->esinE * sin(dE);

  // std::cout<<"Updating velocity: "<<std::endl<<std::endl;

  calc_fdot_gdot(body,fdot,gdot, dE, fp);

  for (int i=0; i<3; i++) {
    body->v0[i] = body->v[i];
    body->v[i] = (*fdot) * body->r0[i] + (*gdot) * body->v[i];
    // std::cout<<std::setprecision(10)<<"New velocity: "<<body->v[i]<<std::endl;
  }

  // body->update_velocity();

  delete[] fdot;
  delete[] gdot;

  // -(pow(*body->a,2.0) / ()
}

void apply_perturbation(Body * b1, Body * b2, Body * b3, double dt) {
  double * acceleration = new double[3];
  double direct = 0.0;
  double indirect = 0.0;
  double * bottom = new double[3];
  double sun = 0.0;

  for (int i=0; i<3; i++) {
    bottom[i] = b1->r[i] - b2->r[i];
  }

  for (int i=0; i<3; i++) {
    direct = (b1->r[i] - b2->r[i])/pow(mag(bottom),3.0);
    indirect = (b2->r[i])/pow(mag(b2->r),3.0);

    // std::cout<<direct;
    // sun = -(*b3->mass)*b1->r[i] / pow(mag(b1->r),3.0);
    acceleration[i] = -(*b2->mass) * (direct - indirect);

    // b1->v[i] = b1->v_list[b1->v_list.size()-1][i] + dt*acceleration[i];
    b1->v[i] += dt*acceleration[i];
  }

  // b1->update_velocity();

  delete[] acceleration;
  delete[] bottom;

}

void update_planet_position(Body * body) {
  double x = *body->a * cos(*body->M);
  double y = *body->a * sin(*body->M);

  body->r[0] = x;
  body->r[1] = y;
  body->r[2] = 0.0;

  body->update_position();


}

#endif
