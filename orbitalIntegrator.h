#ifndef ORBITALINTEGRATOR_H
#define ORBITALINTEGRATOR_H

#include "body.h"

double pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303819644288109756659334461284756;

void cross(double * axb, double * a, double * b);

double mag(double * array);

void calc_orb_elements(Body * body, double mu);

void step_M(Body *, double dt);

double newton_raphson(double M, double e);

double kepler(double *, double *, double *);

double kepler_deriv(double, double, double e);

void calc_E(Body *, double);

// TO DO

void calc_fg(Body *, double *f, double *g, double t1, double t2);

void calc_fdot_gdot(Body * body, double *fdot, double * gdot, double, double);

void update_position(Body *, double dt);

void update_velocity(Body *, double force, double dt);

void apply_perturbation(Body *, Body *, Body *, double);

void update_planet_position(Body *);

#endif
