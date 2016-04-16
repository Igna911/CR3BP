#include "body.h"
#include "orbitalIntegrator.h"
#include "vector"
#include <math.h>
#include <iostream>
#include <iomanip>
#include <string>



int main(void)
{
	double pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303819644288109756659334461284756;
	double m2 = 1e-3;
	double m1 = 1.0 - m2;
	double t1;
	double t2;
	double planet_a = 1.0;
	int iters = 1;
	// double particle_a[] = {0.98,0.96,0.94,0.92,0.90,0.88,0.86,0.84,0.82,0.80};
	double B = 0.99999;
	double A = 0.8;
	// double particle_a_new = pow(1.0/pow(4.0/3.0,2.0),1./3.0);
//	double particle_a_new = 0.9;//0.7631428283688879;////0.8254818122236567;// -- 3/4 resonance!
	// 0.7631428283688879 -- 3/2 resonance;
	double particle_a_new = (A+B)/2.0;
	double particle_a_old = 0.0;
	double P = 2*pi;
	double dt = P/8000.0;
	double integrate_time = 1e4;
	bool stable = true;
	bool bisect = true;
	int count = 0;
	int count2 = 0;
	int count_old = 0;

	int A_count = 563754;
	int B_count = 866117;

	double inc = 0.0000001;

	std::string planet_name = "planet_pos.txt";
	std::string particle_name = "particle_pos.txt";
	std::string star_name = "star_pos.txt";

	double * r_planet;
	double * r_particle;

	double * v_planet;
	double * v_particle;


	// if (!bisect) particle_a_new = B;
	// else particle_a_new = 0.5*(A+B);

	while (bisect) {//fabs(particle_a_new - particle_a_old) > 0.05) {
		stable = true;

		// if (!bisect) particle_a_new -= inc;

		r_planet = new double[3];  // Allocate n ints and save ptr in a.
		r_planet[0] = planet_a;
		r_planet[1] = 0.0;
		r_planet[2] = 0.0;

		r_particle = new double[3];  // Allocate n ints and save ptr in a.
		r_particle[0] = particle_a_new;
		r_particle[1] = 0.0;
		r_particle[2] = 0.0;

		v_planet = new double[3];  // Allocate n ints and save ptr in a.
		v_planet[0] = 0.0;
		v_planet[1] = 0.0;
		v_planet[2] = 0.0;

		v_particle = new double[3];  // Allocate n ints and save ptr in a.
		v_particle[0] = 0.0;
		v_particle[1] = sqrt(1.0*m1/(mag(r_particle)));//sqrt(1.0*m2/(mag(r_planet) - mag(r_particle)));
		// v_particle[1] = sqrt(1.0*m2/(mag(r_planet) - mag(r_particle)));
		v_particle[2] = 0.0;

		Body * particle = new Body(0, particle_a_new, 0.0, 0, 0, r_particle, v_particle, particle_name);
		Body * planet = new Body(m2, planet_a, 0, 0, 0, r_planet, v_planet, planet_name);
		Body * star = new Body(m1, 0.0, 0, 0, 0, r_planet, v_planet, star_name);

		*planet->n = 1.0;
		// *particle->n =1.0;

		t1 = 0.0;
		t2 = 0.0;
		count = 0;

		// calc_orb_elements(particle, 1.0-m2);
		while (t2 < integrate_time*P) {
		// while (count < 10) {

			// calculate orbital elements of particle
			calc_orb_elements(particle, 1.0-m2);

			// Increment time
			t1 = t2;
			t2 += 0.5*dt;

			// Step planet forward half a time step

			step_M(planet,0.5*dt);

			update_planet_position(planet);

			// Step particle forward half a time step

			// step_M(particle,0.5*dt);
			// *particle->E = *particle->M;
			// *particle->E0 = *particle->M0;

			calc_E(particle, 0.5*dt);

			// calc_orb_elements(particle, 1.0);
			// std::cout<<*particle->E<<"\t"<<*particle->M<<std::endl;

			// Find new particle position and velocity

			update_position(particle,0.5*dt);
			update_velocity(particle, 0.0, 0.0);

			// for (int i = 0; i<3; i++) std::cout<<particle->v[i]<<"\t";
			// std::cout<<std::endl;

			// // std::cout<<std::setprecision(20)<<*particle->n<<std::endl;
			//
			// // return 0;
			// // Apply perturbation to particle
			//
			// // std::cout<<*particle->a;
			//
			apply_perturbation(particle,planet,star,dt);
			// update_velocity(particle, 0.0, 0.0);

			// Update orbital elements with new perturbed velocity

			calc_orb_elements(particle, 1.0 - m2);

			// Increment time
			t1 = t2;
			t2 += 0.5*dt;

			// Step planet forward half a time step
			step_M(planet,0.5*dt);
			update_planet_position(planet);


			// Step particle forward half a time step
			// step_M(particle,0.5*dt);
			// *particle->E = *particle->M;
			calc_E(particle, 0.5*dt);

			// *particle->E = *particle->M;
			// *particle->E0 = *particle->M0;

			// Update particle position and velocity with new orbital elements
			update_position(particle,0.5*dt);
			update_velocity(particle, 0.0, 0.0);



			count ++;

			//mag(particle->r) > mag(planet->r) ||
			if (isnan(particle->r[0]) || mag(particle->r) > mag(planet->r) || *particle->e >= 1.0) {
				std::cout<<"FAILED AT "<<count<<std::endl;
				// std::cout<<"e: "<<*particle->e<<std::endl;
				// std::cout<<"a: "<<*particle->a<<std::endl;
				// std::cout<<"M: "<<*particle->M<<std::endl;
				// std::cout<<"E: "<<*particle->E<<std::endl;
				// std::cout<<"x: "<<particle->r[0]<<std::endl;
				// std::cout<<"y: "<<particle->r[1]<<std::endl;
				stable = false;
				// return 0;
				break;
			}

//			if (count%40 == 0) {
//				planet->output_position();
//				particle->output_position();
//
//
//				// std::cout<<std::setprecision(20)<<fabs(*particle->M - *particle->M0)<<std::endl;
//				// std::cout<<t2/P<<std::endl;
//			}

			if (count%100000 == 0) {
				// std::cout<<std::setprecision(10)<<(double)*particle->a<<"\t"<<std::setprecision(10)<<*particle->esinE<<std::endl;
				std::cout<<std::setprecision(5)<<"Complete: "<<t2/(integrate_time*P) * 100<<"% \t"<<t2/P<<" years"<<std::endl;
			}

		}

		std::cout<<"e: "<<*particle->e<<std::endl;
		std::cout<<"a: "<<*particle->a<<std::endl;
		std::cout<<"M: "<<*particle->M<<std::endl;
		std::cout<<"E: "<<*particle->E<<std::endl;
		std::cout<<"x: "<<particle->r[0]<<std::endl;
		std::cout<<"y: "<<particle->r[1]<<std::endl;

		// return 0;

		particle_a_old = particle_a_new;

		 if (stable) {
//			std::cout<<"Stable: "<<particle_a_new<<" at count: "<<count<<std::endl;
//			if (bisect) {
			A = particle_a_new;
			std::cout<<"Orbit a: "<<particle_a_new<<" is stable."<<std::endl;

			particle_a_new = 0.5*(A+B);
			std::cout<<"Moving to orbit a="<<particle_a_new<<std::endl<<std::endl;
//			}

			//
		 }
		 else {
			std::cout<<"Unstable: "<<particle_a_new<<" at count: "<<count<<std::endl;
//			if (bisect) {
//				if (count > B_count) {
			B = particle_a_new;
			std::cout<<"Orbit a: "<<particle_a_new<<" is unstable."<<std::endl;
			particle_a_new = 0.5*(A+B);
			std::cout<<"Moving to orbit a="<<particle_a_new<<std::endl<<std::endl;
//				}
//				else {
//					A = particle_a_new;
//					std::cout<<"Orbit a: "<<particle_a_new<<" is unstable."<<std::endl;
//					particle_a_new = 0.5*(A+B);
//					std::cout<<"Moving to orbit a="<<particle_a_new<<std::endl<<std::endl;
//				}
			}

//		 }

		count2 ++;

		 delete particle;
		 delete planet;
		 delete star;

		 delete[] r_planet;
		 delete[] r_particle;
		 delete[] v_planet;
		 delete[] v_particle;

		// particle_a_new += 0.5;
		// delete[] particle;
	}
  return 0;
}
