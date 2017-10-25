
/**
 * @file 	problem.c
 * @brief 	Gap Opening Simulation implemented with MPI
 * @author 	REBOUND: Hanno Rein <hanno@hanno-rein.de>
 * @author	Problem: Erika Nesvold <erika.r.nesvold@nasa.gov>
 * 
 * @section 	LICENSE
 * Copyright (c) 2011 Hanno Rein, Shangfei Liu
 *
 * This file is part of rebound.
 *
 * rebound is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * rebound is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "main.h"
#include "particle.h"
#include "boundaries.h"
#include "collision_resolve.h"
#include "output.h"
#include "communication_mpi.h"
#include "tree.h"
#include "tools.h"
#include "display.h"

extern int Nmax;

// Define the number of bins in the size distribution (and extrapolated size distirbution)
// of the superparticles
#define numbins 21 // Number of bins in the size distribution
#define extra 30   // Number of additional bins in the extrapolated size distribution

const int nowave = 1; // Flag to remove size distribution wave by extrapolation
const double logbinsize = 0.1; // bin size for size distribution
double powlogbinsize;	// calculated in problem_init
const int shift = 3;		// for building size distribution
const int rho = 3.0E3;		// Density (kg m^-3)

double sizebins[numbins];		// size bins
double bigbins[numbins+extra];	// extrapolated size bins
double mass[numbins];			// mass bins
double bigmass[numbins+extra];	// extrapolated mass bins
double Qd[numbins+extra];		// minimum projectile kinetic energy
double Ecol[numbins+extra][numbins+extra];  // collision energy
double Qsuper[numbins+extra];// threshold for supercatastrophic fragmentation
double Mlr[numbins+extra][numbins+extra]; // mass of largest fragment
double powbigbins[numbins+extra]; // used to normalize fragment distribution
double powsizebins[numbins];	// used to normalize fragment distribution
double fragdistsum[numbins];	// used to normalize fragment distribution

#include "smack.c"

void problem_init(int argc, char* argv[]){
	// Setup constants
	N_active = 2;
	N_tree_fixed = 2;
	boxsize = 10;	// AU
    root_nx = 2;
    root_ny = 2;
	dt = 0.1/0.159; // Timestep (years->code units)
	tmax = (1.0e6+1)*dt;
	
	init_box();
    
    // Change collision_resolve routing from default.
    collision_resolve = collision_resolve_single_fragment;

	// Initial conditions for star
	struct particle star;
	star.x 		= 0; star.y 	= 0; star.z	= 0;
	star.vx 	= 0; star.vy 	= 0; star.vz 	= 0;
	star.ax 	= 0; star.ay 	= 0; star.az 	= 0;
	star.m 		= 1.0;
	star.r		= 1.0*0.0046491;	// Rsun
	
    // Add the star to the particle array, but not to the tree
    Nmax += 128;
    particles = realloc(particles,sizeof(struct particle)*Nmax);
    particles[N] = star;
    N++;
    
    // Initial conditions for planet
    struct particle planet;
    float pa = 2.5;
    float pe = 0.0;
    float pi = 0.0;
    float pomega = 0.0;
    float pOMEGA = 0.0;
    float pM = 0.0;
    float pmass = 0.1*0.001;   // Mjup
    planet = tools_orbit2p(pa,pe,pi,pomega,pOMEGA,pM,star.m,pmass);
    planet.m = pmass;
    planet.r = 1.0*0.0004673;   // Rjup
    planet.number = 1;
    planet.ncol = 0;
    particles_add(planet);
	
	// Orbital parameters of disk
	float amin = 4.0;   // Semi-major axis
	float amax = 1.0;
	float emin = 0.0;   // Eccentricity
	float emax = 0.2;
	float imin = 0.0;   // Inclination
	float imax = emax/2.; 
	float Orange = 2*M_PI;	// Longitude of ascending node
	float orange = 2*M_PI;  // Argument of periapse
	float Mrange = 2*M_PI;  // Mean anomaly
	double a,e,i,OMEGA,omega,M;

    // Superparticle parameters
    int numSP = 200;    // Number of superparticles
    float rSP = pow(10.0,-(2.1));	// Superparticle radius (AU)
    float height = 0.107681;        // Height of disk (AU)
    float factor = 3*height/(4*rSP); // Superparticle filling factor
    float initOD = 1.0e-2; // Initial optical depth
    float Aeach = initOD*M_PI*rSP*rSP*1.5e11*1.5e11/factor; // Cross-sectional area of each superparticle (m^2)
    double p0 = 2.7; // Initial size distribution index
    float distsum = 0.0;         // Used for calculating initial superparticle size distribution
    float C;               // Coefficient of initial superparticle size distribution
    double sdist[numbins]; // Superparticle size distribution
    double spmass = 0.0;         // Superparticle mass
    
    // Create size bin arrays
    for (int i=0;i<numbins;i++){
        sizebins[i] = pow(10,i*logbinsize-shift); // Size bins (m)
    }
    if (nowave) { // Extrapolate to smaller sizes to remove size distribution wave
        for (int i=0; i<numbins+extra; i++) {
            bigbins[i] = pow(10,i*logbinsize-shift-extra*logbinsize); // Longer size bin array (m) for extrapolations
        }
    }
    powlogbinsize = 1.-pow(10.,-3.*logbinsize);
    
    // Collision energy pre-calculations
    const double Scon = 3.5e3;	// Strength regime coefficient
    const double Gcon = 3.0e-8;	// Gravity regime coefficient
    const double spow = -0.38;	// Strength regime exponent
    const double gpow = 1.36;	// Gravity regime exponent
    const double fKE = 1.0;		// Fraction of kinetic energy used
    {
        double bigvolume[numbins+extra];	// extrapolated volume bins
        for (int i=0; i<numbins+extra; i++) {
            if (i>=extra) {
                mass[i-extra] = rho*(4./3)*M_PI*pow(sizebins[i-extra]/2.,3); // small array of mass bins (kg)
            }
            bigvolume[i] = (4./3)*M_PI*pow(bigbins[i]/2.,3); // large array of volume bins (m^3)
            bigmass[i] = rho*bigvolume[i]; // large array of mass bins (kg)
            double Qs = Scon*pow(bigbins[i]*100./2,spow); // internal binding energy
            double Qg = Gcon*rho*pow(bigbins[i]*100./2,gpow); // gravitational binding energy
            Qd[i] = (1./fKE)*(Qs+Qg); // minimum projectile kinetic energy
            Qsuper[i] = -2.0*Qd[i]*(0.1-1); // supercatastrophic collision critical energy
        }
    }
    
    // Fragment distribution pre-caculations
    const double n0 = 2.8;			// Index of fragment distribution
    for (int j = 0; j < numbins+extra; j++) {
        powbigbins[j] = pow(bigbins[j],3.-n0);
    }
    for (int j = 0; j < numbins; j++) {
        powsizebins[j] = pow(sizebins[j],-1.*n0);
    }
    
    for (int lr = 0;  lr < numbins; lr++) {
        for (int k = 0; k <= (lr+extra); k++) {
            fragdistsum[lr] += powbigbins[k];
        }
    }
    
    // Fill size distribution
    for (int j=0; j<numbins; j++) {
        distsum += pow(bigbins[j],-1*p0)*pow(bigbins[j],2);
    }
    C = (4*Aeach)/(M_PI*distsum);	// Normalize size distribution
    for (int i=0; i<numbins; i++) {
        sdist[i] = C*pow(sizebins[i],-p0); // Fill size distribution
        spmass += (M_PI/6)*rho*pow(sizebins[i],3)*sdist[i]; // Keep track of total superparticle mass
    }
    spmass /= 1.98892e30; // Convert superparticles mass to solar units
	
    // Setup particle structures
#ifdef MPI
    for(int N=0;N<numSP/mpi_num;N++) {
#else
    for(int N=0;N<numSP;N++) {
#endif
        struct particle p;
        // Grab orbital parameters from uniform distributions
        a = tools_uniform(amin,amax);
        e = tools_uniform(emin,emax);
        i = tools_uniform(imin,imax);
        OMEGA = tools_uniform(0,Orange);
        omega = tools_uniform(0,orange);
        M = tools_uniform(0,Mrange);
        
        p = tools_orbit2p(a,e,i,omega,OMEGA,M,star.m,spmass); // Convert orbital parameters to Cartesian
        p.m = spmass; // Set superparticle mass
        p.r = rSP; // Set superparticle radius
        p.number = N + N_active; // Set superparticle number
#ifdef MPI
        p.number = N + N_active + mpi_id*numSP/mpi_num;
#endif
        p.ncol = 0; // Superparticle starts with zero collisions
        for (int i=0; i<numbins; i++) {
            p.sdist[i] = sdist[i]; // Set superparticle size distribution
        }
        particles_add(p); // Add superparticle to list of particles
        
#ifdef MPI
        communication_mpi_distribute_particles();
#endif
    }
    
    // IF YOU WANT TO RESTART YOUR MPI SIMULATION
    // Comment out the entire "Setup particle structures" for-loop above,
    // Comment out any lines adding planets above (i.e., particles_add(planet)),
    // And uncomment the following two lines:
    // input_smack_binary("restart.bin");
    // communication_mpi_distribute_particles();

}

void problem_inloop(){
}

void problem_output(){
	output_timing();
	if (output_check(dt*1.0e4)) {
	  output_smack_sizedist("sizedist.txt");
	  output_smack_orbits("orbits.txt");
	  output_smack_coordinates("coordinates.txt");
	}
}

void problem_finish(){
}
