/************************************************************************************
 * This file contains all the functions needed for SMACK, including:
 * 
 * tools_linefit: fits a straight line to data
 * tools_orbit2p: Converts orbital elements to Cartesian coordinates
 * input_smack_binary: Reads a binary file into the particle struct
 * output_smack_coordinates: Outputs a text file containing Cartesian coordinates
 * output_smack_orbits: Outputs a text file containing orbital elements
 * output_smack_sizedist: Outputs a text file containing size distributions
 * output_smack_binary: Updates the state of the system in a binary file
 * output_smack_collisions: Outputs a text file containing the parameters of dust produced
 * collision_resolve_single_fragment: The main SMACK algorithm
 *
 * Author: Erika Nesvold
 *************************************************************************************/

/************************************************************************************
 * This function fits a set of numbers to a straight line y=ax+b, returning a and b
 * @author: Erika Nesvold
 * @param x[] x data
 * @param y[] y data
 * @param size Length of data
 * @return line
 *************************************************************************************/
extern long collisions_Nlog;
struct line {
	double slope;
	double intercept;
};

struct line tools_linefit(double x[], double y[], int size);

struct line tools_linefit(double x[], double y[], int size){
	
	struct line l;
	
	double totalx = 0.0;
	double totaly = 0.0;
	double totalxy = 0.0;
	double totalxx = 0.0;
	for (int i=0; i<size; i++) {
		const double log10x = log10(x[i]);
		const double log10y = log10(y[i]);
		totalx += log10x;
		totaly += log10y;
		totalxy += log10x*log10y;
		totalxx += log10x*log10x;
	}
	double avgx = totalx/size;
	double avgy = totaly/size;
	double avgxy = totalxy/size;
	double avgxx = totalxx/size;
	
	double numerator = avgxy - (avgx*avgy);
	double denominator = avgxx - (avgx*avgx);
	l.slope = numerator/denominator;
	l.intercept = avgy-(l.slope*avgx);
	
	return l;
}

/**************************************************************************************
 * This function calculates x,y,z,vx,vy,vz given a,e,i,omega,OMEGA,M
 * @author: Erika Nesvold
 * @param a semi-major axis
 * @param e eccentricity
 * @param i inclination
 * @param omega argument of periapse
 * @param OMEGA longitude of ascendning node
 * @param M mean anomaly
 * @param Ms mass of star
 * @param Mp mass of particle
 * @return p particle struct will Cartesian coords 
 **************************************************************************************/

struct particle tools_orbit2p(double a, double e, double i, double omega, double OMEGA, double M, double Ms, double Mp);

struct particle tools_orbit2p(double a, double e, double i, double omega, double OMEGA, double M, double Ms, double Mp){
	
	struct particle p;
	// Calculate eccentric anomaly
	double error = 1.0E-6;
	double E;
	if (M<M_PI) {
		E = M + e/2.;
	} else {
		E = M - e/2.;
	}
	double ratio = 1.0;
	while (fabs(ratio)>error) {
		ratio = (E-e*sin(E)-M)/(1-e*cos(E));
		E = E - ratio;
	}
	
	// Calculate x and y before rotation
	double x;
	double y;
	x = a*(cos(E)-e);
	y = a*sqrt(1-pow(e,2))*sin(E);
	
	// Rotate through angles
	double P1[3][3] = {{cos(omega),-1*sin(omega),0.},{sin(omega),cos(omega),0.},{0.,0.,1.}};
	double P2[3][3] = {{1.,0.,0.},{0.,cos(i),-1*sin(i)},{0.,sin(i),cos(i)}};
	double P3[3][3] = {{cos(OMEGA),-1*sin(OMEGA),0.},{sin(OMEGA),cos(OMEGA),0.},{0.,0.,1.}};
	double P4[3][3];
	double Q[3][3];
	for (int i=0; i<3; i++) {
		for (int j=0; j<3; j++) {
			P4[i][j]=P3[i][0]*P2[0][j]+P3[i][1]*P2[1][j]+P3[i][2]*P2[2][j];
		}
	}
	for (int i=0; i<3; i++) {
		for (int j=0; j<3; j++) {
			Q[i][j]=P4[i][0]*P1[0][j]+P4[i][1]*P1[1][j]+P4[i][2]*P1[2][j];
		}
	}
	double R[3][1];
	for (int i=0; i<3; i++) {
		R[i][0] = Q[i][0]*x+Q[i][1]*y+Q[i][2]*0.;
	}
	p.x = R[0][0];
	p.y = R[1][0];
	p.z = R[2][0];
	
	// Calculate velocities
	double mu = Ms+Mp;
	double xdot = -1*sqrt(mu*a)*sin(E)/sqrt(x*x+y*y);
	double ydot = sqrt(mu*a)*sqrt(1-e*e)*cos(E)/sqrt(x*x+y*y);
	double V[3][1];
	for (int i=0; i<3; i++) {
		V[i][0] = Q[i][0]*xdot+Q[i][1]*ydot+Q[i][2]*0.;
	}
	p.vx = V[0][0];
	p.vy = V[1][0];
	p.vz = V[2][0];
	
	p.m = 0;
	p.ax = 0;
	p.ay = 0;
	p.az = 0;
#ifndef COLLISIONS_NONE
	p.r = 0;
	p.lastcollision=0;
#endif
	
	return p;
}

/*********************************************************************************
 * Reads a binary file into the particle struct
 * @author: Erika Nesvold
 * @param filename Output filename
 *********************************************************************************/

void input_smack_binary(char* filename);

void input_smack_binary(char* filename){
#ifdef MPI
	char filename_mpi[1024];
	sprintf(filename_mpi,"%s_%d",filename,mpi_id);
	FILE* inf = fopen(filename_mpi,"rb"); 
#else // MPI
	FILE* inf = fopen(filename,"rb"); 
#endif // MPI
	long bytes = 0;
	int _N;
	bytes += fread(&_N,sizeof(int),1,inf);
	bytes += fread(&t,sizeof(double),1,inf);
#ifdef MPI
	fprintf(stderr,"Found %d particles in file '%s'. \n",_N,filename_mpi);
#else // MPI
	fprintf(stderr,"Found %d particles in file '%s'. \n",_N,filename);
#endif // MPI
	for (int i=0;i<_N;i++){
		struct particle p;
		bytes += fread(&p,sizeof(struct particle),1,inf);
#ifdef COLLISIONS_NONE
		if (p.m < 1.0) {
		  particles_add(p);
		}
#else
		if (p.number > 0) {
			particles_add(p);
		}
#endif
	}
	fclose(inf);
	fprintf(stderr,"%ld bytes read. Restarting at time t=%f\n",bytes,t);
}

/*********************************************************************************
 * Appends the coordinates of the superparticles to an ASCII file.
 * @author: Erika Nesvold
 * @param filename Output filename
 *********************************************************************************/
void output_smack_coordinates(char* filename);

void output_smack_coordinates(char* filename){
#ifdef MPI
	char filename_mpi[1024];
	sprintf(filename_mpi,"%s_%d",filename,mpi_id);
	FILE* of = fopen(filename_mpi,"a"); 
#else // MPI
	FILE* of = fopen(filename,"a"); 
#endif // MPI
	for (int i=0;i<N;i++){
		struct particle p = particles[i];
#ifdef COLLISIONS_NONE
		fprintf(of,"%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",t,p.x,p.y,p.z,p.vx,p.vy,p.vz,p.m);
#else
		fprintf(of,"%e\t%d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",t,(int)particles[i].number,p.x,p.y,p.z,p.vx,p.vy,p.vz,p.m);
#endif
	}
	fclose(of);
}

/*********************************************************************************
 * Appends the orbital elements of the superparticles to an ASCII file.
 * @author: Erika Nesvold
 * @param filename Output filename
 *********************************************************************************/
void output_smack_orbits(char* filename);

void output_smack_orbits(char* filename){
#ifdef MPI
	char filename_mpi[1024];
	sprintf(filename_mpi,"%s_%d",filename,mpi_id);
	FILE* of = fopen(filename_mpi,"a"); 
#else // MPI
	FILE* of = fopen(filename,"a"); 
#endif // MPI
	for (int i=0;i<N;i++){
		struct orbit o = tools_p2orbit(particles[i],particles[0]);
#ifdef COLLISIONS_NONE
		fprintf(of,"%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",t,o.a,o.e,o.inc,o.Omega,o.omega,o.l,particles[i].m);
#else
		fprintf(of,"%e\t%d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%d\n",t,(int)particles[i].number,o.a,o.e,o.inc,o.Omega,o.omega,o.l,particles[i].m,(int)particles[i].ncol);
#endif
	}
	fclose(of);
}

/*********************************************************************************
 * Appends the size distributions of the superparticles to an ASCII file.
 * @author: Erika Nesvold
 * @param filename Output filename
 *********************************************************************************/
void output_smack_sizedist(char* filename);

void output_smack_sizedist(char* filename){
#ifndef COLLISIONS_NONE
#ifdef MPI
	char filename_mpi[1024];
	sprintf(filename_mpi,"%s_%d",filename,mpi_id);
	FILE* of = fopen(filename_mpi,"a");
#else // MPI
	FILE* of = fopen(filename,"a");
#endif //MPI
	for (int i=N_active; i<N; i++) {
		fprintf(of,"%e\t%d\t%d\t",t,(int)particles[i].number,(int)particles[i].ncol);
		for (int j=0; j<numbins; j++) {
			fprintf(of,"%e\t",particles[i].sdist[j]);
		}
		fprintf(of,"\n");
	}
	fclose(of);
#endif
}

/*********************************************************************************
 * Updates the state of the system in a binary file.
 * @author: Erika Nesvold
 * @param filename Output filename
 *********************************************************************************/
void output_smack_binary(char* filename);

void output_smack_binary(char* filename){
#ifdef MPI
	char filename_mpi[1024];
	sprintf(filename_mpi,"%s_%d",filename,mpi_id);
	FILE* of = fopen(filename_mpi,"wb"); 
#else // MPI
	FILE* of = fopen(filename,"wb"); 
#endif // MPI
	if (of==NULL){
		printf("\n\nError while opening file '%s'.\n",filename);
		return;
	}
	int _N = N-1;
	fwrite(&_N,sizeof(int),1,of);
	fwrite(&t,sizeof(double),1,of);
	for (int i=1;i<N;i++){
		struct particle p = particles[i];
		fwrite(&(p),sizeof(struct particle),1,of);
	}
	fclose(of);
}

/*********************************************************************************
 * Outputs the parameters of dust produced in collisions
 * @author: Erika Nesvold
 * @param filename Output filename
 *********************************************************************************/

void output_smack_collisions(char* filename, float time, float x, float y, float z, float vx, float vy, float vz, float m);

void output_smack_collisions(char* filename, float time, float x, float y, float z, float vx, float vy, float vz, float m){
#ifdef MPI
  char filename_mpi[1024];
  sprintf(filename_mpi,"%s_%d",filename,mpi_id);
  FILE* of = fopen(filename_mpi,"a");
#else //MPI
  FILE* of = fopen(filename,"a");
#endif // MPI
  fprintf(of,"%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",t,x,y,z,vx,vy,vz,m);
  fclose(of);
}

/*********************************************************************************
 * Calculates the new size distributions and trajectories of the superparticles
 * @author: Erika Nesvold
 * @param c Collision structure pointing to two particle structs
 *********************************************************************************/

void collision_resolve_single_fragment(struct collision c);

double mfp = 0.0; /* Average mean free path */

void collision_resolve_single_fragment(struct collision c){	

#ifndef COLLISIONS_NONE
	if (particles[c.p1].number*particles[c.p2].number == 0) {
		return;	// Don't collide with star
	}
	if (particles[c.p1].number == 1 || particles[c.p2].number ==1){
	  // TODO: remove superparticle
	  return; // Don't collide with planet
	}
	struct particle p1 = particles[c.p1];
	struct particle p2;
#ifdef MPI
	int isloc = communication_mpi_rootbox_is_local(c.ri);
	if (isloc==1){
#endif // MPI
		p2 = particles[c.p2];
#ifdef MPI
	}else{
		int root_n_per_node = root_n/mpi_num;
		int proc_id = c.ri/root_n_per_node;
		p2 = particles_recv[proc_id][c.p2];
	}
#endif // MPI
	
	if (p1.lastcollision==t || p2.lastcollision==t) return;
	
	if (p1.number < N_active) {
		particles[c.p2].x += 2*boxsize;
		fprintf(stderr,"Planet Collision\n");
		return;
	}
	if (p2.number < N_active) {
		particles[c.p1].x += 2*boxsize;
		fprintf(stderr,"Planet Collision\n");
		return;
	}
	
	// Calculate new size distributions
	
	// Constants
	const double eqzero = 1.0e-2;
	const double mconv = 1.98892e30; // conversion from solar units to kg
	const double dconv = 1.5e11;	  // conversion from AU to m
	const double vconv = 29866.;	  // conversion to m/s
	const double vol1 = (4./3)*M_PI*pow(p1.r*dconv,3);	// Volume of swarm 1 (m^3)
	const double vol2 = (4./3)*M_PI*pow(p2.r*dconv,3);  // Volume of swarm 2 (m^3)
	const double eta = -1.5;                      // supercatastrophic largest remnant exponent
	const double eta18 = pow(1.8,eta);
	
	// Velocity calculations
	struct ghostbox gb = c.gb;
	const double x1 = p1.x + gb.shiftx; // x-coord of SP 1
	const double y1 = p1.y + gb.shifty; // y-coord of SP 1
	const double z1 = p1.z + gb.shiftz; // z-coord of SP 1
	const double x2 = p2.x; // x-coord of SP 2
	const double y2 = p2.y; // y-coord of SP 2
	const double z2 = p2.z; // z-coord of SP 2
	const double rp = p1.r+p2.r; // Sum of SP radii
	if (rp*rp < (x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1)) return; // if the SPs are not overlapping, return
	const double vx1 = p1.vx + gb.shiftvx; // x-vel of SP 1
	const double vy1 = p1.vy + gb.shiftvy; // y-vel of SP 1
	const double vz1 = p1.vz + gb.shiftvz; // z-vel of SP 1
	const double vx2 = p2.vx; // x-vel of SP 2
	const double vy2 = p2.vy; // y-vel of SP 2
	const double vz2 = p2.vz; // z-vel of SP 2
	if ((vx2-vx1)*(x2-x1)+(vy2-vy1)*(y2-y1)+(vz2-vz1)*(z2-z1)>0) return; // if the SPs are not approaching each other, return

	// Calculate center-of-mass
	const double m1 = p1.m; // mass of SP 1
	const double m2 = p2.m; // mass of SP 2
	const double cmx = (m1*vx1 + m2*vx2)/(m1 + m2); // x-coord of center of mass
	const double cmy = (m1*vy1 + m2*vy2)/(m1 + m2); // y-coord of center of mass
	const double cmz = (m1*vz1 + m2*vz2)/(m1 + m2); // z-coord of center of mass

	// Translate into center-of-mass frame
	const double vxa = vx1 - cmx; // x-vel of SP 1 in CoM frame
	const double vya = vy1 - cmy; // y-vel of SP 1 in CoM frame
	const double vza = vz1 - cmz; // z-vel of SP 1 in CoM frame
	const double vxb = vx2 - cmx; // x-vel of SP 2 in CoM frame
	const double vyb = vy2 - cmy; // y-vel of SP 2 in CoM frame
	const double vzb = vz2 - cmz; // z-vel of SP 2 in CoM frame
	const double va = sqrt(vxa*vxa + vya*vya + vza*vza); // vel of SP 1 in CoM frame
	const double vb = sqrt(vxb*vxb + vyb*vyb + vzb*vzb); // vel of SP 2 in CoM frame
	
	// Calculate pathlengths
	double relv = sqrt(pow(p1.vx-p2.vx,2)+pow(p1.vy-p2.vy,2)+pow(p1.vz-p2.vz,2));	// relative velocity of superparticles
	double pathlength1 = relv*(t-p1.lastcollision); // distance travelled by SP 1 since last collision 
	double pathlength2 = relv*(t-p2.lastcollision); // distance travelled by SP 2 since last collision 
	relv = relv*vconv; // relative velocity needs to be in m/s now
	double newrelv = relv; // relative velocity for each loop
	pathlength1 = pathlength1*dconv; // needs to be in m now
	pathlength2 = pathlength2*dconv; // needs to be in m now
	
	// Correct pathlengths for superparticles from outside the disk
	if (pathlength1 > 5.*mfp && collisions_Nlog > 10) {
		pathlength1 = 5.*mfp;
		// fprintf(stderr,"Corrected pathlength1\n");
	}
	if (pathlength2 > 5.*mfp && collisions_Nlog > 10) {
		pathlength2 = 5.*mfp;
		//fprintf(stderr,"Corrected pathlength2\n");
	}
	
	// Initialize arrays for loss, survivors, and fragments
	double survivors1[numbins]; // survivors in SP 1
	double survivors2[numbins]; // survivors in SP 2
	
	// Check for maximum optical depth
	double od1[numbins+extra]; // total optical depth for each bin in SP 1 due to SP 2
	double od2[numbins+extra]; // total optical depth for each bin in SP 2 due to SP 1
	double od1arr[numbins+extra][numbins+extra];	// [target][projectile]
	double od2arr[numbins+extra][numbins+extra];	// Prob of target in # hit by projectile in other
	double Eloss1 = 0.0; // total energy loss in SP 1
	double Eloss2 = 0.0; // total energy loss in SP 2
	double newdist1[numbins]; // new size distribution for SP 1
	double newdist2[numbins]; // new size distribution for SP 2
	// Start with new size distributions equal to initial size distributions
	for (int i=0; i<numbins; i++) {
		newdist1[i] = p1.sdist[i];
		newdist2[i] = p2.sdist[i];
	}
	
	double newpath1 = 0.0; // readjusted path length for SP 1
	double newpath2 = 0.0; // readjusted path length for SP 2
	double rempath1 = pathlength1; // remaining pathlength for SP 1
	double rempath2 = pathlength2; // remaining pathlength for SP 2
	int last = 0; // Flag for last loop
	
	struct line line1; // fit for extrapolated size distribution for SP 1
	struct line line2; // fit for extrapolated size distribution for SP 2
	double bigdist1[numbins+extra]; // extrapolated size distribution for SP 1
	double bigdist2[numbins+extra]; // extrapolated size distribution for SP 2
	
	double newm1,newm2; // new masses of SP 1 and SP 2
	double totm1,totm2; // new mass + lost dust mass for SP 1 and SP 2
	double KE; // total kinetic energy of SP 1 and SP 2
	double vc,vd,vxc,vyc,vzc,vxd,vyd,vzd; // new velocities for SP 1 and SP 2 in CoM frame
	double frac1,frac2; // fractional change in speed for SP 1 and SP 2
	
	int extrapolate1 = 0; // flag for whether to extrapolate SP 1
	int extrapolate2 = 0; // flag for whether to extrapolate SP 2
	int fitlength = numbins;    // number of bins to use for extrapolation
	double sdist1[fitlength];    // subset of SP 1 size dist to use for extrapolation
	double sdist2[fitlength];    // subset of SP 2 size dist to use for extrapolation
	double sizebins1[fitlength]; // subset of size bins to use for extrapolation
	double sizebins2[fitlength]; // subset of size bins to use for extrapolation
	
	double fragmentdist[numbins][numbins+extra][numbins];	// fragment distribution [target] [projectile] [fragments]
	
	int empty1, empty2, empty;
	
	/////////////////////////////////////////////////////////////////////////////////
	// Main Loop
	/////////////////////////////////////////////////////////////////////////////////
	
	int numloop = 0; // number of loops
	do {  
		numloop++; // track number of loops
		
		///////////////////////////////////////////////////////////////////////////////
		// Extrapolate size distributions
		// Check for empty bins
		if (nowave) {
			extrapolate1 = 1;
			extrapolate2 = 1;
			for (int i=0; i<fitlength; i++){
				sdist1[i] = newdist1[i];
				sdist2[i] = newdist2[i];
				sizebins1[i] = sizebins[i];
				sizebins2[i] = sizebins[i];
				if (sdist1[i] <= 0) {
					extrapolate1 = 0; // If one of these bins is empty, don't extrapolate SP 1
				}
				if (sdist2[i] <= 0) {
					extrapolate2 = 0; // If one of these bins is empty, don't extrapolate SP 2
				}
			}
		}
		
		// Extrapolate size distribution for SP 1 if necessary
		if (extrapolate1) {
			line1 = tools_linefit(sizebins1,sdist1,fitlength); // Fit to size distribution
			if (isnan(line1.slope)) {
			  fprintf(stderr,"\nFit line with slope = %f with fitlength %d to:\n",line1.slope,fitlength);
			  fprintf(stderr,"sizebins:\n");
			  for (int j=0; j<numbins; j++){
			    fprintf(stderr,"%e\t",sizebins1[j]);
			  }
			  fprintf(stderr,"\nsdist:\n");
			  for (int j=0; j<numbins; j++){
			    fprintf(stderr,"%e\t",sdist1[j]);
			  }
			  fprintf(stderr,"\n");
			}
			for (int i=0; i<extra; i++) {
				bigdist1[i] = pow(10.,line1.intercept)*pow(bigbins[i],line1.slope); // Build extrapolated distribution
			}
			// Error checking
			if (isnan(bigdist1[0])) {
				fprintf(stderr,"\nCollision between %.0f and %.0f yielded nan in %.0f -- bigdist\n",particles[c.p1].number,particles[c.p2].number,particles[c.p1].number);
				for (int j=0; j<numbins; j++) {
					fprintf(stderr,"%e ",newdist1[j]);
				}
				fprintf(stderr,"\nslope = %f\n",line1.slope);
				exit(0);
			}
		} else {
			for (int i=0; i<extra; i++) {
				bigdist1[i] = 0.0; // Don't build extrapolated distribution if not necessary
			}
		}
		
		
		// Extrapolate size distribution for SP 2 if necessary
		if (extrapolate2) {
			line2 = tools_linefit(sizebins2,sdist2,fitlength); // Fit to size distribution
			if (isnan(line2.slope)) {
			  fprintf(stderr,"\nFit line with slope = %f with fitlength %d to:\n",line2.slope,fitlength);
			  fprintf(stderr,"sizebins:\n");
			  for (int j=0; j<numbins; j++){
			    fprintf(stderr,"%e\t",sizebins2[j]);
			  }
			  fprintf(stderr,"\nsdist:\n");
			  for (int j=0; j<numbins; j++){
			    fprintf(stderr,"%e\t",sdist2[j]);
			  }
			  fprintf(stderr,"\n");
			}
			for (int i=0; i<extra; i++) {
				bigdist2[i] = pow(10.,line2.intercept)*pow(bigbins[i],line2.slope); // Build extrapolated distribution
			}
			// Error checking
			if (isnan(bigdist2[0])) {
				fprintf(stderr,"\nCollision between %.0f and %.0f yielded nan in %.0f -- bigdist\n",particles[c.p1].number,particles[c.p2].number,particles[c.p2].number);
				for (int j=0; j<numbins; j++) {
					fprintf(stderr,"%e ",newdist2[j]);
				}
				fprintf(stderr,"\nslope = %f\n",line2.slope);
				exit(0);
			}
			
		} else {
			for (int i=0; i<extra; i++) {
				bigdist2[i] = 0.0; // Don't build extrapolated distribution if not necessary
			}
		}
		// Fill in the rest of the extrapolated distribution with the original distribution
		for (int i=extra; i<numbins+extra; i++) {
			bigdist1[i] = p1.sdist[i-extra]; 
			bigdist2[i] = p2.sdist[i-extra];
		}
		
		////////////////////////////////////////////////////////////////////////////
		// Calculate optical depth
		const int loopcount = numbins+extra;
		//#pragma omp parallel for 
		for (int i=0; i<loopcount; i++) {
			od1[i] = 0.0;
			od2[i] = 0.0;
			for (int j=0; j<numbins+extra; j++) {
				Ecol[i][j] = 0.5*bigmass[i]*bigmass[j]*newrelv*newrelv/(bigmass[i]+bigmass[j]);
				double Qcol = Ecol[i][j]/bigmass[i];  // collision energy/mass of target
				// Calculate size of largest remnant from collisional energy
				if (Qcol >= Qsuper[i]) {
					Mlr[i][j] = bigmass[i]*(0.1/eta18)*pow(Qcol/Qd[i],eta);
				} else {
					Mlr[i][j] = bigmass[i]*(-0.5*((Qcol/Qd[i])-1)+0.5);
				}
				if (Mlr[i][j] < 0.0) {
					fprintf(stderr,"Error -- Mlr is less than zero \nEnding program\n");
					exit(0);
				}
				//if (Mlr[i][j]/bigmass[i] < powlogbinsize) {
				double term = (bigbins[j]+bigbins[i]);
				double term1 = term*term*(M_PI/4.)/vol1; // Why do I divide by volume here and not in the following lines?
				double term2 = term*term*(M_PI/4.)/vol2; // Don't question it. It had to be this way.
                od1arr[i][j] = bigdist2[j]*term1*rempath1;
                od2arr[i][j] = bigdist1[j]*term2*rempath2;
				od1[i] += od1arr[i][j];
				od2[i] += od2arr[i][j];
			}
		}
		double maxod = 0.0; // max optical depth
		for (int i=0; i<numbins+extra; i++) {
			// Calculate max optical depth
			if (od1[i] >= maxod && i >= extra) {
				maxod = od1[i];
			}
			if (od2[i] >= maxod && i >= extra) {
				maxod = od2[i];
			}
		}
		
		// If max optical depth is greater than 1, only go that fraction of the pathlength
		if (maxod > 1) {
			newpath1 = 1.0*rempath1/maxod; // readjust pathlength
			newpath2 = 1.0*rempath2/maxod; // readjust pathlength
			// Readjust optical depths
			for (int i=0; i<numbins+extra; i++) {
				for (int j=0; j<numbins+extra; j++) {
					od1arr[i][j] = 1.0*od1arr[i][j]/maxod;
					od2arr[i][j] = 1.0*od2arr[i][j]/maxod;
				}
				od1[i] = 1.0*od1[i]/maxod;
				od2[i] = 1.0*od2[i]/maxod;
			}
			// Otherwise, finish the loop
		} else {
			newpath1 = rempath1;
			newpath2 = rempath2;
			last = 1;
		}
		
		///////////////////////////////////////////////////////////////////////////////////////
		// Calculate fragment distributions
		// Loop through each target
		
		//#pragma omp parallel for 
		for (int i = 0; i < numbins; i++) {
			// Loop through each projectile
			for (int j = 0; j < numbins+extra; j++) {
				int lr = 0; // index of largest fragment
				// Find index of largest fragment
				while (mass[lr] < Mlr[i+extra][j]) {
					lr ++;
				}
				if (lr >= numbins) {
					fprintf(stderr,"Warning! Largest remnant outside array.\n");
					fprintf(stderr,"Mlr = %e\n",Mlr[i+extra][j]);
					fprintf(stderr,"bigmass = %e\n",bigmass[i+extra]);
					fprintf(stderr,"Diff = %e\n",Mlr[i+extra][j]-bigmass[i+extra]);
					exit(0);
				}
				// Normalize fragment distribution equal the mass of the target
				double B = 6.*mass[i]/(M_PI*rho*fragdistsum[lr]);
				for (int k = 0; k < numbins; k++) {
					if (k <= lr) {
						fragmentdist[i][j][k] = B*powsizebins[k];
					} else {
						fragmentdist[i][j][k] = 0.0;
					}
				}
			}
		}
		
		double fragments1[numbins]; // fragments from SP 1
		double fragments2[numbins]; // frgments from SP 2
		// Initialize fragment arrays
		for (int i=0; i<numbins; i++) {
			fragments1[i] = 0.0;
			fragments2[i] = 0.0;
		}
		
		/////////////////////////////////////////////////////////////////////////////
		// Calculate survivors and fragments
		for (int i=0; i<numbins; i++) {
			double EL1 = 0.0;  // energy loss per bin in SP 1
			double EL2 = 0.0;  // energy loss per bin in SP 2
			// Calculate loss from each bin
			double loss1 = od1[i+extra]*newdist1[i];
			double loss2 = od2[i+extra]*newdist2[i];
			for (int j=0; j<numbins+extra; j++) {
				// Add to energy loss
				// Each planetesimal in each collision loses half the collisional energy
				EL1 += 0.5*Ecol[i][j]*loss1;
				EL2 += 0.5*Ecol[i][j]*loss2;
			}
			// Total energy loss
			Eloss1 += EL1;
			Eloss2 += EL2;
			// Number of survivors
			survivors1[i] = newdist1[i] - loss1;
			survivors2[i] = newdist2[i] - loss2;
			// Add fragments to smaller bins
			for (int j = 0; j < numbins+extra; j++) {
				double prefac1 = od1arr[i+extra][j]*newdist1[i];
				double prefac2 = od2arr[i+extra][j]*newdist2[i];
				for (int k = 0; k < i; k++) {
					fragments1[k] += prefac1*fragmentdist[i][j][k];
					fragments2[k] += prefac2*fragmentdist[i][j][k];
				}
			}
		}
		// Calculate new size distributions 
		for (int i=0; i<numbins; i++) {
			newdist1[i] = survivors1[i]+fragments2[i];
			newdist2[i] = survivors2[i]+fragments1[i];
		}
		// Calculate remaining pathlength
		rempath1 = rempath1 - newpath1; 
		rempath2 = rempath2 - newpath2;
		
		// Calculate new total mass of each SP
		newm1 = 0.0;
		newm2 = 0.0;
		for (int i=0; i<numbins; i++) {
			newm1 += mass[i]*newdist1[i]/mconv; // now in solar masses
			newm2 += mass[i]*newdist2[i]/mconv;
		}
		
		empty1 = 1;
		for (int i=0; i<numbins; i++) {
			if (newdist1[i] > eqzero) {
				empty1 = 0;
			} 
		}
		empty2 = 1;
		for (int i=0; i<numbins; i++) {
			if (newdist2[i] > eqzero) {
				empty2 = 0;
			} 
		}
		empty = 0;
		if (empty1 || empty2) {
			empty = 1;
		}
		
	} while (! last && !empty && numloop < 1001);  // While the optical depth is still greater than 1 and there is pathlength left, restart loop

	/////////////////////////////////////
	/////////////////////////////////////
	////////////END OF LOOP//////////////
	/////////////////////////////////////
	/////////////////////////////////////
	
	if (empty1) {
		particles[c.p1].x += 2.*boxsize;
		fprintf(stderr,"Superparticle %d is empty\n",(int)particles[c.p1].number);
	}
	if (empty2) {
		particles[c.p1].x += 2.*boxsize;
		fprintf(stderr,"Superparticle %d is empty\n",(int)particles[c.p2].number);
	}
	
	// Calculate new kinetic energy (in CoM frame for the two swarms)
	KE = 0.5*m1*va*va*mconv*vconv*vconv + 0.5*m2*vb*vb*mconv*vconv*vconv - (Eloss1+Eloss2); // Joules
	// Error checking
	if (KE < 0.0) {
		fprintf(stderr,"\nERROR -- More energy used in collisions than available\n");
		fprintf(stderr,"\nEloss1 = %e\n",Eloss1);
		fprintf(stderr,"Eloss2 = %e\n",Eloss2);	
		fprintf(stderr,"Einitial = %e\n",(0.5*m1*va*va+0.5*m2*vb*vb)*mconv*vconv*vconv);
		fprintf(stderr,"KE = %e\n",KE);
		exit(0);
		KE = 0.0;
	}
	
	// Adjust velocities to conserve energy and momentum
	totm1 = (newm1/(newm1+newm2))*(m1+m2);
	totm2 = (newm2/(newm1+newm2))*(m1+m2);
	vc = sqrt(totm2*2*KE/(totm1*(totm1 + totm2)*mconv)); // m/s
	vd = (totm1/totm2)*vc; // m/s
	frac1 = va/vc;
	frac2 = vb/vd;
	if (frac1*frac2 == 0) {
		vxc = vxa;
		vyc = vya;
		vzc = vza;
		vxd = vxb;
		vyd = vyb;
		vzd = vzb;
	} else {
		vxc = (vxa/frac1)/vconv; // now all in REBOUND velocity units
		vyc = (vya/frac1)/vconv;
		vzc = (vza/frac1)/vconv;
		vxd = (vxb/frac2)/vconv;
		vyd = (vyb/frac2)/vconv;
		vzd = (vzb/frac2)/vconv;
	}
	
	// Calculate new relative velocity 
	newrelv = sqrt(pow(vxc-vxd,2)+pow(vyc-vyd,2)+pow(vzc-vzd,2));
	newrelv *= vconv;	  
	
	if (isnan(newrelv) || numloop < 0) {
		fprintf(stderr,"va = [%e,%e,%e]\n",vxa,vya,vza);
		fprintf(stderr,"vb = [%e,%e,%e]\n",vxb,vyb,vzb);
		fprintf(stderr,"pathlength1 = %e\n",pathlength1);
		fprintf(stderr,"pathlength2 = %e\n",pathlength2);
		fprintf(stderr,"mfp*7 = %e\n",mfp*7);
		fprintf(stderr,"Collision between %d (%d) and %d (%d)\n",(int)particles[c.p1].number,(int)particles[c.p1].ncol,(int)particles[c.p2].number,(int)particles[c.p2].ncol);
		fprintf(stderr,"Olddist1:\n");
		for (int i=0; i<numbins; i++) {
			fprintf(stderr,"%e,\t",p1.sdist[i]);
		}
		fprintf(stderr,"\n");
		fprintf(stderr,"Olddist2:\n");
		for (int i=0; i<numbins; i++) {
			fprintf(stderr,"%e,\t",p2.sdist[i]);
		}
		fprintf(stderr,"\n");
		fprintf(stderr,"Newdist1:\n");
		for (int i=0; i<numbins; i++) {
			fprintf(stderr,"%e\t",newdist1[i]);
		}
		fprintf(stderr,"\n");
		fprintf(stderr,"Newdist2:\n");
		for (int i=0; i<numbins; i++) {
			fprintf(stderr,"%e\t",newdist2[i]);
		}
		fprintf(stderr,"\n");
		fprintf(stderr,"newrelv = %e\n",newrelv);
		fprintf(stderr,"numloop = %d\n",numloop);
		exit(0);
	}
	
	// Translate back to original coordinate system
	double vx3 = vxc + cmx; // x-vel of SP 1 in REBOUND frame
	double vy3 = vyc + cmy; // y-vel of SP 1 in REBOUND frame
	double vz3 = vzc + cmz; // z-vel of SP 1 in REBOUND frame
	double vx4 = vxd + cmx; // x-vel of SP 2 in REBOUND frame
	double vy4 = vyd + cmy; // y-vel of SP 2 in REBOUND frame
	double vz4 = vzd + cmz; // z-vel of SP 2 in REBOUND frame
	
	// Keep track of number of collisions and mean free path
	collisions_Nlog++;
	
	// Don't count corrected pathlengths in mfp
	if (pathlength1 == 5*mfp) {
		pathlength1 = mfp;
	}
	if (pathlength2 == 5*mfp) {
		pathlength2 = mfp;
	}
	mfp = ((mfp*(collisions_Nlog-1)*2)+pathlength1+pathlength2)/(collisions_Nlog*2);
	
	// Apply the changes to the particles.
#ifdef MPI
	if (isloc==1){
#endif // MPI
		particles[c.p2].vx = vx4;
		particles[c.p2].vy = vy4;
		particles[c.p2].vz = vz4;
		particles[c.p2].m = newm2;
		particles[c.p2].lastcollision = t;
		for (int i=0; i<numbins; i++) {
			particles[c.p2].sdist[i] = newdist2[i];
		}
		particles[c.p2].ncol += 1;
#ifdef MPI
	}
#endif // MPI
	particles[c.p1].vx = vx3;
	particles[c.p1].vy = vy3;
	particles[c.p1].vz = vz3;
	particles[c.p1].m = newm1;
	particles[c.p1].lastcollision = t; 
	for (int i=0; i<numbins; i++) {
		particles[c.p1].sdist[i] = newdist1[i];
	}
	particles[c.p1].ncol += 1;
#endif // COLLISIONS_NONE
}

