/**
 * @file jetOnPool.c
 * @brief This file contains the simulation code for a VE liquid jet falling onto a VE pool. 
 * @author Vatsal Sanjay
 * @version 1.0
 * These are viscoelastic simulations!
 * @date Jan 6, 2025
*/

#include "embed.h"
#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#define logFile "log-jetOnly-test.dat"
#include "tension.h"

#define tsnap (1e-2)

// Error tolerancs
#define fErr (1e-3)                                 // error tolerance in f1 VOF
#define KErr (1e-4)                                 // error tolerance in VoF curvature calculated using heigh function method (see adapt event)
#define VelErr (1e-2)                               // error tolerances in velocity -- Use 1e-2 for low Oh and 1e-3 to 5e-3 for high Oh/moderate to high J
#define AErr (1e-3)                                 // error tolerance in VoF curvature calculated using heigh function method (see adapt event)

#define epsilon (8e-2)
#define R2(x,y,z) (sq(y) + sq(z))
// #define R2NonLam(x,y,z) pow(1 - sqrt(R2(x,y,z)),1/7)*8/7


int MAXlevel;
// We -> Weber number
// Re_s -> Solvent Reynolds number
// muR -> air-solvent viscosity ratio

double We, Re_s, muR, tmax;
char nameOut[80], dumpFile[80];

// boundary conditions
f[left] = neumann(0.0);
u.n[left] = dirichlet(f[]*2e0*(1-R2(x,y,z)));
u.t[left] = dirichlet(0.0);


u.n[right] = neumann(0.0);
u.t[right] = neumann(0.0);
uf.n[right] = neumann(0.0);
uf.t[right] = neumann(0.0);
p[right]   = dirichlet(0.0);

u.n[top] = neumann(0.0);
u.t[top] = neumann(0.0);
uf.n[top] = neumann(0.0);
uf.t[top] = neumann(0.0);
p[top]   = dirichlet(0.0);

vertex scalar phi[];

int main(int argc, char const *argv[]) {


  // Values taken from the terminal
  L0 = 1e1; //atof(argv[1]);
  MAXlevel = 10; //atoi(argv[2]);
  We = 1e2; //atof(argv[3]);
  Re_s = 1e3; //atof(argv[4]);
  muR = 2e-2;

  tmax = 1e2; //atof(argv[5]);

  init_grid (1 << 8);

  // Create a folder named intermediate where all the simulation snapshots are stored.
  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);
  // Name of the restart file. See writingFiles event.
  sprintf (dumpFile, "restart");


  rho1 = 1., rho2 = 1e-3;
  mu1 = 1e0/Re_s, mu2 = muR/Re_s;

  f.sigma = 1.0/We;

  run();

}

event init (t = 0) {
  if (!restore (file = dumpFile)){
    refine(R2(x,y,z) < 1.5e0 && level < MAXlevel);
    fraction (f, intersection(1-R2(x,y,z), 1e0-x));
    foreach(){
      u.x[] = x < 1e0 ? f[]*2e0*(1-R2(x,y,z)): 0.;
      u.y[] = 0.0;
    }
  }
  #if EMBED
  foreach_vertex(){
    phi[] = -(y - 4e0);
  }
  fractions (phi, cs, fs);
  #endif
}

/**
## Adaptive Mesh Refinement
*/
event adapt(i++){
  scalar KAPPA[];
  curvature(f, KAPPA);
  adapt_wavelet ((scalar *){f, u.x, u.y, KAPPA},
      (double[]){fErr, VelErr, VelErr, KErr},
      MAXlevel, 4);
  
  #if EMBED
  foreach_vertex(){
    phi[] = -(y - 4e0);
  }
  fractions (phi, cs, fs);
  #endif
}

/**
## Dumping snapshots
*/
event writingFiles (t = 0; t += tsnap; t <= tmax) {
  dump (file = dumpFile);
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump(file=nameOut);
}

/**
## Ending Simulation
*/
event end (t = end) {
  if (pid() == 0)
    fprintf(ferr, "Level %d, We %2.1e, Re_s %2.1e, muR %2.1e\n", MAXlevel, We, Re_s, muR);
}

/**
## Log writing
*/
event logWriting (i++) {

  double ke = 0.;
  foreach (reduction(+:ke)){
    ke += (2*pi*y)*(0.5*rho(f[])*(sq(u.x[]) + sq(u.y[])))*sq(Delta);
  }

  static FILE * fp;
  if (pid() == 0) {
    const char* mode = (i == 0) ? "w" : "a";
    fp = fopen(logFile, mode);
    if (fp == NULL) {
      fprintf(ferr, "Error opening log file\n");
      return 1;
    }

    if (i == 0) {
      fprintf(ferr, "Level %d, We %2.1e, Re_s %2.1e, muR %2.1e\n", MAXlevel, We, Re_s, muR);
      fprintf(ferr, "i dt t ke\n");
      fprintf(fp, "Level %d, We %2.1e, Re_s %2.1e, muR %2.1e\n", MAXlevel, We, Re_s, muR);
      fprintf(fp, "i dt t ke\n");
    }

    fprintf(fp, "%d %g %g %g\n", i, dt, t, ke);
    fprintf(ferr, "%d %g %g %g\n", i, dt, t, ke);

    fflush(fp);
    fclose(fp);
  }

  if(ke < -1e-10) return 1;

  if (i > 1e1 && pid() == 0) {
    if (ke > 1e2 || ke < 1e-8) {
      const char* message = (ke > 1e2) ? 
        "The kinetic energy blew up. Stopping simulation\n" : 
        "kinetic energy too small now! Stopping!\n";
      
      fprintf(ferr, "%s", message);
      
      fp = fopen(logFile, "a");
      fprintf(fp, "%s", message);
      fflush(fp);
      fclose(fp);
      
      dump(file=dumpFile);
      return 1;
    }
  }

}
