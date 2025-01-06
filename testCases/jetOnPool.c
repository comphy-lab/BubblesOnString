/**
 * @file jetOnPool.c
 * @brief This file contains the simulation code for a VE liquid jet falling onto a VE pool. 
 * @author Vatsal Sanjay
 * @version 1.0
 * These are viscoelastic simulations!
 * @date Jan 6, 2025
*/

#include "axi.h"
#include "navier-stokes/centered.h"
#include "log-conform-viscoelastic-scalar-2D.h"

#define FILTERED // Smear density and viscosity jumps
#include "two-phaseVE.h"
// #include "two-phase.h"

#define logFile "log-jetOnPool_ViscoElastic.dat"

#include "navier-stokes/conserving.h"
#include "tension.h"
// #include "reduced.h"

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
// Oh -> Solvent Ohnesorge number
// Oha -> air Ohnesorge number
// Bo -> Bond number
// De -> Deborah number
// Ec -> Elasto-capillary number

double We, Re, MuR, Bo, Wi, El, tmax;
char nameOut[80], dumpFile[80];

// boundary conditions
f[left] = dirichlet(y > 1. + epsilon ? 0.0 : y < 1. - epsilon ? 1.0 : 0.5 * (1.0 + tanh((1e0 - R2(x,y,z)) / epsilon)));
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


int main(int argc, char const *argv[]) {


  // Values taken from the terminal
  L0 = 4e1; //atof(argv[1]);
  MAXlevel = 10; //atoi(argv[2]);
  We = 1e2; //atof(argv[3]);
  Re = 1e3; //atof(argv[4]);
  MuR = 2e-2;
  tmax = 1e2; //atof(argv[5]);

  Wi = 1e0; //atof(argv[6]);
  El = 1e-2; //atof(argv[7]);

  init_grid (1 << 8);

  // Create a folder named intermediate where all the simulation snapshots are stored.
  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);
  // Name of the restart file. See writingFiles event.
  sprintf (dumpFile, "restart");


  rho1 = 1., rho2 = 1e-3;
  mu1 = 1e0/Re, mu2 = MuR/Re;

  lambda1 = Wi, lambda2 = 0.;
  G1 = El, G2 = 0.;
  // G.x = Bo/We;
  f.sigma = 1.0/We;

  run();

}

event init (t = 0) {
  if (!restore (file = dumpFile)){
    refine(x < epsilon && level < MAXlevel);
    fraction (f, union(1-R2(x,y,z)-x/epsilon, x-4e0));
    foreach(){
      u.x[] = x < 1e0 ? f[]*2e0*(1-R2(x,y,z)): 0.;
      u.y[] = 0.0;
    }
  }
}

/**
## Adaptive Mesh Refinement
*/
event adapt(i++){
  scalar KAPPA[];
  curvature(f, KAPPA);
  adapt_wavelet ((scalar *){f, u.x, u.y, KAPPA, A11, A22, AThTh, A12},
      (double[]){fErr, VelErr, VelErr, KErr, AErr, AErr, AErr, AErr},
      MAXlevel, 4);
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
    fprintf(ferr, "Level %d, We %2.1e, Re %2.1e, MuR %2.1e, Wi %2.1e, El %2.1e\n", MAXlevel, We, Re, MuR, Wi, El);
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
      fprintf(ferr, "Level %d, We %2.1e, Re %2.1e, MuR %2.1e, Wi %2.1e, El %2.1e\n", MAXlevel, We, Re, MuR, Wi, El);
      fprintf(ferr, "i dt t ke\n");
      fprintf(fp, "Level %d, We %2.1e, Re %2.1e, MuR %2.1e, Wi %2.1e, El %2.1e\n", MAXlevel, We, Re, MuR, Wi, El);
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
