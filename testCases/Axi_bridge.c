/**
 Viscoelastic bridge

A viscoelastic pinch-off characterizes by the "beads on a string"
structure.  In the context of the viscoelastic bridges, Bhat et
al. (2010) conclude that this structure is due to inertia although
the delay in the pinch-off caused by the viscoelastic stresses allows
to this structure to persist. */


#include "axi.h"
#include "navier-stokes/centered.h"
// #define FILTERED 
#include "two-phase.h"
//#include "log-conform.h"
#include "tension.h"
#include "log-conform.h"
// #include "./fene-p2.h"
// #include "./log-conform_loc.h"
//#include "tag.h"

#define RHO_r 0.01
#define OHN 0.2
#define BETA 0.85 //(t > 5. ? 1.-0.15*(1.-exp(-10.*(t-5.))) : 1.)
#define MU_r 0.01
#define DEB 1.
#define LAM 5.
#define ws 15000.
#define Ls 4.
#define ytw (1+ws)
#define ybw  1.
#define domain_size 100. //16.*Ls
#define Fr 0.38
#define Amp 0.4
#define rho1b 1.
#define rho2b rho1b*RHO_r
#define rho_s rho1b/10.
#define mu_s mu1*1.e2
//#define FILTERED 1.
//#define dp_l (t > 0. ? 6.4*(1.+Amp*sin(2.*pi*Fr*(t-0.))) : 6.4)
#define dp_l 6.4*(1.+Amp*sin(2.*pi*Fr*t))


double rc_nozzle = 0.1;

double r_curv_nozzle_1 = 0.1;
double r_curv_nozzle_2 = 0.05;


//double dp_l = 10.; //*(1.+0.2*sin(2.*pi*Fr*t))

int LEVEL=10; //13;

scalar ffilt[],fsfilt[],Sxxfilt[];
vertex scalar phi[], phis[];
scalar csibm[];
scalar dtaupxx[],dtaupyy[],dtaupyx[],dtauqq[];

scalar mask_dtau[], mask_dtau_loc[];
scalar taupxx[];
scalar dp[];
face vector fsibm[];

scalar full_k[];

//face vector muv[];
scalar lambdav[], mupv[];

#define  Uavg 5 //4.//*(1.+Amp*sin(2.*pi*Fr*t))

p[left]=neumann(y < ytw ? dp_l : 0.);
pf[left]=neumann(y < ytw ? dp_l : 0.);
csibm[left] = neumann(0.);
u.n[left] = neumann(0.); //dirichlet(2.*Uavg*(1.-sq(y))*f[]*(1.-csibm[])); //x
u.t[left] = neumann(0.);
f[left] = neumann(0.);//(y <1.5); //(1.+0.*Delta) ? 1. : 0.;


p[right] = dirichlet(0.);
pf[right] = dirichlet(0.);
u.t[right]   = neumann(0.);
u.n[right]   = neumann(0.);
f[right]  = neumann(0.);
csibm[right] = neumann(0.);

uf.n[bottom] = 0.;

uf.n[top] = 0.; //dirichlet(0.);
uf.t[top] = 0.;
u.t[top] = dirichlet(0.);
p[top] = neumann(0.);

#define mup(f)  ((1.-BETA)*OHN*clamp(f,0,1))
#define lambda(f)  (DEB*clamp(f,0,1))
//# define rho3(f,csibm) (clamp(1.-csibm,0,1)*((clamp(f,0.,1.)*(rho1b - rho2b) + rho2b)-rho_s) + rho_s)
//# define mu3(f,csibm) (clamp(1.-csibm,0,1)*((clamp(f,0.,1.)*(mu1 - mu2) + mu2)-mu_s) + mu_s)


void create_ibm (scalar csibm) {

    foreach_vertex(){
      if (y > ybw + r_curv_nozzle_2) {
        phis[] = -(x-Ls);
      }
      else {
          if (x < Ls - r_curv_nozzle_1) {
              phis[] = (y - (ybw));
          }
       else {
              phis[] = 1. - sqrt(sq(x-(Ls - r_curv_nozzle_1))/sq(r_curv_nozzle_1)+sq(y-(ybw+r_curv_nozzle_2))/sq(r_curv_nozzle_2)) ;
          }
      }
    }
 
    foreach_vertex()
    if (x > Ls) {
      if (fabs(y-8.) < fabs(phis[]) || y > 8.) {
        phis[] = (y - 8.);
      }
    }

    boundary ({phis});
    fractions (phis, csibm, fsibm);
    // fractions_cleanup (csibm,fsibm);
    // boundary ({csibm,fsibm});

     csibm.refine = csibm.prolongation = fraction_refine;
     boundary ({csibm});
}

void filter_vof (scalar floc, int niter, scalar ftemp) {
  
  foreach(){
      ftemp[] = floc[];
  }
  boundary ({ftemp});

  for (int ii=1;ii<=niter;ii++) {
    scalar ftemp2 = ftemp;
    foreach(){
      // ftemp2[] = (4.*ftemp[] +
      //   2.*(ftemp[0,1] + ftemp[0,-1] + ftemp[1,0] + ftemp[-1,0]) +
      //   ftemp[-1,-1] + ftemp[1,-1] + ftemp[1,1] + ftemp[-1,1])/16.; //0.5*ftemp[] + 0.125*(ftemp[-1,0] + ftemp[1,0] + ftemp[0,1] + ftemp[0,-1]);
    ftemp2[] = (6.*ftemp[] + 0.5*(ftemp[0,1] + ftemp[0,-1] + ftemp[1,0] + ftemp[-1,0]))/8.; 
    }
    ftemp = ftemp2;

    ftemp[left] = neumann(0.);
    ftemp.refine = ftemp.prolongation = refine_bilinear;
    boundary ({ftemp});
  }

}

int main() {

 DT = 1e-5;
 TOLERANCE = 1e-6;
  size (domain_size);

  init_grid (1 << 8); //-2);
  // L2 = sq(30.);

  rho1 = 1.; rho2 =rho1*RHO_r;
  mu1 = BETA*OHN; mu2 = mu1*MU_r;
  lambda = lambdav;
  mup = mupv;
  f.sigma = 1.;
  run();
}


event init (t = 0) {

  if (!restore (file = "restart")) {

    //refine(x < 1.2*Ls && y < 1.2*ytw && level < LEVEL);
   refine(x < 1.2*Ls && level < LEVEL);

    DT = 1e-5;
    TOLERANCE = 1e-6;

    foreach_vertex() {
         phi[] = -(x - 3.95);
    }
    fractions (phi,f);
    boundary({f});

  foreach(){
  //    tau_p.x.x[] = 20*sq(y)*f[];
  //    tau_p.x.y[] = 2*sq(y)*f[];
  //    tau_p.y.x[] = 2*sq(y)*f[];
  //    u.x[] = 2.*Uavg*(1.-sq(y))*f[]*(1.-csibm[]);
    //   p[] = f[]*dp_l*(Ls-0.1*rc_nozzle-x);
    p[] = f[]*dp_l*(3.95-x);
    }
    
  // boundary ((scalar *){tau_p,tau_qq,u,p});


  }
  // boundary({f});

   create_ibm(csibm);
    
    foreach_dimension () {
      tau_p.x.x.refine = tau_p.x.y.refine = refine_bilinear;
      tau_p.x.x.prolongation = tau_p.x.y.prolongation = refine_bilinear;
      tau_p.x.x.gradient = tau_p.x.y.gradient = minmod;
    }
    tau_qq.prolongation = refine_bilinear;
    tau_qq.refine = refine_bilinear;
    tau_qq.gradient = minmod;

    u.x.gradient = minmod;
    u.y.gradient = minmod;

    unrefine(y > 8. || x>99.9); //3.999);

    // boundary ((scalar *){u,tau_p,tau_qq});

//}
}

/*event dp_increase (t=2; t+=1.; t<=15) {
  double dp_step = 1.;
  dp_l +=  dp_step;
//  p[left] = neumann(dp_l);
} */

event ibm (i++) {

  create_ibm(csibm);
  scalar fcsibm[];
  foreach() {
    fcsibm[] = csibm[];
  }
  fcsibm.refine = fcsibm.prolongation = refine_bilinear;
  // boundary({fcsibm});
  filter_vof(csibm,0,fcsibm);

  foreach() {
    u.x[] *= (1.-csibm[]);
    u.y[] *= (1.-csibm[]);
    
     foreach_dimension() {
     if (csibm[] > 0.1)
     {
       tau_p.x.x[] *= (1.-csibm[]) ;//min(tau_p.x.x[],100.);
       tau_p.y.x[] *= (1.-csibm[]) ;
     //  f[] = min(1.,f[]+csibm[]);
     }
     }
  }
  foreach_face() {
    double csf = (csibm[] + csibm[-1])/2.;
    uf.x[] *= (1.-csf); //fsibm.x[]);
  }

  // boundary((scalar *){uf,u,tau_p});

}

event projection (i++) {

  create_ibm(csibm);

   double U_IBM = 0.;

  face vector src_ibm[], uf_loc[];
  face vector av = a;

  foreach_face() {
      uf_loc.x[] = fm.x[]*(face_value(u.x,0) + dt*a.x[]);
      //src_ibm.x[] = 0.5*(scs[]+scs[-1])/dt*(U_IBM - uf_loc.x[]);
      //src_ibm.x[] = 0.5*(csibm[] + csibm[-1])/dt*(U_IBM - uf_loc.x[]);
      double csf = (csibm[] + csibm[-1])/2.; //fsibm.x[]; //(csibm[] + csibm[-1])/2.;
      src_ibm.x[] = csf/dt*(U_IBM - uf_loc.x[]);
      av.x[] += 1./fm.x[]*src_ibm.x[];
  }
  // boundary((scalar *){av,a});

//  trash({uf});
  foreach_face() {
      uf.x[] = fm.x[]*(face_value(u.x,0) + dt*a.x[]);
  }
  // boundary((scalar *){uf,av,a});
}


event properties (i++) {

  create_ibm(csibm);


 /* mask_dtau.refine = mask_dtau.prolongation = fraction_refine;
  mask_dtau_loc.refine = mask_dtau_loc.prolongation = fraction_refine;*/

//    foreach(){
//      if(x>25. & x<45.){
// 	    full_k[] = f[]/(y+0.001);	
//      }
//    }
  // scalar ftemp[];
  // foreach() {
  //   ftemp[] = f[];
  // }
  // ftemp.refine = ftemp.prolongation = refine_bilinear;
  // boundary({ftemp});
  // filter_vof(f,1,sf);
  foreach() {
    mupv[] = mup(f[]);
//    mupv[] *= (1.-csibm[]);
    lambdav[] = lambda(f[]);
//    lambdav[] *= (1.-csibm[]);
  }

  mupv.prolongation = refine_bilinear;
  lambdav.prolongation = refine_bilinear;
  // boundary((scalar *){mupv,lambdav});

}

/**
The grid will be refined with succesives *refine*. This strategy is not
very good (in fact, I think it is bad). However, since I only want to
test the axi part of th viscoelastic solver comparing with Bhat et
al. I have not worked on it. */

// event adapt (i+=100) {

// scalar s = tau_p.x.x;
// foreach() {
//     ffilt[] = pow((f[1,0]-f[-1,0])/(2.*Delta),2.) + pow((f[0,1]-f[0,-1])/(2.*Delta),2.);
//     fsfilt[] = pow((csibm[1,0]-csibm[-1,0])/(2.*Delta),2.) + pow((csibm[0,1]-csibm[0,-1])/(2.*Delta),2.);
//     Sxxfilt[] = pow((s[1,0]-s[-1,0])/(2.*Delta),2.) + pow((s[0,1]-s[0,-1])/(2.*Delta),2.);
//  }
//   ffilt.prolongation = refine_bilinear;
//   fsfilt.prolongation = refine_bilinear;
//   Sxxfilt.prolongation = refine_bilinear;

//   // boundary({ffilt,fsfilt,Sxxfilt});
//   //adapt_wavelet({fsfilt,u}, (double[]){0.001,0.01,0.01}, LEVEL);
//   // adapt_wavelet({ffilt,fsfilt,u}, (double[]){0.001,0.001,0.01,0.01}, LEVEL);
//   adapt_wavelet({ffilt,u,tau_p}, (double[]){0.001,0.01,0.01,0.001,0.001,0.01,0.01}, LEVEL);   
//   //refine(y>0.95 && y < 1.05 && level< LEVEL+1);
//  // refine (y<1.6 && level<LEVEL);
//   unrefine(y > 8. || x>99.9); //3.999);

// //   double xmin = HUGE; 
// //   //double xmax = HUGE;
// //   foreach (reduction (min:xmin))
// //     if (f[] > 0 && f[] < 1 && y<domain_size/pow(2.,LEVEL)) {
// //       coord p;
// //       coord n = mycs (point, f);
// //       double alpha = plane_alpha(f[],n);
// //       plane_area_center (n, alpha, &p);
// //       double posx  = (x + Delta*p.x); 
// //       if (posx < xmin)
// //         xmin = posx;
// //     }

// //   double yhmin = HUGE; 
// //   double xhmin = HUGE; 
// //   foreach (reduction (min:yhmin)){
// //     if (f[] > 0 && f[] < 1 && x<(xmin-1.)) {
// //       coord p;
// //       coord n = mycs (point, f);
// //       double alpha = plane_alpha(f[],n);
// //       plane_area_center (n, alpha, &p);
// //       double posy  = (y + Delta*p.y); 
// //       double posx = (x + Delta*p.x); 
// //       if (posy < yhmin){
// //         yhmin = posy;
// //         xhmin = posx;
// //      }
// //     }
// //   }

// // if(t>14.01){
// // //  refine(x<=(xhmin+20.) && x>=(xhmin-10) && x>5*Ls && y < 0.15  && level < 14);
// // refine(x>=25. && x<=45. && y < 0.15  && level < 14);
// // }

// }

event snapshot (t=0; t+=1e-2; t<=60) {
  char name0[80];
  sprintf(name0, "snapshot-%4.3f",t);
  scalar pid[];
  foreach() {
    pid[]  =fmod(pid()*(npe() + 37), npe());
    taupxx[] = tau_p.x.x[];
  }
  // boundary({pid});
  p.nodump = false;
  dump(name0);
}

event logfile (i++) {

  if (i == 0)
    fprintf (ferr,
             "t dt mgp.i mgpf.i mgu.i grid->tn perf.t perf.speed \n");
  fprintf (ferr, "%g, %g, %d, %d, %d, %ld, %g, %g \n",
           t, dt, mgp.i, mgpf.i, mgu.i,
           grid->tn, perf.t, perf.speed);
}

/*
event uprofile (t += 0.1) {
  char name1[80];
  sprintf (name1, "ux1_%g_%d.dat",t,pid());
  FILE * fp1 = fopen (name1, "w");
  int np = 15;
  double hh = 1.1/np;
  for (int j = 0; j <= np; j++)
    fprintf(fp1,"%g %g %g %g %g\n", hh*j, interpolate(u.x, 0.01*Ls, hh*j), interpolate(tau_p.x.x, 0.01*Ls, hh*j),interpolate(tau_p.y.x, 0.01*Ls, hh*j),interpolate(tau_p.y.y, 0.01*Ls, hh*j));
  fclose(fp1);

  char name2[80];
  sprintf (name2, "ux2_%g_%d.dat",t,pid());
  FILE * fp2 = fopen (name2, "w");
  for (int j = 0; j <= np; j++)
    fprintf(fp2,"%g %g %g %g %g\n", hh*j, interpolate(u.x, 0.25*Ls, hh*j), interpolate(tau_p.x.x, 0.25*Ls, hh*j),interpolate(tau_p.y.x, 0.25*Ls, hh*j),interpolate(tau_p.y.y, 0.25*Ls, hh*j));
  fclose(fp2);

  char name3[80];
  sprintf (name3, "ux3_%g_%d.dat",t,pid());
  FILE * fp3 = fopen (name3, "w");
  for (int j = 0; j <= np; j++)
    fprintf(fp3,"%g %g %g %g %g\n", hh*j, interpolate(u.x, 0.5*Ls, hh*j), interpolate(tau_p.x.x, 0.5*Ls, hh*j),interpolate(tau_p.y.x, 0.5*Ls, hh*j),interpolate(tau_p.y.y, 0.5*Ls, hh*j));
  fclose(fp3);

  char name4[80];
  sprintf (name4, "ux4_%g_%d.dat",t,pid());
  FILE * fp4 = fopen (name4, "w");
  for (int j = 0; j <= np; j++)
    fprintf(fp4,"%g %g %g %g %g\n", hh*j, interpolate(u.x, 0.75*Ls, hh*j), interpolate(tau_p.x.x, 0.75*Ls, hh*j),interpolate(tau_p.y.x, 0.75*Ls, hh*j),interpolate(tau_p.y.y, 0.75*Ls, hh*j));
  fclose(fp4);
 }*/


