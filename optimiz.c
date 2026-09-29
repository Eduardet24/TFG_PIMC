/*optimiz.c*/

#include <stdio.h>
#include "quantities.h"
#include "main.h"
#include "utils.h"
#include "trial.h"
#include "optimiz.h"
#include "display.h"
#include "randnorm.h"
#include "vmc.h"
#include "memory.h"
#include "compatab.h"
#include MATHINCLUDE 

DOUBLE opt_par[10][1]; // array with pointers to parameters to be optimized
DOUBLE opt_gradient[10];
int opt_N; // number of parameters
DOUBLE opt_U;

DOUBLE opt_Etarget = 2.70; // target velue of the energy
DOUBLE opt_sigma;

/******************************* Optimization *******************************/
void Optimization(void) {
  int i=0;

  Message("\n\nStarting automatic variational optimization of parameters\n");

  LoadCoordinatesOptimization();

  E = Energy(&Epot, &Ekin, &EFF, &Edamping, &Eint, &Eext);
  OptimizationEnergy(opt_Etarget, &E, &opt_sigma);
  Message("  Starting energy per particle is %.15e\n", E/(double)N);

  for(i=0; i<Niter*blck; i++) {
    Message("\n");
#ifdef CRYSTAL
#ifdef CRYSTAL_WIDTH_ARRAY
    //Message("E= %g parameters: %.3g %.3g %.3g %.3g %.3g  %.3g %.3g %.3g %.3g %.3g", E, Crystal.Rz[0], Crystal.Rz[1], Crystal.Rz[2], Crystal.Rz[3], Crystal.Rz[4], Crystal.z[0], Crystal.z[1], Crystal.z[2], Crystal.z[3], Crystal.z[4]);
    //Message("E= %g sigma= %g parameters: (%.3g %.3g) (%.3g %.3g) (%.3g %.3g) (%.3g %.3g) (%.3g %.3g)", E, opt_sigma, Crystal.z[0], Crystal.Rz[0], Crystal.z[1], Crystal.Rz[1], Crystal.z[2], Crystal.Rz[2], Crystal.z[3], Crystal.Rz[3], Crystal.z[4], Crystal.Rz[4]);
    //Message("E= %g sigma= %g parameters: (%.3g %.3g)", E, opt_sigma, Crystal.z[0], Crystal.Rz[0]);
    // Crystal.z[i]
    //opt_gradient[0] = OptimizationGradient(&Crystal.Rz[0]);

    Message("E= %.8g sigma= %.5g ", E/(double)N, opt_sigma);

    //ucr += Crystal_dot_weight_j*Exp(-Crystal_dot_Rx_j*dr[0]*dr[0]-Crystal_dot_Ry_j*dr[1]*dr[1]-Crystal_dot_Rz_j*dr[2]*dr[2]);

//0 - 0. 0. 0. 0. 0. $alpha 1.
//1 - 0. 0. 0. 0. 0. $sho $aho
//2 - 0. 0. 0. 0. 0. $s0 $a0
//3 - 0. 0. $r1 0. 0. $s1 $a1 
//4 - 0. 0. $mr1 0. 0. $s1 $a1 
//5 - 0. 0. $r2 0. 0. $s2 $a2
//6 - 0. 0. $mr2 0. 0. $s2 $a2
//7 - 0. 0. $r3 0. 0. $s3 $a3 
//8 - 0. 0. $mr3 0. 0. $s3 $a3\n";

    Message("p0=(%.3g,", Crystal.Rz[0]);
    OptimizationAnnealing(&Crystal.Rz[0]);
    Message("%.3g) ", Crystal.weight[0]);
    OptimizationAnnealing(&Crystal.weight[0]);

    Message("p1=(%.3g,", Crystal.Rz[3]);
    OptimizationAnnealing(&Crystal.Rz[3]);
    Message("%.3g) ", Crystal.weight[3]);
    OptimizationAnnealing(&Crystal.weight[3]);

    Message("m1=(%.3g,", Crystal.Rz[4]);
    OptimizationAnnealing(&Crystal.Rz[4]);
    Message("%.3g) ", Crystal.weight[4]);
    OptimizationAnnealing(&Crystal.weight[4]);

    Crystal.Rz[3] = Crystal.Rz[4] = 0.5*(Crystal.Rz[3] + Crystal.Rz[4]);
    Crystal.weight[3] = Crystal.weight[4] = 0.5*(Crystal.weight[3] + Crystal.weight[4]);

    Message("p2=(%.3g,", Crystal.Rz[5]);
    OptimizationAnnealing(&Crystal.Rz[5]);
    Message("%.3g) ", Crystal.weight[5]);
    OptimizationAnnealing(&Crystal.weight[5]);

    Message("m2=(%.3g,", Crystal.Rz[6]);
    OptimizationAnnealing(&Crystal.Rz[6]);
    Message("%.3g) ", Crystal.weight[6]);
    OptimizationAnnealing(&Crystal.weight[6]);

    Crystal.Rz[5] = Crystal.Rz[6] = 0.5*(Crystal.Rz[5] + Crystal.Rz[6]);
    Crystal.weight[5] = Crystal.weight[6] = 0.5*(Crystal.weight[5] + Crystal.weight[6]);

    //OptimizationAnnealing(&Crystal.Rz[0]);
    //OptimizationAnnealing(&Crystal.Rz[1]);
    //OptimizationAnnealing(&Crystal.Rz[2]);
    //OptimizationAnnealing(&Crystal.Rz[3]);
    //OptimizationAnnealing(&Crystal.Rz[4]);

    //OptimizationAnnealing(&Crystal.z[0]);
    //OptimizationAnnealing(&Crystal.z[1]);
    //OptimizationAnnealing(&Crystal.z[2]);
    //OptimizationAnnealing(&Crystal.z[3]);
    //OptimizationAnnealing(&Crystal.z[4]);
    //OptimizationAnnealing(&orbital_weight2);
    //OptimizationAnnealing(&orbital_weight3);
#endif
#endif
    if(video) Picture();
  }
}

/************************ Optimization Gradient *****************************/
DOUBLE OptimizationGradient(double *par) {
  DOUBLE Ep;
  DOUBLE par_store;
  DOUBLE gradient;
  DOUBLE dx = 1e-3;
  DOUBLE sigma, sigma_p;

  par_store = *par;

  OptimizationEnergy(opt_Etarget, &E, &sigma);

  *par *= 1. + dx;

  OptimizationEnergy(opt_Etarget, &Ep, &sigma_p);

  //gradient = (Ep-E)/(*par-par_store); // energy gradient
  gradient = (sigma_p-sigma)/sigma; // sigma gradient

  //Message("E= %e gradient %.3e, par %e -> ", E, gradient, *par);
  *par = par_store - gradient*1e-1;
/*  if(gradient>0)
    *par = par_store/1.1;
  else
    *par = par_store*1.1;*/
  //Message("%e ", *par);

  E = Ep;
  opt_sigma = sigma_p;

  return gradient;
}

/************************ Optimization Gradient *****************************/
void OptimizationAnnealing(double *par) {
  DOUBLE Ep, dE;
  DOUBLE sigmap, dsigma;
  DOUBLE par_store;
  DOUBLE xi;
  DOUBLE dx,dy;
  static DOUBLE T = 1e-15; // effective temperature
  static DOUBLE dt = 0.01;
  DOUBLE dOpt;

  par_store = *par;

  OptimizationEnergy(opt_Etarget, &E, &opt_sigma);

  RandomNormal(&dx, &dy, 0., sqrt(dt));
  *par *= fabs(1+dx);

  OptimizationEnergy(opt_Etarget, &Ep, &sigmap);

  dE = Ep - E;
  dsigma = sigmap - opt_sigma;

  //dOpt = dE; // use energy for optimization 
  dOpt = dsigma; // use dispersion for optimization 

  if(dOpt <= 0.) {
    accepted++;
    E = Ep;
    opt_sigma = sigmap;
  }
  else {
    xi = Random();
    if(xi<exp(-dOpt/T)) {
      accepted++;
      E = Ep;
      opt_sigma = sigmap;
    }
    else {
      rejected++;
      *par = par_store;
    }
  }

  T *= 0.99;
}


/********************************* Energy ************************************/
// returns energy of the system averaged over all walkers
// and variance, relative to Etarget
void OptimizationEnergy(double Etarget, double *Emean, double *sigma) {
  DOUBLE dEpot, dEkin, dEff, Edamping;
  DOUBLE weight, Weight = 0;
  DOUBLE Ep;
  int w;
  //DOUBLE Epot, Ekin, EFF;
  //DOUBLE Etot = 0.;

  *Emean = 0.; // energy
  *sigma = 0.; // variance

  //Epot = Ekin = EFF = 0.;
  for(w=0; w<Nwalkers; w++) {
    //E = W[w].Epot + W[w].Ekin; // old energy
    Ep = WalkerEnergy(&W[w], &dEpot, &dEkin, &dEff, &Edamping, &Eint, &Eext); // new energy

    //weight = 1.;
    weight = exp(2.*(U(W[w])-W[w].U));
  
    //Etot += E; // mean "old" energy
    Weight += weight;
    *Emean += Ep * weight; // mean "new" energy
    *sigma += (Ep-Etarget)*(Ep-Etarget);

    //Epot += dEpot * weight;
    //Ekin += dEkin * weight;
    //EFF  += dEff * weight;
  }

  //Etot /= (DOUBLE) (N * Nwalkers);
  *Emean /= (DOUBLE) (N * Nwalkers); // new energy integral
  Weight /= (DOUBLE) (N * Nwalkers); // new norm
  *Emean /= Weight;

  //Epot /= (DOUBLE) (N * Nwalkers); // new energy integral
  //Ekin /= (DOUBLE) (N * Nwalkers); // new energy integral
  //EFF /= (DOUBLE) (N * Nwalkers); // new energy integral

  //Epot /= Weight;
  //Ekin /= Weight;
  //EFF /= Weight;

  
  // -= Etot; // difference
  *sigma /= (DOUBLE) (N * Nwalkers);
  *sigma = sqrt(*sigma);
}

/****************************** Save Coordinates ************************/
int SaveCoordinatesOptimization(void) {
  FILE *out;
  int w, i;
  static int walkers_saved = 0;

  out = fopen("3Dprev.dat", walkers_saved?"a":"w");
  if(out == NULL) {
    perror("\nError:");
    Warning("can't write to 3Dprev.dat file\n");
    return 1;
  }

  //fprintf(out, "%i\n", Nwalkers);
  walkers_saved += Nwalkers;
  for(w=0; w<Nwalkers; w++) {
#ifdef TRIAL_1D
    for(i=0; i<N; i++) fprintf(out, "%.7f ", W[w].z[i]);
    fprintf(out, "\n");
#else
    for(i=0; i<N; i++) fprintf(out, "%.7f %.7f %.7f\n", W[w].x[i], W[w].y[i], W[w].z[i]);
#endif
  }
  fclose(out);

  out = fopen("3DprevNw.dat", "w");
  if(out == NULL) {
    perror("\nError:");
    Warning("can't write to 3DprevNw.dat file\n");
    return 1;
  }

  fprintf(out, "%i\n", walkers_saved);
  fclose(out);

  return 0;
}

/****************************** Load Coordinates ************************/
int LoadCoordinatesOptimization(void) {
  FILE *in;
  int w, i;

  Message("  Loading stored coordinates ...\n");

  in = fopen("3DprevNw.dat", "r");
  if(in == NULL) {
    perror("\nError:");
    Error("can't open 3DprevNw.dat file\n");
  }

  fscanf(in, "%i\n", &Nwalkers);
  fclose(in);

  NwalkersMax = Nwalkers;

  Message("  Allocating memory for %i walkers ... ", Nwalkers);
  AllocateWalkers();
  Message("done\n  Loading walkers ...");

  in = fopen("3Dprev.dat", "r");
  if(in == NULL) {
    perror("\nError:");
    Warning("can't open 3Dprev.dat file\n");
    return 1;
  }

  for(w=0; w<Nwalkers; w++) {
    //fscanf(in, "%" LF, &W[w].E);
    for(i=0; i<N; i++) 
#ifdef TRIAL_1D
    fscanf(in, " %" LF, &W[w].z[i]);
#else
    fscanf(in, "%"LF " %" LF " %" LF "\n", &W[w].x[i], &W[w].y[i], &W[w].z[i]);
#endif
    W[w].U = U(W[w]); // initialize the weight
  }
  fclose(in);

  Message("  done\n");

  return 0;
}
