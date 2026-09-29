/*tvmc.c*/

#include <stdio.h>
#include <math.h>
#include <complex.h>
//#undef I
#include "main.h"
#include "compatab.h"
#include "utils.h"
#include "memory.h"
#include "rw.h"

#ifdef TRIAL_TONKS_TRAP_tVMC // f(x) = 1 - x/a  + (c3R + i c3I) x^2 + (c4R + i c4I) x^3 + (c5R + i c5I) x^4
DOUBLE InterpolateExactRef(DOUBLE x) {
// real part of the Jastrow term
  DOUBLE f2body, sumR;

  f2body = 1.-x/a;
  sumR = cR[3]*x*x + cR[4]*x*x*x + cR[5]*x*x*x*x;
  return f2body + sumR;
}

DOUBLE InterpolateExactImf(DOUBLE x) {
// imaginary part of the Jastrow term
  DOUBLE sumI;

  sumI = cI[3]*x*x + cI[4]*x*x*x + cI[5]*x*x*x*x;

  return sumI;
}

DOUBLE dfdc(DOUBLE x, int i) {
// real nad imaginary part of d f / dc(i)
// sumR = cR[3]*x*x + cR[4]*x*x*x + cR[5]*x*x*x*x;
  if(i == 2)
    return x*x;
  else if(i == 3)
    return x*x*x;
  else if(i == 4)
    return x*x*x*x;
  else
    return 0;
}

DOUBLE InterpolateExactU(DOUBLE x) {
// U = (1/2) ln [f*f]
// f*f = InterpolateExactRef(x)**2 +InterpolateExactImf(x)**2
  DOUBLE f2body, sumR, sumI, fabs2;

  f2body = 1.-x/a;
  sumR = cR[3]*x*x + cR[4]*x*x*x + cR[5]*x*x*x*x;
  sumI = cI[3]*x*x + cI[4]*x*x*x + cI[5]*x*x*x*x;

  fabs2 = f2body*f2body + 2.*f2body*sumR + sumR*sumR + sumI*sumI;

  return 0.5*Log(fabs2);
}

DOUBLE InterpolateExactRefp(DOUBLE x) {
// real part of f'(x)
  DOUBLE f2body, sumR;

  f2body = -1/a;
  sumR = 2.*cR[3]*x + 3.*cR[4]*x*x + 4.*cR[5]*x*x*x;
  return f2body + sumR;
}

DOUBLE InterpolateExactImfp(DOUBLE x) {
// imaginary part of f'(x)
  DOUBLE sumI;

  sumI = 2.*cI[3]*x + 3.*cI[4]*x*x + 4.*cI[5]*x*x*x;

  return sumI;
}

DOUBLE InterpolateExactRefpp(DOUBLE x) {
// real part of f"(x)
  DOUBLE sumR;

  sumR = 2.*cR[3] + 6.*cR[4]*x + 12.*cR[5]*x*x;
  return sumR;
}

DOUBLE InterpolateExactImfpp(DOUBLE x) {
// imaginary part of f"(x)
  DOUBLE sumI;

  sumI = 2.*cI[3] + 6.*cI[4]*x + 12.*cI[5]*x*x;

  return sumI;
}

// Fp = f' / f
DOUBLE InterpolateExactFp(DOUBLE x) {
  return InterpolateExactRefp(x) / InterpolateExactRef(x);
}

DOUBLE InterpolateExactE(DOUBLE x) {
//1D Eloc = [-f"/f] + (f'/f)^2
  double Fp = InterpolateExactFp(x);

  return -InterpolateExactRefpp(x) / InterpolateExactRef(x) + Fp*Fp;
}
#endif

/*********************** t-VMC measurement *******************************/
#define NPARTICLES_MAX 1000
#define NPARAMETERS 10
#define NOBSERVABLES 10

double ReMatrix[NOBSERVABLES][NPARAMETERS] = {0};
double ImMatrix[NOBSERVABLES][NPARAMETERS] = {0};
double ReVector[NOBSERVABLES] = {0};
double ImVector[NOBSERVABLES] = {0};
//int Matrix_times_measured = 0;
// observable k: sum_i |x_i - xcm|^k

// wf. exp( - alpha x^2 - i beta x^2)
int MeasuretVMC(int iter) {
#ifdef TRIAL_TONKS_TRAP_tVMC
  FILE *out;
  static int first_time = ON;
  static int times_measured = 0;
  int w, i, j, k, l;
  double xij, x2walker;
  //double unity_av, x2av, x2sum2av, ReH, ReHx2, Rex2H, ImH, ImHx2, Imx2H;
  double ReD1[NPARTICLES_MAX];
  double ImD1[NPARTICLES_MAX];
  double ReD2[NPARTICLES_MAX];
  double ImD2[NPARTICLES_MAX];
  double ReD2T[NPARTICLES_MAX];
  double ImD2T[NPARTICLES_MAX];
  double O[NOBSERVABLES]; // observables
  double zCM;// , x2;

  if(tvmcNpar>NPARAMETERS) Error(" tVMC : increase number of parameters NPARAMETERS and recompile");
  if(tvmcNobs>NOBSERVABLES) Error(" tVMC : increase number of parameters NOBSERVABLES and recompile");

//#pragma omp parallel for
  for(w=0; w<Nwalkers; w++) {
    // empty arrays
    ArrayEmpty1D(ReD1, i, N);
    ArrayEmpty1D(ImD1, i, N);
    ArrayEmpty1D(ReD2, i, N);
    ArrayEmpty1D(ImD2, i, N);
    ArrayEmpty1D(ReD2T, i, N);
    ArrayEmpty1D(ImD2T, i, N);
    //Matrix_times_measured++;
    times_measured++;

    // calculate observables
    // O_k = <|r_i - r_{CM}|^k>
    zCM = 0.; // calculate center of mass
    for(i=0; i<N; i++) zCM += W[w].z[i];
    zCM /= (DOUBLE) N;
    xij = fabs(W[w].z[0] - W[w].z[1]);
    for(k=0; k<tvmcNobs; k++) {
      O[k] = 0.;
      //for(i=0; i<N; i++) O[k] += pow(fabs(W[w].z[i] - zCM), k);
      O[k] += pow(xij, k); // two particles
    }

    // average x^2
    x2walker = 0.;
    for(i=0; i<N; i++) x2walker += W[w].z[i]*W[w].z[i];

    for(i=0; i<N; i++) {
      ReD1[i] = -2.*cR[2]*W[w].z[i]; // Drift Force (nabla Psi) / Psi
      ImD1[i] = -2.*cI[2]*W[w].z[i];
      ReD2[i] = -2.*cR[2]; // gradient of the Drift Force [(nabla Psi) / Psi]
      ImD2[i] = -2.*cI[2];
      for(j=0; j<N; j++) {
        if(j != i) {
          DOUBLE RefJ, ImfJ, RefpJ, ImfpJ, RefppJ, ImfppJ, ReFpJ, ImFpJ, s;
          xij = fabs(W[w].z[i] - W[w].z[j]);
          RefJ = InterpolateExactRef(xij);
          ImfJ = InterpolateExactImf(xij);
          RefpJ = InterpolateExactRefp(xij);
          ImfpJ = InterpolateExactImfp(xij);
          RefppJ = InterpolateExactRefpp(xij);
          ImfppJ = InterpolateExactImfpp(xij);
          s = ((W[w].z[i] > W[w].z[j]) ? (1.) : (-1.)); // sign
          // f(r) = xij - a -> ReD1[i] += ((W[w].z[i] > W[w].z[j])?(1.):(-1.)) / (xij - a);
          // (a + i b) * (c - id) = ac + bd + i(bc - ad)
          // (RefpJ + i ImfpJ) / (RefJ + i ImfJ) :
          ReFpJ = s*(RefpJ*RefJ + ImfpJ*ImfJ) / (RefJ*RefJ + ImfJ*ImfJ);
          ImFpJ = s*(ImfpJ*RefJ - RefpJ*ImfJ) / (RefJ*RefJ + ImfJ*ImfJ);
          ReD1[i] += ReFpJ;
          ImD1[i] += ImFpJ;

          // f(r) = xij - a ->  ReD2[i] -= 1./ ((xij - a)*(xij - a));
          // (ReFp + i ImFp)(ReFp + i ImFp) = ReFp*ReFp - ImFp*ImFp + i 2 ReFp ImRp
          // (RefppJ + i ImfppJ) / (RefJ + i ImfJ) :
          ReD2[i] += (RefppJ*RefJ + ImfppJ*ImfJ) / (RefJ*RefJ + ImfJ*ImfJ) - (ReFpJ*ReFpJ - ImFpJ*ImFpJ);
          ImD2[i] += (ImfppJ*RefJ - RefppJ*ImfJ) / (RefJ*RefJ + ImfJ*ImfJ) - 2.*ReFpJ*ImFpJ;

          for(k=0; k<tvmcNobs; k++) { // other parameters
            for(l=2; l<tvmcNpar; l++) {
              DOUBLE ReFc;
              ReFc = dfdc(xij, l);

              // dfdc / fJ
              // a / (b + ic) = a(b - ic)/(b^2+c^2) 
              // 0.5 is introduced to avoid double counting
              ReMatrix[k][l] += 0.5*O[k] * ReFc*RefJ / (RefJ*RefJ + ImfJ*ImfJ);
              ImMatrix[k][l] -= 0.5*O[k] * ReFc*ImfJ / (RefJ*RefJ + ImfJ*ImfJ);
            }
          }
        }
      }
    }

    for(i=0; i<N; i++) {
      double ReH, ImH;
      ReD2T[i] = ReD2[i] + ReD1[i]*ReD1[i] - ImD1[i]*ImD1[i];
      ImD2T[i] = ImD2[i] + 2.*ReD1[i]*ImD1[i];

      ReH = -0.5*ReD2T[i] + 0.5*omega_z * omega_z * W[w].z[i] * W[w].z[i];
      ImH = -0.5*ImD2T[i];

      for(k=0; k<tvmcNobs; k++) {
        ReVector[k] += ReH * O[k];
        ImVector[k] += ImH * O[k];
      }

      //ReEkin += -0.5*ReD2T[i]; // needed for test
      //ImEkin += -0.5*ImD2T[i];

      //ReHx2 += -1. - 2.*W[w].z[i]*ReD1[i] + x2walker*(-0.5*ReD2T[i] + 0.5*W[w].z[i]*W[w].z[i]);
      //ImHx2 += -2.*W[w].z[i]*ImD1[i] + x2walker*(-0.5*ImD2T[i]);

      //Rex2H += x2walker*(-0.5*ReD2T[i] + 0.5*W[w].z[i]*W[w].z[i]);
      //Imx2H += x2walker*(-0.5*ImD2T[i]);
    }
    for(k=0; k<tvmcNobs; k++) {
      ReMatrix[k][0] += O[k]; // zeroth parameter is the normalization
      ReMatrix[k][1] += O[k]*x2walker; // first parameter exp( c x2walker)
    }
  }
  //unity_av /= (double) (Nwalkers);
  //x2av /= (double) (Nwalkers);
  //x2sum2av /= (double) (Nwalkers);
  //ReH /= (double) (Nwalkers);
  //ImH /= (double) (Nwalkers);
  //ReHx2 /= (double) (Nwalkers);
  //ImHx2 /= (double) (Nwalkers);
  //Rex2H /= (double) (Nwalkers);
  //Imx2H /= (double) (Nwalkers);

  //if(first_time == ON) fprintf(out, "#1 <x2> 2 <x2>^2 3-4 <H> 5-6 <Hx^2> 7-8 <x^2H>\n");
  //fprintf(out, "%.15" LE " %.15" LE " %.15" LE " %.15" LE " %.15" LE  " %.15" LE  " %.15" LE " %.15" LE "\n", x2av, x2sum2av, ReH, ImH, ReHx2, ImHx2, Rex2H, Imx2H);
  //fclose(out);

  if(iter % Niter == 0) {
    Fopen(out, "tvmcmatrixRe.dat", "w", "tvmcmatrixRe.dat");
    for(k=0; k<tvmcNobs; k++) {
      for(l=0; l<tvmcNpar; l++) {
        //fprintf(out, "%.15" LE "+I*%.15" LE " ", ReMatrix[k][l] / (DOUBLE) Matrix_times_measured, ImMatrix[k][l] / (DOUBLE) Matrix_times_measured);
        fprintf(out, "%.15" LE " ", ReMatrix[k][l] / (DOUBLE) times_measured);
      }
      fprintf(out, "\n");
    }
    fclose(out);

    Fopen(out, "tvmcmatrixIm.dat", "w", "tvmcmatrixIm.dat");
    for(k=0; k<tvmcNobs; k++) {
      for(l=0; l<tvmcNpar; l++) {
        //fprintf(out, "%.15" LE "+I*%.15" LE " ", ReMatrix[k][l] / (DOUBLE) Matrix_times_measured, ImMatrix[k][l] / (DOUBLE) Matrix_times_measured);
        fprintf(out, "%.15" LE " ", ImMatrix[k][l] / (DOUBLE) times_measured);
      }
      fprintf(out, "\n");
    }
    fclose(out);

    Fopen(out, "tvmcvectorRe.dat", "w", "tvmcvectorRe.dat");
    for(k=0; k<tvmcNobs; k++) {
      //fprintf(out, "%.15" LE "+I*%.15" LE " ", ReMatrix[k][l] / (DOUBLE) Matrix_times_measured, ImMatrix[k][l] / (DOUBLE) Matrix_times_measured);
      fprintf(out, "%.15" LE "\n", ReVector[k] / (DOUBLE) times_measured);
    }
    fclose(out);

    Fopen(out, "tvmcvectorIm.dat", "w", "tvmcvectorIm.dat");
    for(k=0; k<tvmcNobs; k++) {
        //fprintf(out, "%.15" LE "+I*%.15" LE " ", ReMatrix[k][l] / (DOUBLE) Matrix_times_measured, ImMatrix[k][l] / (DOUBLE) Matrix_times_measured);
      fprintf(out, "%.15" LE "\n", ImVector[k] / (DOUBLE) times_measured);
    }
    fclose(out);
  }

  first_time = OFF;

#endif
  return 0;
}