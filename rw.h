/*rw.h*/

#ifndef _RW_H_
#define _RW_H_

#include "main.h"

int LoadParameters(void);
int LoadParticleCoordinates(void);
void LoadParametersSk(void);

void SaveWaveFunction(struct Grid *G, char *file_name, DOUBLE min, DOUBLE max);
void SaveImpurityWaveFunction(void);
int SaveCoordinates(char *file_particles);
int SaveEnergyVMC(DOUBLE E, DOUBLE EFF);
int SaveEnergyDMC(DOUBLE E);
int SaveEnergyCls(DOUBLE E);
int SaveOBDM(void);
int SaveOBDMpure(void);
int SaveTBDM(void);
int SaveOBDMMatrix(void);
int SaveTBDMMatrix(void);
int SavePairDistribution(void);
int SaveRadialDistribution(void);
int SaveRadialZDistribution(void);
int LoadRadialZDistribution(void);
int SaveR2(DOUBLE f1, DOUBLE f2);
int SaveZ2(DOUBLE f1, DOUBLE f2);
int SavePureR2(DOUBLE f);
int SavePureZ2(DOUBLE f);
int SaveSuperfluidDensity(void);
/*int SavePureDistribution(void);
int SaveCoordinates(void);
int SaveControlVariables(void);*/
void SaveBlockData(int);
void SaveMeanR2(void);
void SaveMeanR2DMC(void);
void SaveOrderParameter(void);
int LoadSuperfluidDensityArray(void);
int SaveSuperfluidDensity5Array(void);
int SaveStaticStructureFactor(void);
int SavePairDistribtuionMatrix(void);
int SaveMomentumDistributionPure(void);
int SaveRadialDistributionMatrix(void);
int LoadCrystalCoordinates(void);
int SaveFormFactor(void);
int SaveDefectBarrier(void);
int SaveMomentumDistributionMatrix(void);
void SaveLattice(void);
void SaveVext(void);
void SaveOneBodyWaveFunction(void);
void SaveThreeBodyWaveFunction(void);
void SaveSpinPolarization(void);
int SaveHyperRadius(void);
int SaveEffectivePotential(void);
void CopyPureCoordinates(void);
int SavePureCoordinates(void);
int SaveEpotPure(void);
int SavePairDistribtuionR2(void);

#define Fopen(out, fname, mode, text) \
  out = fopen(fname, mode);\
  if(out == NULL) {        \
    perror("\nError:");    \
    Warning("Can't load/save %s from/to file %s\n", text, fname);\
  }

#endif
