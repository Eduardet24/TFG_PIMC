/*pimc.h*/

#ifndef _PIMC_H_
#define _PIMC_H_


// one component
void PIMCMoveOneByOne(int w);
void PIMCPseudopotentialMoveOneByOne(int w);
#ifdef PIMC_HARMONIC_NORMAL_MODES
void PIMCMoveHarmonicNormalModes(void);
#endif

struct Worm {
  DOUBLE x_ira;
  DOUBLE y_ira;
  DOUBLE z_ira;
  DOUBLE x_masha;
  DOUBLE y_masha;
  DOUBLE z_masha;
  int w;
  int i;
};

extern struct Worm worm;

#ifdef WORM_ALGORITHM
void InitWorm(void);
void FreeWorm(void);
void PIMCWormMoves(void);
void WormDiagnostic(void);
void SaveWormGapHistogram(const char *filename);
void MeasureWormPermutations(void);
void SaveWormPermutations(const char *filename);
void MeasureWormOBDM(void);
void SaveWormOBDM(const char *filename);
void MeasureWormSuperfluid(void);
DOUBLE ComputeWormSuperfluidFraction(long long samples, const DOUBLE A_sum[3], const DOUBLE A2_sum[3], const DOUBLE Icl_sum[3], DOUBLE W2_sum);
void SaveWormSuperfluid(const char *filename);
#endif

#endif
