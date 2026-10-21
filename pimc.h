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

#endif
