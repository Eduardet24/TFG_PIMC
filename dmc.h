/*DMC.h*/

#ifndef _DMC_H_
#define _DMC_H_

#define DMC_LINEAR_ONE_SYMMETRIC_MOVE

void DMCQuadraticMove(int w);
void DMCLinearMoveDrift(int w);
void DMCLinearMoveGaussian(int w);
void DMCLinearMetropolisMove(int w);
void DMCQuadraticMetropolisMove(int w);
void DMCQuarticMove(int w);
void DMCPseudopotentialMove(int w);
void DMCLinearMoveOne(int w);

void GaussianJump(int w, DOUBLE dt);
void GaussianJumpWalker(int w, DOUBLE dt);

void DriftJump(DOUBLE dt, int w);
void DriftJump1(DOUBLE dt, int w);
void DriftJump2(DOUBLE dt, int w);

void DriftForce(DOUBLE **F, DOUBLE **R, int w);

void Branching(void);
void BranchingPseudopotential(void);

void CopyVectorToWalker(struct Walker* W, DOUBLE** R, int *spin);
void CopyWalkerToVector(DOUBLE** R, struct Walker W);
void CopyVector(DOUBLE** out, DOUBLE** in);
DOUBLE EnergyO(void);
DOUBLE EnergyODMC(void);

void CalculateDriftForceAndKineticEnergyNumerically(struct Walker *walker, DOUBLE *Ekin);
void WalkersSort(void); // put all walkers into order: ALIVE, update dead_walkers_storage
#endif
