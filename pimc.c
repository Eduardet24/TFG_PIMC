/** @file pimc.c
 *  @brief Implementation of the PIMC algorithm.
*/
/*pimc.c*/
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <memory.h>
#include "main.h"
#include "randnorm.h"
#include "rw.h"
#include "compatab.h"
#include "memory.h"
#include "vmc.h"
#include "mymath.h"
#include "pimc.h"
#include "pigs.h"
#include "utils.h"
#include MATHINCLUDE

#ifdef WORM_ALGORITHM
static int *worm_ring_visited = NULL;
static int *worm_ring_p = NULL;
static int *worm_ring_w = NULL;
static DOUBLE *worm_ring_x = NULL;
static DOUBLE *worm_ring_y = NULL;
static DOUBLE *worm_ring_z = NULL;
static int worm_ring_cap_n = 0;
static int worm_ring_cap_total = 0;

static int EnsureWormRingScratch(int n_particles, int total_beads) {
  if(n_particles > worm_ring_cap_n) {
    int *new_vis = (int*) realloc(worm_ring_visited, (size_t)n_particles * sizeof(int));
    if(new_vis == NULL) return 0;
    worm_ring_visited = new_vis;
    worm_ring_cap_n = n_particles;
  }
  if(total_beads > worm_ring_cap_total) {
    int *new_p = (int*) realloc(worm_ring_p, (size_t)total_beads * sizeof(int));
    int *new_w = (int*) realloc(worm_ring_w, (size_t)total_beads * sizeof(int));
    DOUBLE *new_x = (DOUBLE*) realloc(worm_ring_x, (size_t)total_beads * sizeof(DOUBLE));
    DOUBLE *new_y = (DOUBLE*) realloc(worm_ring_y, (size_t)total_beads * sizeof(DOUBLE));
    DOUBLE *new_z = (DOUBLE*) realloc(worm_ring_z, (size_t)total_beads * sizeof(DOUBLE));
    if(new_p != NULL) worm_ring_p = new_p;
    if(new_w != NULL) worm_ring_w = new_w;
    if(new_x != NULL) worm_ring_x = new_x;
    if(new_y != NULL) worm_ring_y = new_y;
    if(new_z != NULL) worm_ring_z = new_z;
    if(new_p == NULL || new_w == NULL || new_x == NULL || new_y == NULL || new_z == NULL) return 0;
    worm_ring_cap_total = total_beads;
  }
  return 1;
}
#endif

#ifdef PIMC_HARMONIC_NORMAL_MODES
/******************* Exact Harmonic PIMC Normal-Mode Draw ******************/
#define NM_CACHE_SLOTS 64

typedef struct {
  int ring_beads;
  DOUBLE tau;
  DOUBLE omega2;
  DOUBLE inv_sqrt_n;
  DOUBLE amp_mode;
  DOUBLE sigma0;
  DOUBLE sigma_nyq;
  int n_modes;
  int cap_modes;
  DOUBLE *sigma_mode;
  DOUBLE *step_cos;
  DOUBLE *step_sin;
} NormalModeCache1D;

static NormalModeCache1D nm_cache[NM_CACHE_SLOTS];
static int nm_cache_next = 0;
static int nm_cache_atexit_registered = 0;
static DOUBLE *nm_plain_ring = NULL;
static int nm_plain_ring_cap = 0;

static void FreeNormalModeCache(void) {
  int i;
  for(i = 0; i < NM_CACHE_SLOTS; i++) {
    if(nm_cache[i].sigma_mode != NULL) free(nm_cache[i].sigma_mode);
    if(nm_cache[i].step_cos != NULL) free(nm_cache[i].step_cos);
    if(nm_cache[i].step_sin != NULL) free(nm_cache[i].step_sin);
    nm_cache[i].sigma_mode = NULL;
    nm_cache[i].step_cos = NULL;
    nm_cache[i].step_sin = NULL;
    nm_cache[i].ring_beads = 0;
    nm_cache[i].cap_modes = 0;
  }
  if(nm_plain_ring != NULL) {
    free(nm_plain_ring);
    nm_plain_ring = NULL;
  }
  nm_plain_ring_cap = 0;
}

static const NormalModeCache1D* GetNormalModeCache(int ring_beads, DOUBLE tau, DOUBLE omega2) {
  int i, mode;
  int n_modes = (ring_beads - 1) / 2;
  NormalModeCache1D *slot = NULL;

  if(!nm_cache_atexit_registered) {
    atexit(FreeNormalModeCache);
    nm_cache_atexit_registered = 1;
  }

  for(i = 0; i < NM_CACHE_SLOTS; i++) {
    if(nm_cache[i].ring_beads == ring_beads &&
       nm_cache[i].tau == tau &&
       nm_cache[i].omega2 == omega2) {
      return &nm_cache[i];
    }
    if(nm_cache[i].ring_beads == 0 && slot == NULL) {
      slot = &nm_cache[i];
    }
  }

  if(slot == NULL) {
    slot = &nm_cache[nm_cache_next];
    nm_cache_next = (nm_cache_next + 1) % NM_CACHE_SLOTS;
  }

  if(n_modes > slot->cap_modes) {
    DOUBLE *new_sig = (DOUBLE*) realloc(slot->sigma_mode, (size_t)n_modes * sizeof(DOUBLE));
    DOUBLE *new_cos = (DOUBLE*) realloc(slot->step_cos, (size_t)n_modes * sizeof(DOUBLE));
    DOUBLE *new_sin = (DOUBLE*) realloc(slot->step_sin, (size_t)n_modes * sizeof(DOUBLE));
    if(new_sig != NULL) slot->sigma_mode = new_sig;
    if(new_cos != NULL) slot->step_cos = new_cos;
    if(new_sin != NULL) slot->step_sin = new_sin;
    if(new_sig == NULL || new_cos == NULL || new_sin == NULL) {
      Error("Out of memory in GetNormalModeCache\n");
    }
    slot->cap_modes = n_modes;
  }

  slot->ring_beads = ring_beads;
  slot->tau = tau;
  slot->omega2 = omega2;
  slot->n_modes = n_modes;
  slot->inv_sqrt_n = 1.0 / Sqrt((DOUBLE)ring_beads);
  slot->amp_mode = Sqrt(2.0 / (DOUBLE)ring_beads);
  slot->sigma0 = 1.0 / Sqrt(tau * omega2);
  slot->sigma_nyq = (ring_beads % 2 == 0) ? (1.0 / Sqrt(4.0 / tau + tau * omega2)) : 0.0;

  for(mode = 1; mode <= n_modes; mode++) {
    DOUBLE step = 2.0 * PI * (DOUBLE)mode / (DOUBLE)ring_beads;
    DOUBLE half_sin = Sin(0.5 * step);
    DOUBLE eigenvalue = 4.0 * half_sin * half_sin / tau + tau * omega2;
    slot->step_cos[mode - 1] = Cos(step);
    slot->step_sin[mode - 1] = Sin(step);
    slot->sigma_mode[mode - 1] = 1.0 / Sqrt(eigenvalue);
  }

  return slot;
}

static void InitHarmonicRingDCAndNyquist(DOUBLE *ring_out, int ring_beads, const NormalModeCache1D *c) {
  int n;
  DOUBLE z0, z1, dc;

  RandomNormal(&z0, &z1, 0., 1.0);
  dc = z0 * c->sigma0 * c->inv_sqrt_n;

  if((ring_beads & 1) == 0) {
    DOUBLE nyq = z1 * c->sigma_nyq * c->inv_sqrt_n;
    DOUBLE even_val = dc + nyq;
    DOUBLE odd_val = dc - nyq;
    for(n = 0; n + 1 < ring_beads; n += 2) {
      ring_out[n] = even_val;
      ring_out[n + 1] = odd_val;
    }
  } else {
    for(n = 0; n < ring_beads; n++) {
      ring_out[n] = dc;
    }
  }
}

static void SampleHarmonicNormalRing1D(DOUBLE *ring_out, int ring_beads, DOUBLE tau, DOUBLE omega2) {
  int n, mode;
  int half_n = (ring_beads - 1) / 2;
  int is_even = ((ring_beads & 1) == 0);
  int half_idx = ring_beads >> 1;
  const NormalModeCache1D *c = GetNormalModeCache(ring_beads, tau, omega2);
  DOUBLE amp = c->amp_mode;

  InitHarmonicRingDCAndNyquist(ring_out, ring_beads, c);

  for(mode = 0; mode < c->n_modes; mode++) {
    DOUBLE sig = c->sigma_mode[mode];
    DOUBLE step_cosine = c->step_cos[mode];
    DOUBLE step_sine = c->step_sin[mode];
    DOUBLE q_cos, q_sin;
    DOUBLE cosine = step_cosine;
    DOUBLE sine = step_sine;

    RandomNormal(&q_cos, &q_sin, 0., sig);
    q_cos *= amp;
    q_sin *= amp;

    ring_out[0] += q_cos;
    if(is_even) {
      ring_out[half_idx] += (mode & 1) ? q_cos : -q_cos;
    }

    for(n = 1; n <= half_n; n++) {
      DOUBLE c_term = q_cos * cosine;
      DOUBLE s_term = q_sin * sine;
      DOUBLE angle;
      ring_out[n] += c_term + s_term;
      ring_out[ring_beads - n] += c_term - s_term;

      angle = cosine * step_cosine - sine * step_sine;
      sine = sine * step_cosine + cosine * step_sine;
      cosine = angle;
      if((n & 31) == 0) {
        DOUBLE inv_r = 0.5 * (3.0 - (cosine * cosine + sine * sine));
        cosine *= inv_r;
        sine *= inv_r;
      }
    }
  }
}

#if defined(MOVE_IN_X) && defined(MOVE_IN_Y)
static void SampleHarmonicNormalRing2D(DOUBLE *ring_x, DOUBLE *ring_y, int ring_beads, DOUBLE tau, DOUBLE omega2) {
  int n, mode;
  int half_n = (ring_beads - 1) / 2;
  int is_even = ((ring_beads & 1) == 0);
  int half_idx = ring_beads >> 1;
  const NormalModeCache1D *c = GetNormalModeCache(ring_beads, tau, omega2);
  DOUBLE amp = c->amp_mode;

  InitHarmonicRingDCAndNyquist(ring_x, ring_beads, c);
  InitHarmonicRingDCAndNyquist(ring_y, ring_beads, c);

  for(mode = 0; mode < c->n_modes; mode++) {
    DOUBLE sig = c->sigma_mode[mode];
    DOUBLE step_cosine = c->step_cos[mode];
    DOUBLE step_sine = c->step_sin[mode];
    DOUBLE qx_cos, qx_sin, qy_cos, qy_sin;
    DOUBLE cosine = step_cosine;
    DOUBLE sine = step_sine;

    RandomNormal(&qx_cos, &qx_sin, 0., sig);
    RandomNormal(&qy_cos, &qy_sin, 0., sig);
    qx_cos *= amp;
    qx_sin *= amp;
    qy_cos *= amp;
    qy_sin *= amp;

    ring_x[0] += qx_cos;
    ring_y[0] += qy_cos;
    if(is_even) {
      ring_x[half_idx] += (mode & 1) ? qx_cos : -qx_cos;
      ring_y[half_idx] += (mode & 1) ? qy_cos : -qy_cos;
    }

    for(n = 1; n <= half_n; n++) {
      DOUBLE cx_term = qx_cos * cosine;
      DOUBLE sx_term = qx_sin * sine;
      DOUBLE cy_term = qy_cos * cosine;
      DOUBLE sy_term = qy_sin * sine;
      DOUBLE angle;

      ring_x[n] += cx_term + sx_term;
      ring_x[ring_beads - n] += cx_term - sx_term;
      ring_y[n] += cy_term + sy_term;
      ring_y[ring_beads - n] += cy_term - sy_term;

      angle = cosine * step_cosine - sine * step_sine;
      sine = sine * step_cosine + cosine * step_sine;
      cosine = angle;
      if((n & 31) == 0) {
        DOUBLE inv_r = 0.5 * (3.0 - (cosine * cosine + sine * sine));
        cosine *= inv_r;
        sine *= inv_r;
      }
    }
  }
}
#endif

void PIMCMoveHarmonicNormalModes(void) {
  int bead, particle;
  int beads = Nwalkers;
  DOUBLE tau;

  if(beads > nm_plain_ring_cap) {
    DOUBLE *new_ring = (DOUBLE*) realloc(nm_plain_ring, (size_t)beads * sizeof(DOUBLE));
    if(new_ring == NULL) Error("Out of memory in PIMCMoveHarmonicNormalModes\n");
    nm_plain_ring = new_ring;
    nm_plain_ring_cap = beads;
  }

  if(N < 1 || beads < 1 || T <= 0.) {
    Error("Harmonic normal-mode sampling requires N>=1, Nwalkers>=1, and T>0\n");
  }
  CaseX(if(omega_x2 <= 0.) Error("Harmonic normal-mode sampling requires omega_x>0\n"));
  CaseY(if(omega_y2 <= 0.) Error("Harmonic normal-mode sampling requires omega_y>0\n"));
  CaseZ(if(omega_z2 <= 0.) Error("Harmonic normal-mode sampling requires omega_z>0\n"));

  tau = 1. / (T * (DOUBLE)beads);

#ifdef WORM_ALGORITHM
  if(worm_next != NULL) {
    if(!EnsureWormRingScratch(N, N * beads)) {
      Error("Out of memory in PIMCMoveHarmonicNormalModes worm scratch\n");
    }
    int *visited = worm_ring_visited;
    int *ring_p = worm_ring_p;
    int *ring_w = worm_ring_w;
    DOUBLE *ring_x = worm_ring_x;
    DOUBLE *ring_y = worm_ring_y;
    DOUBLE *ring_z = worm_ring_z;
    (void)ring_x;
    (void)ring_y;
    (void)ring_z;
    memset(visited, 0, (size_t)N * sizeof(int));

    for(particle = 0; particle < N; particle++) {
      if(visited[particle]) continue;
      if(worm_open != 0 && worm_next[0][particle] < 0) {
        visited[particle] = 1;
        continue;
      }

      /* Rastrear el anillo desde slice 0; en Sector G ignorar cualquier tramo abierto */
      int is_closed = 1;
      int ring_beads = 0;
      int curr_p = particle;
      int curr_w = 0;
      do {
        if(curr_w == 0) visited[curr_p] = 1;
        ring_p[ring_beads] = curr_p;
        ring_w[ring_beads] = curr_w;
        ring_beads++;
        curr_p = worm_next[curr_w][curr_p];
        if(curr_p < 0) {
          is_closed = 0;
          break;
        }
        curr_w = (curr_w + 1 < beads) ? curr_w + 1 : 0;
      } while(!(curr_w == 0 && curr_p == particle) && ring_beads < N * beads);

      if(!is_closed || !(curr_w == 0 && curr_p == particle)) continue;

#if defined(MOVE_IN_X) && defined(MOVE_IN_Y)
      if(omega_x2 == omega_y2) {
        SampleHarmonicNormalRing2D(ring_x, ring_y, ring_beads, tau, omega_x2);
      } else {
        SampleHarmonicNormalRing1D(ring_x, ring_beads, tau, omega_x2);
        SampleHarmonicNormalRing1D(ring_y, ring_beads, tau, omega_y2);
      }
#else
      CaseX(SampleHarmonicNormalRing1D(ring_x, ring_beads, tau, omega_x2));
      CaseY(SampleHarmonicNormalRing1D(ring_y, ring_beads, tau, omega_y2));
#endif
      CaseZ(SampleHarmonicNormalRing1D(ring_z, ring_beads, tau, omega_z2));

      /* Asignar las coordenadas calculadas a las cuentas físicas */
      for(int n = 0; n < ring_beads; n++) {
        int rw = ring_w[n];
        int rp = ring_p[n];
#ifdef MOVE_IN_X
        W[rw].x[rp] = ring_x[n];
#else
        W[rw].x[rp] = 0.;
#endif
#ifdef MOVE_IN_Y
        W[rw].y[rp] = ring_y[n];
#else
        W[rw].y[rp] = 0.;
#endif
#ifdef MOVE_IN_Z
        W[rw].z[rp] = ring_z[n];
#else
        W[rw].z[rp] = 0.;
#endif
      }
    }

    return;
  }
#endif

  for(particle = 0; particle < N; particle++) {
#ifdef MOVE_IN_X
    SampleHarmonicNormalRing1D(nm_plain_ring, beads, tau, omega_x2);
    for(bead = 0; bead < beads; bead++) W[bead].x[particle] = nm_plain_ring[bead];
#else
    for(bead = 0; bead < beads; bead++) W[bead].x[particle] = 0.;
#endif

#ifdef MOVE_IN_Y
    SampleHarmonicNormalRing1D(nm_plain_ring, beads, tau, omega_y2);
    for(bead = 0; bead < beads; bead++) W[bead].y[particle] = nm_plain_ring[bead];
#else
    for(bead = 0; bead < beads; bead++) W[bead].y[particle] = 0.;
#endif

#ifdef MOVE_IN_Z
    SampleHarmonicNormalRing1D(nm_plain_ring, beads, tau, omega_z2);
    for(bead = 0; bead < beads; bead++) W[bead].z[particle] = nm_plain_ring[bead];
#else
    for(bead = 0; bead < beads; bead++) W[bead].z[particle] = 0.;
#endif
  }
}
#endif

/************************** PIMC Move One By One ****************************/
static void PIMCMoveSingleBeadParticle(int w, int i, int w1, int w2,
                                       DOUBLE tau, DOUBLE inv_tau, DOUBLE inv_2tau, DOUBLE step_sigma,
                                       int *local_accepted, int *local_rejected) {
#ifndef INTERACTION_ABSENT
  int j;
#endif
  DOUBLE du = 0.;
  DOUBLE xi;
  DOUBLE x = 0., y = 0., z = 0.;
  DOUBLE xp = 0., yp = 0., zp = 0.;
  DOUBLE dr[3];
  DOUBLE xm, ym, zm;

  CaseX(x = W[w].x[i]);
  CaseY(y = W[w].y[i]);
  CaseZ(z = W[w].z[i]);
#ifdef WORM_ALGORITHM
  int prev_p = (worm_prev != NULL) ? worm_prev[w][i] : i;
  int next_p = (worm_next != NULL) ? worm_next[w][i] : i;
  if(prev_p < 0 && next_p < 0) return; /* detached bead in open-chain sector */
#else
  int prev_p = i;
  int next_p = i;
#endif

#ifdef TRIAL_3D
  {
    DOUBLE dx, dy, dz;
    RandomNormal(&dx, &dy, 0., step_sigma);
    RandomNormal(&dz, &xi, 0., step_sigma);
    xp = x + dx;
    yp = y + dy;
    zp = z + dz;
  }
#elif defined(TRIAL_2D)
  {
    DOUBLE dx, dy;
    RandomNormal(&dx, &dy, 0., step_sigma);
    xp = x + dx;
    yp = y + dy;
  }
#else
  {
    DOUBLE dx, dz;
    RandomNormal(&dz, &dx, 0., step_sigma);
    zp = z + dz;
  }
#endif

  ReduceToTheBoxXYZ(&xp, &yp, &zp);

  du += tau * (Vext(x, y, z) - Vext(xp, yp, zp));

  /* kinetic energy contribution (two links for interior beads, one link for Ira/Masha) */
  if(prev_p >= 0 && next_p >= 0) {
    DOUBLE r2_new, r2_old;
    xm = ym = zm = 0.;
    dr[0] = dr[1] = dr[2] = 0.;
    CaseX(xm = 0.5 * (W[w1].x[prev_p] + W[w2].x[next_p]); dr[0] = xp - xm);
    CaseY(ym = 0.5 * (W[w1].y[prev_p] + W[w2].y[next_p]); dr[1] = yp - ym);
    CaseZ(zm = 0.5 * (W[w1].z[prev_p] + W[w2].z[next_p]); dr[2] = zp - zm);
    r2_new = FindNearestImage(&dr[0], &dr[1], &dr[2]);

    CaseX(dr[0] = x - xm);
    CaseY(dr[1] = y - ym);
    CaseZ(dr[2] = z - zm);
    r2_old = FindNearestImage(&dr[0], &dr[1], &dr[2]);
    du -= (r2_new - r2_old) * inv_tau;
  } else {
    int w_conn = (prev_p >= 0) ? w1 : w2;
    int p_conn = (prev_p >= 0) ? prev_p : next_p;
    DOUBLE r2_new, r2_old;
    xm = ym = zm = 0.;
    dr[0] = dr[1] = dr[2] = 0.;
    CaseX(xm = W[w_conn].x[p_conn]; dr[0] = xp - xm);
    CaseY(ym = W[w_conn].y[p_conn]; dr[1] = yp - ym);
    CaseZ(zm = W[w_conn].z[p_conn]; dr[2] = zp - zm);
    r2_new = FindNearestImage(&dr[0], &dr[1], &dr[2]);

    CaseX(dr[0] = x - xm);
    CaseY(dr[1] = y - ym);
    CaseZ(dr[2] = z - zm);
    r2_old = FindNearestImage(&dr[0], &dr[1], &dr[2]);
    du -= (r2_new - r2_old) * inv_2tau;
  }

#ifndef INTERACTION_ABSENT
  for(j = 0; j < N; j++) {
#ifdef WORM_ALGORITHM
    if(j != i && (worm_prev == NULL || worm_next == NULL || worm_prev[w][j] >= 0 || worm_next[w][j] >= 0)) {
#else
    if(j != i) {
#endif
      DOUBLE r2;
      dr[0] = xp - W[w].x[j];
      dr[1] = yp - W[w].y[j];
      dr[2] = zp - W[w].z[j];
      r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
      du -= tau * InteractionEnergy(sqrt(r2));

      dr[0] = x - W[w].x[j];
      dr[1] = y - W[w].y[j];
      dr[2] = z - W[w].z[j];
      r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
      du += tau * InteractionEnergy(sqrt(r2));
    }
  }
#endif

  if(du > 0.) {
    CaseX(W[w].x[i] = xp);
    CaseY(W[w].y[i] = yp);
    CaseZ(W[w].z[i] = zp);
    (*local_accepted)++;
  } else {
    xi = Random();
    if(xi < Exp(du)) {
      CaseX(W[w].x[i] = xp);
      CaseY(W[w].y[i] = yp);
      CaseZ(W[w].z[i] = zp);
      (*local_accepted)++;
    } else {
      (*local_rejected)++;
    }
  }
}

/** @brief Perform a single PIMC move for one bead and return local accept/reject counts.
 *  @param w Index of the bead to be moved.
 *  @param accepted_out Output count of accepted single-particle moves.
 *  @param rejected_out Output count of rejected single-particle moves.
 */
static void PIMCMoveOneByOneBead(int w, int *accepted_out, int *rejected_out) {
  int i;
  int w1 = (w > 0) ? (w - 1) : (Nwalkers - 1);
  int w2 = (w + 1 < Nwalkers) ? (w + 1) : 0;
  DOUBLE tau = 1. / (T * (DOUBLE)Nwalkers);
  DOUBLE inv_tau = 1. / tau;
  DOUBLE inv_2tau = 0.5 * inv_tau;
  DOUBLE step_sigma = Sqrt(dt_vmc);
  int local_accepted = 0;
  int local_rejected = 0;

  for(i = 0; i < N; i++) {
    PIMCMoveSingleBeadParticle(w, i, w1, w2, tau, inv_tau, inv_2tau, step_sigma,
                               &local_accepted, &local_rejected);
  }

  *accepted_out = local_accepted;
  *rejected_out = local_rejected;
}

#ifdef WORM_ALGORITHM
/** @brief Perform single-bead Metropolis updates only along the active open worm chain in Sector G. */
void PIMCMoveOpenWormChain(void) {
  int curr_w, curr_p, max_steps, steps;
  int local_accepted = 0;
  int local_rejected = 0;
  DOUBLE tau, inv_tau, inv_2tau, step_sigma;

  if(Nwalkers <= 0 || N <= 0) return;
  if(worm_open == 0 || worm_next == NULL || worm_prev == NULL ||
     worm_bead_masha < 0 || worm_particle_masha < 0) {
    PIMCMoveOneByOneSweep();
    return;
  }

  tau = 1. / (T * (DOUBLE)Nwalkers);
  inv_tau = 1. / tau;
  inv_2tau = 0.5 * inv_tau;
  step_sigma = Sqrt(dt_vmc);

  curr_w = worm_bead_masha;
  curr_p = worm_particle_masha;
  max_steps = N * Nwalkers;

  for(steps = 0; steps < max_steps && curr_p >= 0; steps++) {
    int w1 = (curr_w > 0) ? (curr_w - 1) : (Nwalkers - 1);
    int w2 = (curr_w + 1 < Nwalkers) ? (curr_w + 1) : 0;
    int next_p = worm_next[curr_w][curr_p];
    PIMCMoveSingleBeadParticle(curr_w, curr_p, w1, w2, tau, inv_tau, inv_2tau, step_sigma,
                               &local_accepted, &local_rejected);
    curr_p = next_p;
    curr_w = w2;
  }

  accepted += local_accepted;
  rejected += local_rejected;
}
#endif

/** @brief Perform a single PIMC move for one bead.
 *  @param w Index of the bead to be moved.
 */
void PIMCMoveOneByOne(int w) {
  int local_accepted = 0;
  int local_rejected = 0;

  PIMCMoveOneByOneBead(w, &local_accepted, &local_rejected);

#pragma omp atomic
  accepted += local_accepted;
#pragma omp atomic
  rejected += local_rejected;
}

/** @brief Perform a race-free parallel PIMC sweep across all beads using checkerboard coloring. */
void PIMCMoveOneByOneSweep(void) {
  int w;
  int acc_sum = 0;
  int rej_sum = 0;

  if(Nwalkers <= 0) return;

  if(Nwalkers < 3) {
    for(w = 0; w < Nwalkers; w++) {
      int local_acc = 0, local_rej = 0;
      PIMCMoveOneByOneBead(w, &local_acc, &local_rej);
      acc_sum += local_acc;
      rej_sum += local_rej;
    }
    accepted += acc_sum;
    rejected += rej_sum;
    return;
  }

  /* For odd Nwalkers, bead Nwalkers-1 neighbors both 0 and Nwalkers-2, so isolate it in a 3rd step. */
  int count_even = (Nwalkers % 2 == 0) ? (Nwalkers / 2) : ((Nwalkers - 1) / 2);
  int count_odd = Nwalkers / 2;

#pragma omp parallel reduction(+:acc_sum, rej_sum)
  {
    int idx;
#pragma omp for
    for(idx = 0; idx < count_even; idx++) {
      int local_acc = 0, local_rej = 0;
      PIMCMoveOneByOneBead(2 * idx, &local_acc, &local_rej);
      acc_sum += local_acc;
      rej_sum += local_rej;
    }
#pragma omp for
    for(idx = 0; idx < count_odd; idx++) {
      int local_acc = 0, local_rej = 0;
      PIMCMoveOneByOneBead(2 * idx + 1, &local_acc, &local_rej);
      acc_sum += local_acc;
      rej_sum += local_rej;
    }
  }

  if(Nwalkers % 2 == 1) {
    int local_acc = 0, local_rej = 0;
    PIMCMoveOneByOneBead(Nwalkers - 1, &local_acc, &local_rej);
    acc_sum += local_acc;
    rej_sum += local_rej;
  }

  accepted += acc_sum;
  rejected += rej_sum;
}


/*************************** PIGS pseudopotential *************************************/
/** @brief Perform a single PIMC move for the pseudopotential.
 *  @param w Index of the walker to be moved.
 */
void PIMCPseudopotentialMoveOneByOne(int w) {
  int i,j;
  int w1,w2; // index of the walker with [w-1] and [w+1]
  DOUBLE dS, dSkin, dSdelta, dSwf;
  DOUBLE dx,dy,dz;
  DOUBLE tau;
  int chain_edge = OFF;
  int up;
  DOUBLE x, y, z;
  DOUBLE xm, ym, zm; // position of the middle point
  //DOUBLE mA; // i.e. 2 mu
  //DOUBLE aA;
#ifdef EXTERNAL_POTENTIAL
  DOUBLE dEpot;
#endif
#ifdef SECURE
  DOUBLE dS_check;
  DOUBLE dSkin_check_old, dSpot_check_old, dSwf_check_old;
  DOUBLE dSkin_check_new, dSpot_check_new, dSwf_check_new;
#endif

  x = y = z = 0.;
  xm = ym = zm = 0.;

  //tau = dt; // imaginary time step
  tau = beta;

  w1 = w-1;
  w2 = w+1;
  if(w == 0) { // treat w=0 and w=Nwalkers-1 separately, according to the cyclic condition
    w1 = Nwalkers-1;
  }
  if(w == Nwalkers-1) { // treat w=0 and w=Nwalkers-1 separately
    w2 = 0;
  }

  // Move particle i of the bead
  for(i=0; i<N; i++) { // move i-th particle to the point rm
    dSkin = dSdelta = dSwf = 0.;
#ifdef EXTERNAL_POTENTIAL
    dEpot = 0.;
#endif

    CaseZ(Wp[w].z[i] = W[w].z[i]); // store coordinates

#ifdef SECURE
   dS_check = -PathIntegralPseudopotentialSdetailed(&dSkin_check_old, &dSpot_check_old, &dSwf_check_old);
#endif

    // Wp is the trial position
    RandomNormal3(&dx, &dy, &dz, 0., Sqrt(dt));
    W[w].z[i] += dz;
#ifdef BC_1DPBC
    ReduceToTheBoxXYZ(&W[w].x[i], &W[w].y[i], &W[w].z[i]);
#endif

  // contribution from the Laplacian, - (R'(i)-R(i+/-1))^2 / 2 + (R(i)-R(i+/-1))^2 / 2
#ifndef INFINITE_MASS_COMPONENT_UP // mA finite
    dSkin += PropagatorLaplacian(mA, tau, W[w].z[i] - W[w1].z[i]);
    dSkin += PropagatorLaplacian(mA, tau, W[w].z[i] - W[w2].z[i]);
    dSkin -= PropagatorLaplacian(mA, tau, Wp[w].z[i] - W[w1].z[i]);
    dSkin -= PropagatorLaplacian(mA, tau, Wp[w].z[i] - W[w2].z[i]);
#endif
#ifdef EXTERNAL_POTENTIAL
  dEpot += VextUp(W[w].x[i], W[w].y[i], W[w].z[i], i); // new external potential energy
  dEpot -= VextUp(Wp[w].z[i],Wp[w].y[i],Wp[w].z[i],i); // old external potential energy
#endif

    // contribution from Pseudopotential
    for(j=0; j<N; j++) if(j != i) { // interaction energy: AA
      dSdelta += PropagatorPseudopotenial1D_ij(mA, aA, tau, W[w].z[i],  W[w].z[j], W[w1].z[i], W[w1].z[j]);
      dSdelta += PropagatorPseudopotenial1D_ij(mA, aA, tau, W[w].z[i],  W[w].z[j], W[w2].z[i], W[w2].z[j]);
      dSdelta -= PropagatorPseudopotenial1D_ij(mA, aA, tau, Wp[w].z[i], W[w].z[j], W[w1].z[i], W[w1].z[j]);
      dSdelta -= PropagatorPseudopotenial1D_ij(mA, aA, tau, Wp[w].z[i], W[w].z[j], W[w2].z[i], W[w2].z[j]);
    }

#ifdef SECURE
    if(isnan(dSdelta)) Warning("NaN weight in PIGS move!\n");
#endif

    dS = dSkin + dSdelta + dSwf; // dSdelta - contribution from delta function, dSwf from the trial w.f.
#ifdef EXTERNAL_POTENTIAL
    dS -= tau*dEpot; // potential energy contribution
#endif
#ifdef SECURE
   dS_check += PathIntegralPseudopotentialSdetailed(&dSkin_check_new, &dSpot_check_new, &dSwf_check_new);
   if(fabs((dS-dS_check)/dS)>1e-5) {
     Message(" dS check (%lf %lf) mismatch dSkin %lf %lf, dSpot %lf %lf\n", dS, dS_check, dSkin, dSkin_check_new-dSkin_check_old, dSdelta, dSpot_check_new-dSpot_check_old);
     SaveCoordinates("check.dat");
   }
#endif

    if(Random() < Exp(dS)) { // The Metropolis code
      W[w].U -= dS;
      accepted_one++;
    }
    else {
      CaseX(W[w].x[i] = Wp[w].x[i]);
      CaseY(W[w].y[i] = Wp[w].y[i]);
      CaseZ(W[w].z[i] = Wp[w].z[i]);
      rejected_one++;
    }
  }
}


#ifdef WORM_ALGORITHM
/***********************************************************************
 *  Algoritmo Worm (Prokof'ev-Svistunov / Dornheim 2014)
 *
 *  Movimientos del Worm en el sector G
 *
 *  - Enlaces:
 *    worm_next[w][i]: índice de partícula en el bead (w+1)%Nwalkers.
 *    worm_prev[w][i]: índice de partícula en el bead (w-1+Nwalkers)%Nwalkers.
 *    Si el enlace está abierto (Ira o Masha), vale -1.
 *
 *  - Sector Z (worm_open == 0):
 *    Todos los enlaces están cerrados.
 *    Movimiento: Open (corta un enlace aleatorio, creando Ira y Masha con gap = 1).
 *
 *  - Sector G (worm_open == 1):
 *    Ira en (worm_bead, worm_particle), Masha en (worm_bead_masha, worm_particle_masha).
 *    Gap: delta_w = (worm_bead_masha - worm_bead + Nwalkers) % Nwalkers.
 *    Movimientos:
 *      Advance/Recede cambian el gap en +/-1; Close reconecta cuando gap == 1.
 *      Para N>1 también se proponen intercambios Swap antes de Masha.
 *      IraSpatialMove desplaza espacialmente la cabeza Ira.
 ***********************************************************************/

static int worm_snapshot_seen[6] = {0, 0, 0, 0, 0, 0};
static const char *worm_snapshot_names[6] = {
  "Open", "Close", "Advance", "Recede", "Swap", "IraSpatialMove"
};

static int *worm_swap_candidates = NULL;
static DOUBLE *worm_swap_weights = NULL;

static DOUBLE worm_obdm_inv_step_x = 0.0;
static DOUBLE worm_obdm_inv_step_r = 0.0;
static DOUBLE worm_obdm_inv_step_s = 0.0;
static DOUBLE worm_obdm_max2 = 0.0;
static DOUBLE worm_obdm_smax2 = 0.0;
static int worm_obdm_geom_grid = -1;
static DOUBLE worm_obdm_geom_max = -1.0;

static void SyncWormOBDMGeometry(void) {
  if(worm_obdm_geom_grid == worm_obdm_grid && worm_obdm_geom_max == worm_obdm_max) {
    return;
  }
  worm_obdm_inv_step_x = (DOUBLE)worm_obdm_grid / (2.0 * worm_obdm_max);
  worm_obdm_inv_step_r = (DOUBLE)worm_obdm_grid / worm_obdm_max;
  worm_obdm_inv_step_s = (DOUBLE)worm_obdm_grid / (2.0 * worm_obdm_max);
  worm_obdm_max2 = worm_obdm_max * worm_obdm_max;
  worm_obdm_smax2 = 4.0 * worm_obdm_max2;
  worm_obdm_geom_grid = worm_obdm_grid;
  worm_obdm_geom_max = worm_obdm_max;
}

static int** AllocContiguous2DInt(int rows, int cols, const char *name) {
  int **arr = (int**) malloc((size_t)rows * sizeof(int*));
  int *data;
  int r;
  if(arr == NULL) Error("Out of memory allocating %s pointers\n", name);
  data = (int*) calloc((size_t)rows * (size_t)cols, sizeof(int));
  if(data == NULL) {
    free(arr);
    Error("Out of memory allocating %s data\n", name);
  }
  for(r = 0; r < rows; r++) {
    arr[r] = data + (size_t)r * (size_t)cols;
  }
  return arr;
}

static void FreeContiguous2DInt(int ***arr) {
  if(*arr != NULL) {
    free((*arr)[0]);
    free(*arr);
    *arr = NULL;
  }
}

static long long** AllocContiguous2DLL(int rows, int cols, const char *name) {
  long long **arr = (long long**) malloc((size_t)rows * sizeof(long long*));
  long long *data;
  int r;
  if(arr == NULL) Error("Out of memory allocating %s pointers\n", name);
  data = (long long*) calloc((size_t)rows * (size_t)cols, sizeof(long long));
  if(data == NULL) {
    free(arr);
    Error("Out of memory allocating %s data\n", name);
  }
  for(r = 0; r < rows; r++) {
    arr[r] = data + (size_t)r * (size_t)cols;
  }
  return arr;
}

static void FreeContiguous2DLL(long long ***arr) {
  if(*arr != NULL) {
    free((*arr)[0]);
    free(*arr);
    *arr = NULL;
  }
}

void FreeWorm(void) {
  if(worm_ring_visited != NULL) {
    free(worm_ring_visited);
    worm_ring_visited = NULL;
  }
  if(worm_ring_p != NULL) {
    free(worm_ring_p);
    worm_ring_p = NULL;
  }
  if(worm_ring_w != NULL) {
    free(worm_ring_w);
    worm_ring_w = NULL;
  }
  if(worm_ring_x != NULL) {
    free(worm_ring_x);
    worm_ring_x = NULL;
  }
  if(worm_ring_y != NULL) {
    free(worm_ring_y);
    worm_ring_y = NULL;
  }
  if(worm_ring_z != NULL) {
    free(worm_ring_z);
    worm_ring_z = NULL;
  }
  worm_ring_cap_n = 0;
  worm_ring_cap_total = 0;
  if(worm_swap_candidates != NULL) {
    free(worm_swap_candidates);
    worm_swap_candidates = NULL;
  }
  if(worm_swap_weights != NULL) {
    free(worm_swap_weights);
    worm_swap_weights = NULL;
  }
  if(worm_gap_histogram != NULL) {
    free(worm_gap_histogram);
    worm_gap_histogram = NULL;
  }
  if(worm_cycle_histogram != NULL) {
    free(worm_cycle_histogram);
    worm_cycle_histogram = NULL;
  }
  FreeContiguous2DLL(&worm_obdm_matrix);
  FreeContiguous2DLL(&worm_obdm_radial);
  FreeContiguous2DLL(&worm_obdm_matrix_all);
  worm_obdm_geom_grid = -1;
  worm_obdm_geom_max = -1.0;
  if(worm_obdm_dist != NULL) {
    free(worm_obdm_dist);
    worm_obdm_dist = NULL;
  }
  if(worm_obdm_dist_all != NULL) {
    free(worm_obdm_dist_all);
    worm_obdm_dist_all = NULL;
  }
  FreeContiguous2DInt(&worm_next);
  FreeContiguous2DInt(&worm_prev);
}

void InitWorm(void) {
  int w, i;
  if(worm_next != NULL) return;

  if(getenv("WORM_DISABLE_SWAP") != NULL && atoi(getenv("WORM_DISABLE_SWAP")) != 0) {
    worm_swap_enabled = 0;
  }
  if(getenv("WORM_SAVE_SNAPSHOTS") != NULL) {
    worm_save_snapshots = (atoi(getenv("WORM_SAVE_SNAPSHOTS")) != 0);
  }

  atexit(FreeWorm);

  {
    int move;
    for(move = 0; move < 6; move++) worm_snapshot_seen[move] = 0;
    if(worm_save_snapshots) {
      FILE *snapshot_file = fopen("outworm_snapshots.dat", "w");
      if(snapshot_file != NULL) {
        fprintf(snapshot_file, "# worm_snapshot_v1\n");
        fclose(snapshot_file);
      }
    }
  }

  worm_next = AllocContiguous2DInt(NwalkersMax, N, "worm_next");
  worm_prev = AllocContiguous2DInt(NwalkersMax, N, "worm_prev");
  for(w = 0; w < NwalkersMax; w++) {
    for(i = 0; i < N; i++) {
      worm_next[w][i] = i;
      worm_prev[w][i] = i;
    }
  }

  worm_swap_candidates = (int*) malloc((size_t)N * sizeof(int));
  worm_swap_weights = (DOUBLE*) malloc((size_t)N * sizeof(DOUBLE));
  if(worm_swap_candidates == NULL || worm_swap_weights == NULL) {
    Error("Out of memory allocating worm_swap scratch\n");
  }

  worm_gap_histogram = (long long*) calloc((size_t)(NwalkersMax + 1), sizeof(long long));
  worm_max_gap_observed = 0;

  worm_cycle_histogram = (long long*) calloc((size_t)(N + 1), sizeof(long long));
  worm_cycle_samples = 0;

  if(worm_obdm_grid < 10) worm_obdm_grid = 100;

  worm_obdm_matrix = AllocContiguous2DLL(worm_obdm_grid, worm_obdm_grid, "worm_obdm_matrix");
  worm_obdm_radial = AllocContiguous2DLL(worm_obdm_grid, worm_obdm_grid, "worm_obdm_radial");
  worm_obdm_matrix_all = AllocContiguous2DLL(worm_obdm_grid, worm_obdm_grid, "worm_obdm_matrix_all");
  worm_obdm_dist = (long long*) calloc((size_t)worm_obdm_grid, sizeof(long long));
  worm_obdm_dist_all = (long long*) calloc((size_t)worm_obdm_grid, sizeof(long long));
  if(worm_gap_histogram == NULL || worm_cycle_histogram == NULL ||
     worm_obdm_dist == NULL || worm_obdm_dist_all == NULL) {
    Error("Out of memory allocating worm histograms\n");
  }
  SyncWormOBDMGeometry();
  worm_obdm_samples_gap1 = 0;
  worm_obdm_samples_all = 0;

  worm_sf_samples = 0;
  worm_sf_block_samples = 0;
  worm_sf_num_blocks = 0;
  worm_sf_rho_block_mean = 0.0;
  worm_sf_rho_block_m2 = 0.0;
  worm_sf_block_weight_sum = 0.0;
  worm_sf_block_weight_sum2 = 0.0;
  worm_sf_W2_sum = 0.0;
  worm_sf_block_W2_sum = 0.0;
  for(i = 0; i < 3; i++) {
    worm_sf_A_sum[i] = 0.0;
    worm_sf_A2_sum[i] = 0.0;
    worm_sf_Icl_sum[i] = 0.0;
    worm_sf_block_A_sum[i] = 0.0;
    worm_sf_block_A2_sum[i] = 0.0;
    worm_sf_block_Icl_sum[i] = 0.0;
  }
}

static int WormSnapshotWrite(int move, const char *phase) {
  FILE *snapshot_file;
  int bead, particle;

  snapshot_file = fopen("outworm_snapshots.dat", "a");
  if(snapshot_file == NULL) return 0;

  fprintf(snapshot_file, "EVENT %s %s %d %d %d %d %d %d %d\n",
          worm_snapshot_names[move], phase, Nwalkers, N, worm_open,
          worm_bead, worm_particle, worm_bead_masha, worm_particle_masha);
  for(bead = 0; bead < Nwalkers; bead++) {
    for(particle = 0; particle < N; particle++) {
#if defined(TRIAL_1D) && !defined(BC_1DPBC_X)
      DOUBLE coord_primary = W[bead].z[particle];
#else
      DOUBLE coord_primary = W[bead].x[particle];
#endif
      fprintf(snapshot_file, "BEAD %d %d %.17g %.17g %d %d\n",
              bead, particle, (double)coord_primary, (double)W[bead].y[particle],
              worm_next[bead][particle], worm_prev[bead][particle]);
    }
  }
  fprintf(snapshot_file, "END\n");
  fclose(snapshot_file);
  return 1;
}

static int WormSnapshotBefore(int move) {
  if(!worm_save_snapshots || move < 0 || move >= 6 || worm_snapshot_seen[move]) return 0;
  return WormSnapshotWrite(move, "before");
}

static void WormSnapshotAfter(int move) {
  if(!worm_save_snapshots) return;
  if(WormSnapshotWrite(move, "after")) worm_snapshot_seen[move] = 1;
}

#define WORM_MAX_GAP 6
#define WORM_PROB_CLOSE_GAP1     0.60
#define WORM_PROB_RECEDE_GAP1    0.25
#define WORM_PROB_SPATIAL_GAP1   0.15
#define WORM_PROB_ADVANCE_GAP_GT1 0.45
#define WORM_PROB_RECEDE_GAP_GT1  0.45
#define WORM_PROB_SPATIAL_GAP_GT1 0.10

#define WORM_PROB_CLOSE_N_GT1    0.40
#define WORM_PROB_SWAP_N_GT1     0.40
#define WORM_PROB_SPATIAL_N_GT1  0.20
#define WORM_PROB_ADVANCE_N_GT1  0.20
#define WORM_PROB_RECEDE_N_GT1   0.20
#define WORM_PROB_RECEDE_N_GAP1  0.20
#define WORM_PROB_SWAP_N_GAP1    0.20
#define WORM_PROB_SPATIAL_N_GAP1 0.20

INLINE DOUBLE WormSegmentDistance2(int w1, int p1, int w2, int p2) {
  DOUBLE dx = W[w2].x[p2] - W[w1].x[p1];
  DOUBLE dy = W[w2].y[p2] - W[w1].y[p1];
  DOUBLE dz = W[w2].z[p2] - W[w1].z[p1];
  return FindNearestImage(&dx, &dy, &dz);
}

INLINE void WormStepPosition(DOUBLE x, DOUBLE y, DOUBLE z,
                             DOUBLE dx, DOUBLE dy, DOUBLE dz,
                             DOUBLE *xp, DOUBLE *yp, DOUBLE *zp) {
  (void)dx;
  (void)dy;
  (void)dz;
  *xp = x;
  *yp = y;
  *zp = z;
  CaseX(*xp += dx);
  CaseY(*yp += dy);
  CaseZ(*zp += dz);
  ReduceToTheBoxXYZ(xp, yp, zp);
}

static void WormOpenMove(void) {
  int w_cut, i_cut, j_cut, w_next;
  DOUBLE r2, tau_local, acc, p_close;

  worm_try_open++;

  w_cut = (int)(Random() * Nwalkers);
  if(w_cut >= Nwalkers) w_cut = Nwalkers - 1;
  i_cut = (int)(Random() * N);
  if(i_cut >= N) i_cut = N - 1;

  w_next = (w_cut + 1) % Nwalkers;
  j_cut = worm_next[w_cut][i_cut];
  if(j_cut < 0) return; /* ya abierto */

  r2 = WormSegmentDistance2(w_cut, i_cut, w_next, j_cut);

  tau_local = 1.0 / (T * (DOUBLE)Nwalkers);
  p_close = (N > 1 && worm_swap_enabled) ? WORM_PROB_CLOSE_N_GT1 : WORM_PROB_CLOSE_GAP1;

  /* Aceptación Metropolis: ganancia del resorte cinético al cortar */
  acc = worm_C * Exp(r2 / (2.0 * tau_local))
    * p_close * (DOUBLE)Nwalkers * (DOUBLE)N;

  if(acc >= 1.0 || Random() < acc) {
    int record_snapshot = WormSnapshotBefore(0);
    worm_open = 1;
    worm_bead = w_cut;
    worm_particle = i_cut;
    worm.w = w_cut;
    worm.i = i_cut;
    worm.x_ira = W[w_cut].x[i_cut];
    worm.y_ira = W[w_cut].y[i_cut];
    worm.z_ira = W[w_cut].z[i_cut];

    worm_bead_masha = w_next;
    worm_particle_masha = j_cut;
    worm.x_masha = W[w_next].x[j_cut];
    worm.y_masha = W[w_next].y[j_cut];
    worm.z_masha = W[w_next].z[j_cut];

    worm_next[w_cut][i_cut] = -1;
    worm_prev[w_next][j_cut] = -1;

    worm_acc_open++;
    if(record_snapshot) WormSnapshotAfter(0);
  }
}

static void WormCloseMove(void) {
  int w_masha_expected;
  DOUBLE r2, tau_local, acc, p_close;

  worm_try_close++;

  w_masha_expected = (worm_bead + 1) % Nwalkers;
  if(worm_bead_masha != w_masha_expected) {
    return; /* solo se puede cerrar si los beads son contiguos (gap == 1) */
  }

  r2 = WormSegmentDistance2(worm_bead, worm_particle, worm_bead_masha, worm_particle_masha);

  tau_local = 1.0 / (T * (DOUBLE)Nwalkers);
  p_close = (N > 1 && worm_swap_enabled) ? WORM_PROB_CLOSE_N_GT1 : WORM_PROB_CLOSE_GAP1;

  /* Aceptación Metropolis: coste del resorte cinético al reconectar */
  acc = (1.0 / worm_C) * Exp(-r2 / (2.0 * tau_local))
    / (p_close * (DOUBLE)Nwalkers * (DOUBLE)N);

  if(acc >= 1.0 || Random() < acc) {
    int record_snapshot = WormSnapshotBefore(1);
    worm_next[worm_bead][worm_particle] = worm_particle_masha;
    worm_prev[worm_bead_masha][worm_particle_masha] = worm_particle;

    worm_open = 0;
    worm_bead = -1;
    worm_particle = -1;
    worm_bead_masha = -1;
    worm_particle_masha = -1;
    worm.w = -1;
    worm.i = -1;
    worm.x_ira = 0.;
    worm.y_ira = 0.;
    worm.z_ira = 0.;
    worm.x_masha = 0.;
    worm.y_masha = 0.;
    worm.z_masha = 0.;

    worm_acc_close++;
    if(record_snapshot) WormSnapshotAfter(1);
  }
}

static void WormAdvanceMove(void) {
  int w_curr, w_next, p_idx, p_next, gap, particle, candidate_count;
  DOUBLE dx = 0., dy = 0., dz = 0., x_new, y_new, z_new, dV, acc, tau_local, proposal_ratio;
  DOUBLE p_advance, p_recede_reverse;

  worm_try_advance++;

  w_curr = worm_bead;
  p_idx = worm_particle;
  w_next = (w_curr + 1) % Nwalkers;

  /* No avanzar si ya estamos adyacentes a Masha (eso corresponde a Close) */
  if(w_next == worm_bead_masha) return;

  gap = (worm_bead_masha - worm_bead + Nwalkers) % Nwalkers;

  /* Recede leaves one detached bead on the slice immediately after Ira. */
  p_next = -1;
  candidate_count = 0;
  for(particle = 0; particle < N; particle++) {
    if(worm_prev[w_next][particle] < 0 && worm_next[w_next][particle] < 0) {
      p_next = particle;
      candidate_count++;
      if(candidate_count > 1) return;
    }
  }
  if(candidate_count != 1) return;

  tau_local = 1.0 / (T * (DOUBLE)Nwalkers);

  /* Muestreo gaussiano del propagador cinético libre desde Ira */
  RandomNormal3(&dx, &dy, &dz, 0., Sqrt(tau_local));
  WormStepPosition(W[w_curr].x[p_idx], W[w_curr].y[p_idx], W[w_curr].z[p_idx],
                   dx, dy, dz, &x_new, &y_new, &z_new);

  /* Cambio en la acción potencial en el nuevo bead */
  dV = Vext(x_new, y_new, z_new);
  p_advance = (N > 1 && worm_swap_enabled) ? WORM_PROB_ADVANCE_N_GT1 : WORM_PROB_ADVANCE_GAP_GT1;
  if(gap == 2) {
    p_recede_reverse = (N > 1 && worm_swap_enabled) ? WORM_PROB_RECEDE_N_GAP1 : WORM_PROB_RECEDE_GAP1;
  } else {
    p_recede_reverse = (N > 1 && worm_swap_enabled) ? WORM_PROB_RECEDE_N_GT1 : WORM_PROB_RECEDE_GAP_GT1;
  }
  proposal_ratio = p_recede_reverse / p_advance;
  acc = Exp(-tau_local * dV) * proposal_ratio;

  if(acc >= 1.0 || Random() < acc) {
    int record_snapshot = WormSnapshotBefore(2);
    W[w_next].x[p_next] = x_new;
    W[w_next].y[p_next] = y_new;
    W[w_next].z[p_next] = z_new;

    worm_next[w_curr][p_idx] = p_next;
    worm_prev[w_next][p_next] = p_idx;

    worm_bead = w_next;
    worm_particle = p_next;
    worm.w = w_next;
    worm.i = p_next;
    worm.x_ira = x_new;
    worm.y_ira = y_new;
    worm.z_ira = z_new;

    worm_next[w_next][p_next] = -1;

    worm_acc_advance++;
    if(record_snapshot) WormSnapshotAfter(2);
  }
}

static void WormRecedeMove(void) {
  int w_curr, w_prev, p_idx, p_prev, gap;
  DOUBLE dV, acc, tau_local, proposal_ratio;
  DOUBLE p_recede, p_advance_reverse;

  worm_try_recede++;

  w_curr = worm_bead;
  p_idx = worm_particle;
  w_prev = (w_curr - 1 + Nwalkers) % Nwalkers;

  gap = (worm_bead_masha - worm_bead + Nwalkers) % Nwalkers;
    p_prev = worm_prev[w_curr][p_idx];
  /* Límite máximo de apertura (M_bar en Dornheim) y protección contra colisión con Masha */
    if(gap >= WORM_MAX_GAP || w_prev == worm_bead_masha || p_prev < 0
      || worm_next[w_prev][p_prev] != p_idx) return;

  tau_local = 1.0 / (T * (DOUBLE)Nwalkers);

  /* Al borrar el bead de Ira, se elimina su potencial */
  dV = Vext(W[w_curr].x[p_idx], W[w_curr].y[p_idx], W[w_curr].z[p_idx]);
  p_advance_reverse = (N > 1 && worm_swap_enabled) ? WORM_PROB_ADVANCE_N_GT1 : WORM_PROB_ADVANCE_GAP_GT1;
  if(gap == 1) {
    p_recede = (N > 1 && worm_swap_enabled) ? WORM_PROB_RECEDE_N_GAP1 : WORM_PROB_RECEDE_GAP1;
  } else {
    p_recede = (N > 1 && worm_swap_enabled) ? WORM_PROB_RECEDE_N_GT1 : WORM_PROB_RECEDE_GAP_GT1;
  }
  proposal_ratio = p_advance_reverse / p_recede;
  acc = Exp(+tau_local * dV) * proposal_ratio;

  if(acc >= 1.0 || Random() < acc) {
    int record_snapshot = WormSnapshotBefore(3);
    worm_next[w_prev][p_prev] = -1;
    worm_prev[w_curr][p_idx] = -1;

    worm_bead = w_prev;
    worm_particle = p_prev;
    worm.w = w_prev;
    worm.i = p_prev;
    worm.x_ira = W[w_prev].x[p_prev];
    worm.y_ira = W[w_prev].y[p_prev];
    worm.z_ira = W[w_prev].z[p_prev];

    worm_acc_recede++;
    if(record_snapshot) WormSnapshotAfter(3);
  }
}

static void WormSwapMove(void) {
  int w_curr, w_next, i_ira, k, p_old, cand_count;
  int *candidates = worm_swap_candidates;
  DOUBLE *weights = worm_swap_weights;
  const int *prev_next;
  const DOUBLE *wx_next, *wy_next, *wz_next;
  DOUBLE x_ira = 0., y_ira = 0., z_ira = 0.;
  DOUBLE x_old = 0., y_old = 0., z_old = 0.;
  DOUBLE sum_weights1, sum_weights2, r, cum, acc;
  DOUBLE inv_2tau, d2;
  int chosen_idx, chosen_k;

  worm_try_swap++;

  if(!worm_swap_enabled || N < 2) return; /* Se requieren swap activo y al menos 2 partículas */
  if(candidates == NULL || weights == NULL) return;

  w_curr = worm_bead;
  i_ira = worm_particle;
  w_next = (w_curr + 1 < Nwalkers) ? (w_curr + 1) : 0;
  inv_2tau = 0.5 * T * (DOUBLE)Nwalkers;

  prev_next = worm_prev[w_next];
  wx_next = W[w_next].x;
  wy_next = W[w_next].y;
  wz_next = W[w_next].z;
  CaseX(x_ira = W[w_curr].x[i_ira]);
  CaseY(y_ira = W[w_curr].y[i_ira]);
  CaseZ(z_ira = W[w_curr].z[i_ira]);
  (void)x_ira; (void)y_ira; (void)z_ira;
  (void)wx_next; (void)wy_next; (void)wz_next;

  /* 1 & 2. Filtrar candidatos con enlace entrante válido y acumular chi_1(k) en una sola pasada */
  cand_count = 0;
  sum_weights1 = 0.0;
  for(k = 0; k < N; k++) {
    if(prev_next[k] >= 0) {
      DOUBLE dx = 0., dy = 0., dz = 0.;
      DOUBLE w_val;
      CaseX(dx = wx_next[k] - x_ira);
      CaseY(dy = wy_next[k] - y_ira);
      CaseZ(dz = wz_next[k] - z_ira);
      d2 = FindNearestImage(&dx, &dy, &dz);
      w_val = Exp(-d2 * inv_2tau);
      candidates[cand_count] = k;
      weights[cand_count] = w_val;
      sum_weights1 += w_val;
      cand_count++;
    }
  }

  if(cand_count == 0 || sum_weights1 <= 0.0) return;

  /* 3. Elegir candidato k con probabilidad P(k) = chi_1(k) / Sigma_1 */
  r = Random() * sum_weights1;
  cum = 0.0;
  chosen_idx = cand_count - 1;
  for(int c = 0; c < cand_count - 1; c++) {
    cum += weights[c];
    if(r <= cum) {
      chosen_idx = c;
      break;
    }
  }
  chosen_k = candidates[chosen_idx];
  p_old = prev_next[chosen_k];

  /* 4. Calcular Sigma_2 para el movimiento inverso desde la nueva cabeza Ira (w_curr, p_old) */
  CaseX(x_old = W[w_curr].x[p_old]);
  CaseY(y_old = W[w_curr].y[p_old]);
  CaseZ(z_old = W[w_curr].z[p_old]);
  sum_weights2 = 0.0;
  for(int c = 0; c < cand_count; c++) {
    DOUBLE dx = 0., dy = 0., dz = 0.;
    k = candidates[c];
    CaseX(dx = wx_next[k] - x_old);
    CaseY(dy = wy_next[k] - y_old);
    CaseZ(dz = wz_next[k] - z_old);
    d2 = FindNearestImage(&dx, &dy, &dz);
    sum_weights2 += Exp(-d2 * inv_2tau);
  }

  if(sum_weights2 <= 0.0) return;

  /* 5. Aceptación Metropolis: A = min(1, Sigma_1 / Sigma_2) */
  acc = sum_weights1 / sum_weights2;

  if(acc >= 1.0 || Random() < acc) {
    int record_snapshot = WormSnapshotBefore(4);
    /* Reconectar enlaces */
    worm_next[w_curr][i_ira] = chosen_k;
    worm_prev[w_next][chosen_k] = i_ira;
    worm_next[w_curr][p_old] = -1;

    /* Actualizar cabeza Ira */
    worm_bead = w_curr;
    worm_particle = p_old;
    worm.w = w_curr;
    worm.i = p_old;
    worm.x_ira = x_old;
    worm.y_ira = y_old;
    worm.z_ira = z_old;

    worm_acc_swap++;
    if(record_snapshot) WormSnapshotAfter(4);
  }
}

static void WormIraSpatialMove(void) {
  int w_curr, w_prev, p_idx;
  DOUBLE dx = 0., dy = 0., dz = 0., x_new, y_new, z_new, dV_kin, dV_ext, dS, acc, tau_local;
  DOUBLE dx_old, dy_old, dz_old, dx_new, dy_new, dz_new, d2_old, d2_new;

  worm_try_iramove++;

  w_curr = worm_bead;
  p_idx = worm_particle;
  w_prev = (w_curr - 1 + Nwalkers) % Nwalkers;

  tau_local = 1.0 / (T * (DOUBLE)Nwalkers);

  /* Desplazamiento gaussiano local */
  RandomNormal3(&dx, &dy, &dz, 0., Sqrt(dt_vmc > 0 ? dt_vmc : tau_local));
  WormStepPosition(W[w_curr].x[p_idx], W[w_curr].y[p_idx], W[w_curr].z[p_idx],
                   dx, dy, dz, &x_new, &y_new, &z_new);

  /* Cambio cinético con el bead anterior */
  int p_in = (worm_prev != NULL && worm_prev[w_curr][p_idx] >= 0) ? worm_prev[w_curr][p_idx] : p_idx;
  dx_old = W[w_curr].x[p_idx] - W[w_prev].x[p_in];
  dy_old = W[w_curr].y[p_idx] - W[w_prev].y[p_in];
  dz_old = W[w_curr].z[p_idx] - W[w_prev].z[p_in];
  d2_old = FindNearestImage(&dx_old, &dy_old, &dz_old);

  dx_new = x_new - W[w_prev].x[p_in];
  dy_new = y_new - W[w_prev].y[p_in];
  dz_new = z_new - W[w_prev].z[p_in];
  d2_new = FindNearestImage(&dx_new, &dy_new, &dz_new);

  dV_kin = (d2_new - d2_old) / (2.0 * tau_local);
  dV_ext = tau_local * (Vext(x_new, y_new, z_new)
                        - Vext(W[w_curr].x[p_idx], W[w_curr].y[p_idx], W[w_curr].z[p_idx]));

  dS = -(dV_kin + dV_ext);
  acc = Exp(dS);

  if(acc >= 1.0 || Random() < acc) {
    int record_snapshot = WormSnapshotBefore(5);
    W[w_curr].x[p_idx] = x_new;
    W[w_curr].y[p_idx] = y_new;
    W[w_curr].z[p_idx] = z_new;
    worm.x_ira = x_new;
    worm.y_ira = y_new;
    worm.z_ira = z_new;
    worm_acc_iramove++;
    if(record_snapshot) WormSnapshotAfter(5);
  }
}

void PIMCWormMoves(void) {
  int sub, gap;
  DOUBLE r;

  if(worm_next == NULL) {
    InitWorm();
  }

  /* Realizar varios intentos por iteración para una buena exploración ergódica */
  for(sub = 0; sub < 5; sub++) {
    if(worm_open == 0) {
      worm_count_Z++;
      WormOpenMove();
    } else {
      worm_count_G++;

      gap = (worm_bead_masha - worm_bead + Nwalkers) % Nwalkers;
      if(gap < NwalkersMax && worm_gap_histogram != NULL) {
        worm_gap_histogram[gap]++;
      }
      if(gap > worm_max_gap_observed) worm_max_gap_observed = gap;

      r = Random();
      if(N > 1 && worm_swap_enabled) {
        if(gap == 1) {
          if(r < WORM_PROB_CLOSE_N_GT1) {
            WormCloseMove();
          } else if(r < WORM_PROB_CLOSE_N_GT1 + WORM_PROB_RECEDE_N_GAP1) {
            WormRecedeMove();
          } else if(r < WORM_PROB_CLOSE_N_GT1 + WORM_PROB_RECEDE_N_GAP1
                    + WORM_PROB_SWAP_N_GAP1) {
            WormSwapMove();
          } else {
            WormIraSpatialMove();
          }
        } else {
          if(r < WORM_PROB_ADVANCE_N_GT1) {
            WormAdvanceMove();
          } else if(r < WORM_PROB_ADVANCE_N_GT1 + WORM_PROB_RECEDE_N_GT1) {
            WormRecedeMove();
          } else if(r < WORM_PROB_ADVANCE_N_GT1 + WORM_PROB_RECEDE_N_GT1
                    + WORM_PROB_SWAP_N_GT1) {
            WormSwapMove();
          } else {
            WormIraSpatialMove();
          }
        }
      } else {
        /* Para N == 1: Advance, Recede, Spatial, Close */
        if(gap == 1) {
          if(r < WORM_PROB_CLOSE_GAP1) {
            WormCloseMove();
          } else if(r < WORM_PROB_CLOSE_GAP1 + WORM_PROB_RECEDE_GAP1) {
            WormRecedeMove();
          } else {
            WormIraSpatialMove();
          }
        } else {
          if(r < WORM_PROB_ADVANCE_GAP_GT1) {
            WormAdvanceMove();
          } else if(r < WORM_PROB_ADVANCE_GAP_GT1 + WORM_PROB_RECEDE_GAP_GT1) {
            WormRecedeMove();
          } else {
            WormIraSpatialMove();
          }
        }
      }
    }
  }
}

void MeasureWormPermutations(void) {
  int p, curr_p, curr_w, L;
  int *visited;

  if(worm_open != 0 || worm_next == NULL) return;

  if(!EnsureWormRingScratch(N, 0)) return;
  visited = worm_ring_visited;
  memset(visited, 0, (size_t)N * sizeof(int));

  worm_cycle_samples++;

  for(p = 0; p < N; p++) {
    if(!visited[p]) {
      L = 0;
      curr_p = p;
      curr_w = 0;
      do {
        if(curr_w == 0) {
          visited[curr_p] = 1;
          L++;
        }
        curr_p = worm_next[curr_w][curr_p];
        curr_w = (curr_w + 1) % Nwalkers;
        if(L > N) break;
      } while(!(curr_w == 0 && curr_p == p));

      if(worm_cycle_histogram != NULL && L <= N) {
        worm_cycle_histogram[L] += L; /* ponderado por número de partículas en ciclo L */
      }
    }
  }
}

void SaveWormPermutations(const char *filename) {
  FILE *f;
  int L;
  DOUBLE total_particles;

  if(worm_cycle_histogram == NULL || worm_cycle_samples == 0) return;
  f = fopen(filename, "w");
  if(f == NULL) return;

  total_particles = (DOUBLE)worm_cycle_samples * (DOUBLE)N;

  fprintf(f, "# cycle_length particle_count fraction\n");
  for(L = 1; L <= N; L++) {
    DOUBLE frac = total_particles > 0. ? (DOUBLE)worm_cycle_histogram[L] / total_particles : 0.;
    fprintf(f, "%d %lld %.6f\n", L, worm_cycle_histogram[L], frac);
  }
  fclose(f);
}

void WormDiagnostic(void) {
  long long total = worm_count_Z + worm_count_G;
  if(total == 0) total = 1;
  double frac_G = (double)worm_count_G / (double)total;
  double rate_open = worm_try_open > 0 ? (double)worm_acc_open / (double)worm_try_open * 100.0 : 0.0;
  double rate_close = worm_try_close > 0 ? (double)worm_acc_close / (double)worm_try_close * 100.0 : 0.0;
  double rate_adv = worm_try_advance > 0 ? (double)worm_acc_advance / (double)worm_try_advance * 100.0 : 0.0;
  double rate_rec = worm_try_recede > 0 ? (double)worm_acc_recede / (double)worm_try_recede * 100.0 : 0.0;
  double rate_swap = worm_try_swap > 0 ? (double)worm_acc_swap / (double)worm_try_swap * 100.0 : 0.0;
  double rate_ira = worm_try_iramove > 0 ? (double)worm_acc_iramove / (double)worm_try_iramove * 100.0 : 0.0;

  double p1 = 0., p2 = 0.;
  if(worm_cycle_samples > 0 && worm_cycle_histogram != NULL) {
    double tot = (double)worm_cycle_samples * (double)N;
    if(tot > 0.) {
      p1 = (double)worm_cycle_histogram[1] / tot;
      if(N >= 2) p2 = (double)worm_cycle_histogram[2] / tot;
    }
  }

  double rhos = 0.0;
  if(worm_sf_samples > 0) {
    rhos = (double)ComputeWormSuperfluidFraction(worm_sf_samples, worm_sf_A_sum, worm_sf_A2_sum, worm_sf_Icl_sum, worm_sf_W2_sum);
  }

  Message("Worm F4: Z=%lld G=%lld (fG=%.2f) | Op: %.1f%% | Cl: %.1f%% | Adv: %.1f%% | Rec: %.1f%% | Sw: %.1f%% | Ira: %.1f%% | P(1)=%.2f P(2)=%.2f | rhos/rho=%.4f\n",
          worm_count_Z, worm_count_G, frac_G,
          rate_open, rate_close, rate_adv, rate_rec, rate_swap, rate_ira,
          p1, p2, rhos);
}

void SaveWormGapHistogram(const char *filename) {
  FILE *f;
  int g;
  if(worm_gap_histogram == NULL) return;
  f = fopen(filename, "w");
  if(f == NULL) return;
  fprintf(f, "# gap count\n");
  for(g = 1; g <= WORM_MAX_GAP && g <= NwalkersMax; g++) {
    fprintf(f, "%d %lld\n", g, worm_gap_histogram[g]);
  }
  fclose(f);
}

void MeasureWormOBDM(void) {
  int gap, ix1, ix2, is;
  DOUBLE x1 = 0., y1 = 0., z1 = 0., x2 = 0., y2 = 0., z2 = 0.;
  DOUBLE dx, dy, dz, s2;

  if(worm_open == 0 || worm_obdm_matrix == NULL) return;

  SyncWormOBDMGeometry();

  gap = worm_bead_masha - worm_bead;
  if(gap < 0) gap += Nwalkers;

  CaseX(x1 = W[worm_bead_masha].x[worm_particle_masha]; x2 = W[worm_bead].x[worm_particle]);
  CaseY(y1 = W[worm_bead_masha].y[worm_particle_masha]; y2 = W[worm_bead].y[worm_particle]);
  CaseZ(z1 = W[worm_bead_masha].z[worm_particle_masha]; z2 = W[worm_bead].z[worm_particle]);

#ifdef TRIAL_1D
#ifdef BC_1DPBC_X
  ix1 = (int)floor((x1 + worm_obdm_max) * worm_obdm_inv_step_x);
  ix2 = (int)floor((x2 + worm_obdm_max) * worm_obdm_inv_step_x);
#else
  ix1 = (int)floor((z1 + worm_obdm_max) * worm_obdm_inv_step_x);
  ix2 = (int)floor((z2 + worm_obdm_max) * worm_obdm_inv_step_x);
#endif
#else
  int iy1, iy2;
  ix1 = (int)floor((x1 + worm_obdm_max) * worm_obdm_inv_step_x);
  ix2 = (int)floor((x2 + worm_obdm_max) * worm_obdm_inv_step_x);
  iy1 = (int)floor((y1 + worm_obdm_max) * worm_obdm_inv_step_x);
  iy2 = (int)floor((y2 + worm_obdm_max) * worm_obdm_inv_step_x);
#ifdef TRIAL_3D
  int iz1 = (int)floor((z1 + worm_obdm_max) * worm_obdm_inv_step_x);
  int iz2 = (int)floor((z2 + worm_obdm_max) * worm_obdm_inv_step_x);
#endif
#endif

  dx = x1 - x2;
  dy = y1 - y2;
  dz = z1 - z2;
  s2 = FindNearestImage(&dx, &dy, &dz);
  is = (s2 < worm_obdm_smax2) ? (int)floor(Sqrt(s2) * worm_obdm_inv_step_s) : worm_obdm_grid;

  worm_obdm_samples_all++;
  if(ix1 >= 0 && ix1 < worm_obdm_grid && ix2 >= 0 && ix2 < worm_obdm_grid) {
    worm_obdm_matrix_all[ix1][ix2]++;
    worm_obdm_matrix_all[ix2][ix1]++;
  }
#ifndef TRIAL_1D
  if(iy1 >= 0 && iy1 < worm_obdm_grid && iy2 >= 0 && iy2 < worm_obdm_grid) {
    worm_obdm_matrix_all[iy1][iy2]++;
    worm_obdm_matrix_all[iy2][iy1]++;
  }
#ifdef TRIAL_3D
  if(iz1 >= 0 && iz1 < worm_obdm_grid && iz2 >= 0 && iz2 < worm_obdm_grid) {
    worm_obdm_matrix_all[iz1][iz2]++;
    worm_obdm_matrix_all[iz2][iz1]++;
  }
#endif
#endif
  if(is >= 0 && is < worm_obdm_grid) {
    worm_obdm_dist_all[is]++;
  }

  if(gap == 1) {
    int ir1, ir2;
#ifdef TRIAL_1D
#ifdef BC_1DPBC_X
    ir1 = (int)floor(fabs(x1) * worm_obdm_inv_step_r);
    ir2 = (int)floor(fabs(x2) * worm_obdm_inv_step_r);
#else
    ir1 = (int)floor(fabs(z1) * worm_obdm_inv_step_r);
    ir2 = (int)floor(fabs(z2) * worm_obdm_inv_step_r);
#endif
#elif defined(TRIAL_2D)
    DOUBLE r1_2 = x1 * x1 + y1 * y1;
    DOUBLE r2_2 = x2 * x2 + y2 * y2;
    ir1 = (r1_2 < worm_obdm_max2) ? (int)floor(Sqrt(r1_2) * worm_obdm_inv_step_r) : worm_obdm_grid;
    ir2 = (r2_2 < worm_obdm_max2) ? (int)floor(Sqrt(r2_2) * worm_obdm_inv_step_r) : worm_obdm_grid;
#else
    DOUBLE r1_2 = x1 * x1 + y1 * y1 + z1 * z1;
    DOUBLE r2_2 = x2 * x2 + y2 * y2 + z2 * z2;
    ir1 = (r1_2 < worm_obdm_max2) ? (int)floor(Sqrt(r1_2) * worm_obdm_inv_step_r) : worm_obdm_grid;
    ir2 = (r2_2 < worm_obdm_max2) ? (int)floor(Sqrt(r2_2) * worm_obdm_inv_step_r) : worm_obdm_grid;
#endif

    worm_obdm_samples_gap1++;
    if(ix1 >= 0 && ix1 < worm_obdm_grid && ix2 >= 0 && ix2 < worm_obdm_grid) {
      worm_obdm_matrix[ix1][ix2]++;
      worm_obdm_matrix[ix2][ix1]++;
    }
#ifndef TRIAL_1D
    if(iy1 >= 0 && iy1 < worm_obdm_grid && iy2 >= 0 && iy2 < worm_obdm_grid) {
      worm_obdm_matrix[iy1][iy2]++;
      worm_obdm_matrix[iy2][iy1]++;
    }
#ifdef TRIAL_3D
    if(iz1 >= 0 && iz1 < worm_obdm_grid && iz2 >= 0 && iz2 < worm_obdm_grid) {
      worm_obdm_matrix[iz1][iz2]++;
      worm_obdm_matrix[iz2][iz1]++;
    }
#endif
#endif
    if(ir1 >= 0 && ir1 < worm_obdm_grid && ir2 >= 0 && ir2 < worm_obdm_grid) {
      worm_obdm_radial[ir1][ir2]++;
      worm_obdm_radial[ir2][ir1]++;
    }
    if(is >= 0 && is < worm_obdm_grid) {
      worm_obdm_dist[is]++;
    }
  }
}

void SaveWormOBDM(const char *filename) {
  FILE *f;
  int i, j;
  DOUBLE step_x, step_r, x_val, xp_val, norm_trace, trace_sum;

  if(worm_obdm_matrix == NULL || worm_obdm_samples_gap1 == 0) return;

  step_x = (2.0 * worm_obdm_max) / (DOUBLE)worm_obdm_grid;
  step_r = worm_obdm_max / (DOUBLE)worm_obdm_grid;

  trace_sum = 0.0;
  for(i = 0; i < worm_obdm_grid; i++) {
    trace_sum += (DOUBLE)worm_obdm_matrix[i][i];
  }
  if(trace_sum <= 0.0) return;
  norm_trace = (DOUBLE)N / (trace_sum * step_x);

  f = fopen(filename, "w");
  if(f == NULL) return;

  fprintf(f, "# WORM_OBDM_MATRIX_GAP1\n");
    fprintf(f, "# grid: %d min: %.6f max: %.6f step: %.6f samples: %lld norm_trace: %.15e N: %d tau: %.15e\n",
      worm_obdm_grid, -worm_obdm_max, worm_obdm_max, step_x, worm_obdm_samples_gap1,
      norm_trace, N, 1.0 / (T * (DOUBLE)Nwalkers));
  fprintf(f, "# x x_prime count rho1_normalized\n");

  for(i = 0; i < worm_obdm_grid; i++) {
    x_val = -worm_obdm_max + ((DOUBLE)i + 0.5) * step_x;
    for(j = 0; j < worm_obdm_grid; j++) {
      xp_val = -worm_obdm_max + ((DOUBLE)j + 0.5) * step_x;
      fprintf(f, "%.6f %.6f %lld %.15e\n",
              x_val, xp_val, worm_obdm_matrix[i][j],
              norm_trace * (DOUBLE)worm_obdm_matrix[i][j]);
    }
    fprintf(f, "\n");
  }
  fclose(f);

  f = fopen("outworm_obdm_dist.dat", "w");
  if(f != NULL) {
    DOUBLE s_max = 2.0 * worm_obdm_max;
    DOUBLE step_s = s_max / (DOUBLE)worm_obdm_grid;
    long long dist_sum = 0;
    for(i = 0; i < worm_obdm_grid; i++) dist_sum += worm_obdm_dist[i];
    fprintf(f, "# WORM_OBDM_DISTANCE_GAP1: samples=%lld dist_sum=%lld\n", worm_obdm_samples_gap1, dist_sum);
    fprintf(f, "# s count p_s\n");
    for(i = 0; i < worm_obdm_grid; i++) {
      DOUBLE s_val = ((DOUBLE)i + 0.5) * step_s;
      DOUBLE p_s = (dist_sum > 0) ? ((DOUBLE)worm_obdm_dist[i] / ((DOUBLE)dist_sum * step_s)) : 0.0;
      fprintf(f, "%.6f %lld %.15e\n", s_val, worm_obdm_dist[i], p_s);
    }
    fclose(f);
  }

  f = fopen("outworm_obdm_endpoint_radii.dat", "w");
  if(f != NULL) {
    DOUBLE radial_sum = 0.0;
    for(i = 0; i < worm_obdm_grid; i++) {
      for(j = 0; j < worm_obdm_grid; j++) {
        radial_sum += (DOUBLE)worm_obdm_radial[i][j];
      }
    }
    fprintf(f, "# WORM_OBDM_ENDPOINT_RADII_GAP1: samples=%lld\n", worm_obdm_samples_gap1);
    fprintf(f, "# r_masha r_ira count p_radii\n");
    for(i = 0; i < worm_obdm_grid; i++) {
      DOUBLE r_masha = ((DOUBLE)i + 0.5) * step_r;
      for(j = 0; j < worm_obdm_grid; j++) {
        DOUBLE r_ira = ((DOUBLE)j + 0.5) * step_r;
        DOUBLE p_radii = radial_sum > 0.0
          ? (DOUBLE)worm_obdm_radial[i][j] / (radial_sum * step_r * step_r)
          : 0.0;
        fprintf(f, "%.6f %.6f %lld %.15e\n",
                r_masha, r_ira, worm_obdm_radial[i][j], p_radii);
      }
      fprintf(f, "\n");
    }
    fclose(f);
  }

  f = fopen("outworm_obdm_all.dat", "w");
  if(f != NULL) {
    fprintf(f, "# WORM_OBDM_MATRIX_ALL_GAPS\n");
    fprintf(f, "# grid: %d min: %.6f max: %.6f step: %.6f samples: %lld\n",
            worm_obdm_grid, -worm_obdm_max, worm_obdm_max, step_x, worm_obdm_samples_all);
    for(i = 0; i < worm_obdm_grid; i++) {
      x_val = -worm_obdm_max + ((DOUBLE)i + 0.5) * step_x;
      for(j = 0; j < worm_obdm_grid; j++) {
        xp_val = -worm_obdm_max + ((DOUBLE)j + 0.5) * step_x;
        fprintf(f, "%.6f %.6f %lld\n", x_val, xp_val, worm_obdm_matrix_all[i][j]);
      }
      fprintf(f, "\n");
    }
    fclose(f);
  }
}

DOUBLE ComputeWormSuperfluidFraction(long long samples, const DOUBLE A_sum[3], const DOUBLE A2_sum[3], const DOUBLE Icl_sum[3], DOUBLE W2_sum) {
  (void)A_sum;
  if(samples <= 0 || T <= 0.0 || N <= 0) return 0.0;
  DOUBLE n_samples = (DOUBLE)samples;

#ifdef BC_ABSENT
#ifdef TRIAL_2D
  DOUBLE A2_mean = A2_sum[2] / n_samples;
  DOUBLE Icl_mean = Icl_sum[2] / n_samples;
  if(Icl_mean <= 0.0) return 0.0;
  return (4.0 * T / Icl_mean) * A2_mean;
#elif defined(TRIAL_3D)
  int axis, valid_axes = 0;
  DOUBLE rhos_sum = 0.0;
  for(axis = 0; axis < 3; axis++) {
    DOUBLE A2_mean = A2_sum[axis] / n_samples;
    DOUBLE Icl_mean = Icl_sum[axis] / n_samples;
    if(Icl_mean > 0.0) {
      rhos_sum += (4.0 * T / Icl_mean) * A2_mean;
      valid_axes++;
    }
  }
  return (valid_axes > 0) ? (rhos_sum / (DOUBLE)valid_axes) : 0.0;
#else
  (void)A2_sum; (void)Icl_sum; (void)W2_sum;
  return 0.0;
#endif
#else
  /* Periodic boundary conditions: winding number estimator rho_s/rho = T * <W^2> / (D * N) */
  (void)A2_sum; (void)Icl_sum;
  DOUBLE W2_mean = W2_sum / n_samples;
  return (T * W2_mean) / ((DOUBLE)D * (DOUBLE)N);
#endif
}

void MeasureWormSuperfluid(void) {
  int w, p;
  (void)w;
  (void)p;

  if(worm_open != 0 || worm_next == NULL || Nwalkers <= 0 || N <= 0 || T <= 0.0) return;

#ifdef BC_ABSENT
#ifdef TRIAL_2D
  DOUBLE Az_twice = 0.0;
  DOUBLE Icl = 0.0;
  for(w = 0; w < Nwalkers; w++) {
    int w_next = (w + 1 < Nwalkers) ? (w + 1) : 0;
    for(p = 0; p < N; p++) {
      int p_next = worm_next[w][p];
      if(p_next < 0) return;
      DOUBLE x = W[w].x[p], y = W[w].y[p];
      Az_twice += x * W[w_next].y[p_next] - y * W[w_next].x[p_next];
      Icl += x * x + y * y;
    }
  }
  DOUBLE Az = 0.5 * Az_twice;
  Icl /= (DOUBLE)Nwalkers;

  worm_sf_samples++;
  worm_sf_A_sum[2] += Az;
  worm_sf_A2_sum[2] += Az * Az;
  worm_sf_Icl_sum[2] += Icl;

  worm_sf_block_samples++;
  worm_sf_block_A_sum[2] += Az;
  worm_sf_block_A2_sum[2] += Az * Az;
  worm_sf_block_Icl_sum[2] += Icl;
#elif defined(TRIAL_3D)
  DOUBLE A[3] = {0.0, 0.0, 0.0};
  DOUBLE Icl[3] = {0.0, 0.0, 0.0};
  int axis;
  for(w = 0; w < Nwalkers; w++) {
    int w_next = (w + 1 < Nwalkers) ? (w + 1) : 0;
    for(p = 0; p < N; p++) {
      int p_next = worm_next[w][p];
      if(p_next < 0) return;
      DOUBLE x1 = W[w].x[p], y1 = W[w].y[p], z1 = W[w].z[p];
      DOUBLE x2 = W[w_next].x[p_next], y2 = W[w_next].y[p_next], z2 = W[w_next].z[p_next];
      A[0] += y1 * z2 - z1 * y2;
      A[1] += z1 * x2 - x1 * z2;
      A[2] += x1 * y2 - y1 * x2;
      Icl[0] += y1 * y1 + z1 * z1;
      Icl[1] += z1 * z1 + x1 * x1;
      Icl[2] += x1 * x1 + y1 * y1;
    }
  }
  worm_sf_samples++;
  worm_sf_block_samples++;
  for(axis = 0; axis < 3; axis++) {
    A[axis] *= 0.5;
    Icl[axis] /= (DOUBLE)Nwalkers;
    worm_sf_A_sum[axis] += A[axis];
    worm_sf_A2_sum[axis] += A[axis] * A[axis];
    worm_sf_Icl_sum[axis] += Icl[axis];
    worm_sf_block_A_sum[axis] += A[axis];
    worm_sf_block_A2_sum[axis] += A[axis] * A[axis];
    worm_sf_block_Icl_sum[axis] += Icl[axis];
  }
#endif
#else
  DOUBLE Wx = 0.0, Wy = 0.0, Wz = 0.0;
  for(w = 0; w < Nwalkers; w++) {
    int w_next = (w + 1 < Nwalkers) ? (w + 1) : 0;
    for(p = 0; p < N; p++) {
      int p_next = worm_next[w][p];
      if(p_next < 0) return;
      DOUBLE dr[3] = {
        W[w_next].x[p_next] - W[w].x[p],
        W[w_next].y[p_next] - W[w].y[p],
        W[w_next].z[p_next] - W[w].z[p]
      };
      FindNearestImage(&dr[0], &dr[1], &dr[2]);
      CaseX(Wx += dr[0]);
      CaseY(Wy += dr[1]);
      CaseZ(Wz += dr[2]);
    }
  }
  DOUBLE W2 = Wx * Wx + Wy * Wy + Wz * Wz;
  worm_sf_samples++;
  worm_sf_W2_sum += W2;
  worm_sf_block_samples++;
  worm_sf_block_W2_sum += W2;
#endif
}

void AccumulateWormSuperfluidBlock(long long block_samples, DOUBLE rho_s_block) {
  DOUBLE w, w_new, delta;
  if(block_samples <= 0) return;

  w = (DOUBLE)block_samples;
  w_new = worm_sf_block_weight_sum + w;
  delta = rho_s_block - worm_sf_rho_block_mean;

  worm_sf_num_blocks++;
  worm_sf_rho_block_mean += (w / w_new) * delta;
  worm_sf_rho_block_m2 += w * delta * (rho_s_block - worm_sf_rho_block_mean);
  worm_sf_block_weight_sum = w_new;
  worm_sf_block_weight_sum2 += w * w;
}

void SaveWormSuperfluid(const char *filename) {
  FILE *f;
  DOUBLE rhos, block_mean = 0.0, block_err = 0.0;
  DOUBLE n_samples;

  if(worm_sf_samples == 0) return;
  f = fopen(filename, "w");
  if(f == NULL) return;

  n_samples = (DOUBLE)worm_sf_samples;
  rhos = ComputeWormSuperfluidFraction(worm_sf_samples, worm_sf_A_sum, worm_sf_A2_sum, worm_sf_Icl_sum, worm_sf_W2_sum);

  if(worm_sf_num_blocks > 0) {
    block_mean = worm_sf_rho_block_mean;
    if(worm_sf_num_blocks > 1 && worm_sf_block_weight_sum > 0.0) {
      DOUBLE denom = worm_sf_block_weight_sum * (worm_sf_block_weight_sum * worm_sf_block_weight_sum - worm_sf_block_weight_sum2);
      if(denom > 0.0 && worm_sf_rho_block_m2 > 0.0) {
        block_err = Sqrt((worm_sf_rho_block_m2 * worm_sf_block_weight_sum2) / denom);
      }
    }
  }

#ifdef BC_ABSENT
#ifdef TRIAL_2D
  fprintf(f, "# WORM_SUPERFLUID_AREA_ESTIMATOR_2D\n");
  fprintf(f, "# N=%d Nwalkers=%d T=%.15g beta=%.15g swap=%d samples=%lld blocks=%lld\n",
          N, Nwalkers, T, 1.0 / T, worm_swap_enabled, worm_sf_samples, worm_sf_num_blocks);
  fprintf(f, "# rho_s_over_rho rho_s_block_mean rho_s_block_err A2_mean A_mean Icl_mean\n");
  fprintf(f, "%.15e %.15e %.15e %.15e %.15e %.15e\n",
          rhos, block_mean, block_err,
          worm_sf_A2_sum[2] / n_samples,
          worm_sf_A_sum[2] / n_samples,
          worm_sf_Icl_sum[2] / n_samples);
#elif defined(TRIAL_3D)
  fprintf(f, "# WORM_SUPERFLUID_AREA_ESTIMATOR_3D\n");
  fprintf(f, "# N=%d Nwalkers=%d T=%.15g beta=%.15g swap=%d samples=%lld blocks=%lld\n",
          N, Nwalkers, T, 1.0 / T, worm_swap_enabled, worm_sf_samples, worm_sf_num_blocks);
  fprintf(f, "# rho_s_over_rho rho_s_block_mean rho_s_block_err Ax2_mean Ay2_mean Az2_mean Icl_x Icl_y Icl_z\n");
  fprintf(f, "%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",
          rhos, block_mean, block_err,
          worm_sf_A2_sum[0] / n_samples, worm_sf_A2_sum[1] / n_samples, worm_sf_A2_sum[2] / n_samples,
          worm_sf_Icl_sum[0] / n_samples, worm_sf_Icl_sum[1] / n_samples, worm_sf_Icl_sum[2] / n_samples);
#endif
#else
  fprintf(f, "# WORM_SUPERFLUID_WINDING_NUMBER\n");
  fprintf(f, "# N=%d Nwalkers=%d T=%.15g beta=%.15g swap=%d samples=%lld blocks=%lld\n",
          N, Nwalkers, T, 1.0 / T, worm_swap_enabled, worm_sf_samples, worm_sf_num_blocks);
  fprintf(f, "# rho_s_over_rho rho_s_block_mean rho_s_block_err W2_mean\n");
  fprintf(f, "%.15e %.15e %.15e %.15e\n",
          rhos, block_mean, block_err, worm_sf_W2_sum / n_samples);
#endif
  fclose(f);
}
#endif /* WORM_ALGORITHM */
