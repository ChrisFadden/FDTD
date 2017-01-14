/** @file updateH.cpp
 *  @brief Update equations for the magnetic field in free space
 *
 *  This implements the loops for updating the magnetic field
 *  in free space, not including materials or PML.
 *
 *  The basic update equation in free space is:
 *
 *  \f[
 *      Hx^{n+\frac{1}{2}}_{i,j+\frac{1}{2}} =
 * Hx^{n-\frac{1}{2}}_{i,j+\frac{1}{2}}
 *      - \frac{\Delta t}{\mu_0 \Delta x}(Ez^{n}_{i,j+1} - Ez^{n}_{i,j})
 *  \f]
 *
 *  \f[
 *      Hy^{n+\frac{1}{2}}_{i+\frac{1}{2},j} =
 * Hy^{n+\frac{1}{2}}_{i+\frac{1}{2},j}
 *      + \frac{\Delta t}{\mu_0 \Delta y}(Ez^{n}_{i+1,j} - Ez^{n}_{i,j})
 *  \f]
 */

#include "update.h"
#include <cstddef>
#include <iostream>

#define Ez_dy curlA(i, j, k)
#define Ey_dz curlB(i, j, k)
#define Ez_dx curlA(i, j, k)
#define Ex_dz curlB(i, j, k)
#define Ey_dx curlA(i, j, k)
#define Ex_dy curlB(i, j, k)

void updateH(field *f, grid *g, PML *p) {

  calcCurlHx(f, g);
  for (loop i = 0; i != nx - 1; ++i) {
    for (loop j = 0; j != ny - 1; ++j) {
      for (loop k = 1; k != nz - 1; ++k) {
        Hx(i, j, k) = Da * Hx(i, j, k) + Db * (Ez_dy - Ey_dz);
      }

      // Hx PML z-direction
      for (loop k = 1; k != pml_size; ++k) {
        psiHx_z1(i, j, k - 1) = p->bh.at(k - 1) * psiHx_z1(i, j, k - 1) -
                                p->ch.at(k - 1) * p->kappa_h.at(k - 1) * Ey_dz;
        Hx(i, j, k) = Hx(i, j, k) + Db * psiHx_z1(i, j, k - 1);
      }

      loop kk = pml_size - 2;
      for (loop k = nz - pml_size; k != nz - 1; ++k) {
        psiHx_z2(i, j, kk) = p->bh.at(kk) * psiHx_z2(i, j, kk) -
                             p->ch.at(kk) * p->kappa_h.at(kk) * Ey_dz;
        Hx(i, j, k) = Hx(i, j, k) + Db * psiHx_z2(i, j, kk);
        kk = kk - 1;
      }
    }

    // Hx PML y-direction
    for (loop k = 1; k != nz - 1; ++k) {
      for (loop j = 0; j != pml_size - 1; ++j) {
        psiHx_y1(i, j, k) = p->bh.at(j) * psiHx_y1(i, j, k) +
                            p->ch.at(j) * p->kappa_h.at(j) * Ez_dy;
        Hx(i, j, k) = Hx(i, j, k) + Db * psiHx_y1(i, j, k);
      }

      loop jj = pml_size - 2;
      for (loop j = ny - pml_size; j != ny - 1; ++j) {
        psiHx_y2(i, jj, k) = p->bh.at(jj) * psiHx_y2(i, jj, k) +
                             p->ch.at(jj) * p->kappa_h.at(jj) * Ez_dy;
        Hx(i, j, k) = Hx(i, j, k) + Db * psiHx_y2(i, jj, k);
        jj = jj - 1;
      }
    }
  }

  calcCurlHy(f, g);
  for (loop j = 0; j != ny - 1; ++j) {
    for (loop i = 0; i != nx - 1; ++i) {
      for (loop k = 1; k != nz - 1; ++k) {
        Hy(i, j, k) = Da * Hy(i, j, k) + Db * (Ez_dx - Ex_dz);
      }

      // Hy PML z-direction
      for (loop k = 1; k != pml_size; ++k) {
        psiHy_z1(i, j, k - 1) = p->bh.at(k - 1) * psiHy_z1(i, j, k - 1) -
                                p->ch.at(k - 1) * p->kappa_h.at(k - 1) * Ex_dz;
        Hy(i, j, k) = Hy(i, j, k) + Db * psiHy_z1(i, j, k - 1);
      }

      loop kk = pml_size - 2;
      for (loop k = nz - pml_size; k != nz - 1; ++k) {
        psiHy_z2(i, j, kk) = p->bh.at(kk) * psiHy_z2(i, j, kk) -
                             p->ch.at(kk) * p->kappa_h.at(kk) * Ex_dz;
        Hy(i, j, k) = Hy(i, j, k) + Db * psiHy_z2(i, j, kk);
        kk = kk - 1;
      }
    }

    // Hy PML x-direction
    for (loop k = 1; k != nz - 1; ++k) {
      for (loop i = 0; i != pml_size - 1; ++i) {
        psiHy_x1(i, j, k) = p->bh.at(i) * psiHy_x1(i, j, k) +
                            p->ch.at(i) * p->kappa_h.at(i) * Ez_dx;
        Hy(i, j, k) = Hy(i, j, k) + Db * psiHy_x1(i, j, k);
      }

      loop ii = pml_size - 2;
      for (loop i = nx - pml_size; i != nx - 1; ++i) {
        psiHy_x2(ii, j, k) = p->bh.at(ii) * psiHy_x2(ii, j, k) +
                             p->ch.at(ii) * p->kappa_h.at(ii) * Ez_dx;
        Hy(i, j, k) = Hy(i, j, k) + Db * psiHy_x2(ii, j, k);
        ii = ii - 1;
      }
    }
  }

  calcCurlHz(f, g);
  for (loop k = 0; k != nz - 1; ++k) {
    for (loop i = 0; i != nx - 1; ++i) {
      for (loop j = 0; j != ny - 1; ++j) {
        Hz(i, j, k) = Da * Hz(i, j, k) + Db * (Ey_dx - Ex_dy);
      }

      // Hz PML y-direction
      for (loop j = 0; j != pml_size - 1; ++j) {
        psiHz_y1(i, j, k) = p->bh.at(j) * psiHz_y1(i, j, k) -
                            p->ch.at(j) * p->kappa_h.at(j) * Ex_dy;
        Hz(i, j, k) = Hz(i, j, k) + Db * psiHz_y1(i, j, k);
      }

      loop jj = pml_size - 2;
      for (loop j = ny - pml_size; j != ny - 1; ++j) {
        psiHz_y2(i, jj, k) = p->bh.at(jj) * psiHz_y2(i, jj, k) -
                             p->ch.at(jj) * p->kappa_h.at(jj) * Ex_dy;
        Hz(i, j, k) = Hz(i, j, k) + Db * psiHz_y2(i, jj, k);
        jj = jj - 1;
      }
    }

    // Hz PML x-direction
    for (loop j = 0; j != ny - 1; ++j) {
      for (loop i = 0; i != pml_size - 1; ++i) {
        psiHz_x1(i, j, k) = p->bh.at(i) * psiHz_x1(i, j, k) +
                            p->ch.at(i) * p->kappa_h.at(i) * Ey_dx;
        Hz(i, j, k) = Hz(i, j, k) + Db * psiHz_x1(i, j, k);
      }

      loop ii = pml_size - 2;
      for (loop i = nx - pml_size; i != nx - 1; ++i) {
        psiHz_x2(ii, j, k) = p->bh.at(ii) * psiHz_x2(ii, j, k) +
                             p->ch.at(ii) * p->kappa_h.at(ii) * Ey_dx;
        Hz(i, j, k) = Hz(i, j, k) + Db * psiHz_x2(ii, j, k);
        ii = ii - 1;
      }
    }
  }

  return;
}
