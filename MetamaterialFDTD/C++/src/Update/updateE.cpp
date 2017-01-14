/** @file updateE.cpp
 *  @brief Update equations for the electric field in free space
 *
 *  This has the loops for updating the electric field
 *  in free space.  This does not include updating the
 *  field in PML regions, or any materials.
 *
 *  The basic update for the electric field in free space
 *  is:
 *  \f[
 *    Ez^{n+1}_{i,j} =  Ez^{n-1}_{i,j} + \frac{\Delta t}{\epsilon_0 \Delta x}
 *    (Hy^{n+\frac{1}{2}}_{i+\frac{1}{2},j} -
 * Hy^{n+\frac{1}{2}}_{i-\frac{1}{2},j})
 *    - \frac{\Delta t}{\epsilon_0 \Delta y}
 *    (Hx^{n+\frac{1}{2}}_{i,j+\frac{1}{2}} -
 * Hx^{n+\frac{1}{2}}_{i,j-\frac{1}{2}})
 *  \f]
 */

#include "update.h"
#include <cstddef>
#include <iostream>

#define Hz_dy curlA(i, j, k)
#define Hy_dz curlB(i, j, k)
#define Hz_dx curlA(i, j, k)
#define Hx_dz curlB(i, j, k)
#define Hy_dx curlA(i, j, k)
#define Hx_dy curlB(i, j, k)

void updateE(field *f, grid *g, PML *p) {

  calcCurlEx(f, g);
  for (loop i = 0; i != nx - 1; ++i) {
    for (loop j = 1; j != ny - 1; ++j) {
      for (loop k = 0; k != nz - 1; ++k) {
        Ex(i, j, k) = Ca * Ex(i, j, k) + Cb * (Hz_dy - Hy_dz);
      }

      // Ex PML z-direction
      for (loop k = 0; k != pml_size; ++k) {
        psiEx_z1(i, j, k) = p->be.at(k) * psiEx_z1(i, j, k) -
                            p->ce.at(k) * p->kappa_e.at(k) * Hy_dz;
        Ex(i, j, k) = Ex(i, j, k) + Cb * psiEx_z1(i, j, k);
      }

      loop kk = pml_size - 1;
      for (loop k = nz - pml_size - 1; k != nz - 1; ++k) {
        psiEx_z2(i, j, kk) = p->be.at(kk) * psiEx_z2(i, j, kk) -
                             p->ce.at(kk) * p->kappa_e.at(kk) * Hy_dz;
        Ex(i, j, k) = Ex(i, j, k) + Cb * psiEx_z2(i, j, kk);
        kk = kk - 1;
      }
    }

    // Ex PML y-direction
    for (loop k = 0; k != nz - 1; ++k) {
      for (loop j = 1; j != pml_size; ++j) {
        psiEx_y1(i, j, k) = p->be.at(j) * psiEx_y1(i, j, k) +
                            p->ce.at(j) * p->kappa_e.at(j) * Hz_dy;
        Ex(i, j, k) = Ex(i, j, k) + Cb * psiEx_y1(i, j, k);
      }

      loop jj = pml_size - 1;
      for (loop j = ny - pml_size; j != ny - 1; ++j) {
        psiEx_y2(i, jj, k) = p->be.at(jj) * psiEx_y2(i, jj, k) +
                             p->ce.at(jj) * p->kappa_e.at(jj) * Hz_dy;
        Ex(i, j, k) = Ex(i, j, k) + Cb * psiEx_y2(i, jj, k);
        jj = jj - 1;
      }
    }
  }

  calcCurlEy(f, g);
  for (loop j = 0; j != ny - 1; ++j) {
    for (loop i = 1; i != nx - 1; ++i) {
      for (loop k = 0; k != nz - 1; ++k) {
        Ey(i, j, k) = Ca * Ey(i, j, k) + Cb * (Hz_dx - Hx_dz);
      }

      // Ey PML z-direction
      for (loop k = 0; k != pml_size; ++k) {
        psiEy_z1(i, j, k) = p->be.at(k) * psiEy_z1(i, j, k) -
                            p->ce.at(k) * p->kappa_e.at(k) * Hx_dz;
        Ey(i, j, k) = Ey(i, j, k) + Cb * psiEy_z1(i, j, k);
      }

      loop kk = pml_size - 1;
      for (loop k = nz - pml_size - 1; k != nz - 1; ++k) {
        psiEy_z2(i, j, kk) = p->be.at(kk) * psiEy_z2(i, j, kk) -
                             p->ce.at(kk) * p->kappa_e.at(kk) * Hx_dz;
        Ey(i, j, k) = Ey(i, j, k) + Cb * psiEy_z2(i, j, kk);
        kk = kk - 1;
      }
    }

    // Ey PML x-direction
    for (loop k = 0; k != nz - 1; ++k) {
      for (loop i = 1; i != pml_size; ++i) {
        psiEy_x1(i, j, k) = p->be.at(i) * psiEy_x1(i, j, k) +
                            p->ce.at(i) * p->kappa_e.at(i) * Hz_dx;
        Ey(i, j, k) = Ey(i, j, k) + Cb * psiEy_x1(i, j, k);
      }

      loop ii = pml_size - 1;
      for (loop i = nx - pml_size; i != nx - 1; ++i) {
        psiEy_x2(ii, j, k) = p->be.at(ii) * psiEy_x2(ii, j, k) +
                             p->ce.at(ii) * p->kappa_e.at(ii) * Hz_dx;
        Ey(i, j, k) = Ey(i, j, k) + Cb * psiEy_x2(ii, j, k);
        ii = ii - 1;
      }
    }
  }

  calcCurlEz(f, g);
  for (loop k = 1; k != nz - 1; ++k) {
    for (loop i = 1; i != nx - 1; ++i) {
      for (loop j = 1; j != ny - 1; ++j) {
        Ez(i, j, k) = Ca * Ez(i, j, k) + Cb * (Hy_dx - Hx_dy);
      }

      // Ez PML y-direction
      for (loop j = 1; j != pml_size; ++j) {
        psiEz_y1(i, j, k) = p->be.at(j) * psiEz_y1(i, j, k) -
                            p->ce.at(j) * p->kappa_e.at(j) * Hx_dy;
        Ez(i, j, k) = Ez(i, j, k) + Db * psiEz_y1(i, j, k);
      }

      loop jj = pml_size - 1;
      for (loop j = ny - pml_size; j != ny - 1; ++j) {
        psiEz_y2(i, jj, k) = p->be.at(jj) * psiEz_y2(i, jj, k) -
                             p->ce.at(jj) * p->kappa_e.at(jj) * Hx_dy;
        Ez(i, j, k) = Ez(i, j, k) + Db * psiEz_y2(i, jj, k);
        jj = jj - 1;
      }
    }

    // Ez PML x-direction
    for (loop j = 1; j != ny - 1; ++j) {
      for (loop i = 1; i != pml_size; ++i) {
        psiEz_x1(i, j, k) = p->be.at(i) * psiEz_x1(i, j, k) +
                            p->ce.at(i) * p->kappa_e.at(i) * Hy_dx;
        Ez(i, j, k) = Ez(i, j, k) + Cb * psiEz_x1(i, j, k);
      }

      loop ii = pml_size - 1;
      for (loop i = nx - pml_size; i != nx - 1; ++i) {
        psiEz_x2(ii, j, k) = p->be.at(ii) * psiEz_x2(ii, j, k) +
                             p->ce.at(ii) * p->kappa_e.at(ii) * Hy_dx;
        Ez(i, j, k) = Ez(i, j, k) + Db * psiEz_x2(ii, j, k);
        ii = ii - 1;
      }
    }
  }
  return;
}
