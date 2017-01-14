/** @file PML_initialize.cpp
 *  @brief Initializes PML parameters
 */

#include "constants.h"
#include "pml.h"
#include <cstddef>
#include <iostream>
#include <math.h>
#define MAGIC_SCALE 0.6

typedef std::size_t loop;

/**
 * @brief Initializes PML parameters
 *
 * There are several variables relevant to the PML.  The PML acts as a complex
 * frequency shifted tensor post-multiplying  the permeability and permittivity.
 * The tensor is diagonal, and is governed by:
 * \f[
 *    S_{xx} = s_y s_z s_x^{-1}
 * \f]
 * \f[
 *    S_{yy} = s_x s_z s_y^{-1}
 * \f]
 * \f[
 *    S_{zz} = s_x s_y s_z^{-1}
 * \f]
 * Where \f$s_{w} = \kappa_{w} + \frac{\sigma_w}{\alpha_w + j \omega
 * \epsilon}\f$.  R(0) is the desired reflection error, m is a
 * constant for the PML, usually around 3. Cpsi are the coefficients
 * for the PML update equations.
 *
 * \f[
 *    \sigma_{max} = \frac{-(m+1) ln(R(0))}{2 \eta d_{pml}} =
 * \frac{0.8(m+1)}{\eta \Delta}
 * \f]
 * \f[
 *    \alpha(i) = \alpha_{max} (\frac{d_{pml} - i}{d_{pml}})^m
 * \f]
 *
 * \f[
 *    \kappa(i) = 1 + (\kappa_{max} - 1) (\frac{i}{d_{pml}})^m
 * \f]
 *
 * \f[
 *    C_{\psi}^{w,w}(i) = e^{-(\frac{\sigma(i)}{\kappa(i)} +
 * \alpha(i))(\frac{\Delta
 * t}{\epsilon_0})}
 * \f]
 *
 * \f[
 *    C_{\psi}^{w,v}(i) = \sigma(i) \frac{(C_{\psi}^{w,w} - 1)}{(\sigma(i) +
 * \kappa(i) \alpha(i)) / \kappa(i)}
 * \f]
 */
void PML::initialize() {

  // PML auxiliary arrays
  PsiHx_y1.resize((nx) * (pml_size - 1) * (nz));
  PsiHx_y2.resize((nx) * (pml_size - 1) * (nz));
  PsiHx_z1.resize((nx) * (ny - 1) * (pml_size - 1));
  PsiHx_z2.resize((nx) * (ny - 1) * (pml_size - 1));

  PsiHy_x1.resize((pml_size - 1) * (ny) * (nz));
  PsiHy_x2.resize((pml_size - 1) * (ny) * (nz));
  PsiHy_z1.resize((nx - 1) * (ny) * (pml_size - 1));
  PsiHy_z2.resize((nx - 1) * (ny) * (pml_size - 1));

  PsiHz_x1.resize((pml_size - 1) * (ny - 1) * (nz - 1));
  PsiHz_x2.resize((pml_size - 1) * (ny - 1) * (nz - 1));
  PsiHz_y1.resize((nx - 1) * (pml_size - 1) * (nz - 1));
  PsiHz_y2.resize((nx - 1) * (pml_size - 1) * (nz - 1));

  PsiEx_y1.resize((nx - 1) * (pml_size) * (nz - 1));
  PsiEx_y2.resize((nx - 1) * (pml_size) * (nz - 1));
  PsiEx_z1.resize((nx - 1) * (ny) * (pml_size));
  PsiEx_z2.resize((nx - 1) * (ny) * (pml_size));

  PsiEy_x1.resize((pml_size) * (ny - 1) * (nz - 1));
  PsiEy_x2.resize((pml_size) * (ny - 1) * (nz - 1));
  PsiEy_z1.resize((nx - 1) * (ny - 1) * (pml_size));
  PsiEy_z2.resize((nx - 1) * (ny - 1) * (pml_size));

  PsiEz_x1.resize((pml_size) * (ny) * (nz));
  PsiEz_x2.resize((pml_size) * (ny) * (nz));
  PsiEz_y1.resize((nx) * (pml_size) * (nz));
  PsiEz_y2.resize((nx) * (pml_size) * (nz));

  be.resize(pml_size);
  ce.resize(pml_size);
  bh.resize((pml_size)-1);
  ch.resize((pml_size)-1);

  kappa_e.resize(pml_size);
  kappa_h.resize((pml_size)-1);
  // Temporary arrays
  array alpha_e(pml_size);
  array sigma_e(pml_size);
  // array kappa_e(pml_size);

  array alpha_h((pml_size)-1);
  array sigma_h((pml_size)-1);
  // array kappa_h((pml_size)-1);

  real sigma_max = 0.75 * (0.8 * ((pml_m) + 1) / (dx * eta0));

  // E pml
  for (loop i = 0; i != pml_size; ++i) {
    sigma_e.at(i) =
        sigma_max * pow((pml_size - 1.0 - i) / (pml_size - 1), pml_m);
    alpha_e.at(i) = alpha_max * pow((i / (pml_size - 1)), pml_ma);
    kappa_e.at(i) =
        1.0 +
        (kappa_max - 1) * pow((pml_size - 1.0 - i) / (pml_size - 1), pml_m);

    be.at(i) =
        exp(-(sigma_e.at(i) / kappa_e.at(i) + alpha_e.at(i)) * dt / eps0);
    ce.at(i) = sigma_e.at(i) * (be.at(i) - 1) /
               (sigma_e.at(i) + kappa_e.at(i) * alpha_e.at(i)) / kappa_e.at(i);
  }

  // H pml
  for (loop i = 0; i != pml_size - 1; ++i) {
    sigma_h.at(i) =
        sigma_max *
        pow(((float)pml_size - 1 - i - 0.5) / (pml_size - 1), pml_m);
    alpha_h.at(i) = alpha_max * pow((i + 1 - 0.5) / (pml_size - 1), pml_ma);
    kappa_h.at(i) = 1.0 +
                    (kappa_max - 1) *
                        (pow((pml_size - 1 - i - 0.5) / (pml_size - 1), pml_m));

    bh.at(i) =
        exp(-(sigma_h.at(i) / kappa_h.at(i) + alpha_h.at(i)) * dt / eps0);
    ch.at(i) = sigma_h.at(i) * (bh.at(i) - 1) /
               (sigma_h.at(i) + kappa_h.at(i) * alpha_h.at(i)) / kappa_h.at(i);
  }

  // Denominator Coefficients
  loop ii = pml_size - 2;
  for (loop i = 0; i < nx - 1; ++i) {

    if (i < pml_size - 1) {
      denHx(i) = 1.0 / (kappa_h.at(i) * dx);
    } else if (i >= nx - pml_size) {
      denHx(i) = 1.0 / (kappa_h.at(ii) * dx);
      ii = ii - 1;
    } else {
      denHx(i) = 1.0 / dx;
    }
  }

  loop jj = pml_size - 2;
  for (loop j = 0; j < ny - 1; ++j) {

    if (j < pml_size - 1) {
      denHy(j) = 1.0 / (kappa_h.at(j) * dy);
    } else if (j >= ny - pml_size) {
      denHy(j) = 1.0 / (kappa_h.at(jj) * dy);
      jj = jj - 1;
    } else {
      denHy(j) = 1.0 / dy;
    }
  }

  loop kk = pml_size - 2;
  for (loop k = 1; k < nz - 1; ++k) {

    if (k < pml_size) {
      denHz(k) = 1.0 / (kappa_h.at(k - 1) * dz);
    } else if (k >= nz - pml_size) {
      denHz(k) = 1.0 / (kappa_h.at(kk) * dz);
      kk = kk - 1;
    } else {
      denHz(k) = 1.0 / dz;
    }
  }

  ii = pml_size - 1;
  for (loop i = 0; i < nx - 1; ++i) {

    if (i < pml_size) {
      denEx(i) = 1.0 / (kappa_e.at(i) * dx);
    } else if (i >= nx - pml_size) {
      denEx(i) = 1.0 / (kappa_e.at(ii) * dx);
      ii = ii - 1;
    } else {
      denEx(i) = 1.0 / dx;
    }
  }

  jj = pml_size - 1;
  for (loop j = 0; j < ny - 1; ++j) {

    if (j < pml_size) {
      denEy(j) = 1.0 / (kappa_e.at(j) * dy);
    } else if (j >= ny - pml_size) {
      denEy(j) = 1.0 / (kappa_e.at(jj) * dy);
      jj = jj - 1;
    } else {
      denEy(j) = 1.0 / dy;
    }
  }

  kk = pml_size - 1;
  for (loop k = 0; k < nz - 1; ++k) {

    if (k < pml_size) {
      denEz(k) = 1.0 / (kappa_e.at(k) * dz);
    } else if (k >= nz - 1 - pml_size) {
      denEz(k) = 1.0 / (kappa_e.at(kk) * dz);
      kk = kk - 1;
    } else {
      denEz(k) = 1.0 / dz;
    }
  }

  return;
}
