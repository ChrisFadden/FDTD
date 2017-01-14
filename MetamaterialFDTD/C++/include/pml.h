#ifndef PML_H
#define PML_H

#include "constants.h"
#include "field.h"
#include "grid.h"
class PML {
public:
  PML(field *Field, grid *Grid) : f(Field), g(Grid) {}
  void initialize();
  void updateH();
  void updateE();

  array PsiHx_y1, PsiHx_y2;
  array PsiHx_z1, PsiHx_z2;

  array PsiHy_x1, PsiHy_x2;
  array PsiHy_z1, PsiHy_z2;

  array PsiHz_x1, PsiHz_x2;
  array PsiHz_y1, PsiHz_y2;

  array PsiEx_y1, PsiEx_y2;
  array PsiEx_z1, PsiEx_z2;

  array PsiEy_x1, PsiEy_x2;
  array PsiEy_z1, PsiEy_z2;

  array PsiEz_x1, PsiEz_x2;
  array PsiEz_y1, PsiEz_y2;

  array bh, ch;
  array be, ce;

  array kappa_e, kappa_h;

  field *f;
  grid *g;
};

// array access macros
#define psiHx_y1(i, j, k)                                                      \
  p->PsiHx_y1.at(((i) * (pml_size - 1) + (j)) * (nz) + (k))
#define psiHx_y2(i, j, k)                                                      \
  p->PsiHx_y2.at(((i) * (pml_size - 1) + (j)) * (nz) + (k))
#define psiHx_z1(i, j, k)                                                      \
  p->PsiHx_z1.at(((i) * (ny - 1) + (j)) * (pml_size - 1) + (k))
#define psiHx_z2(i, j, k)                                                      \
  p->PsiHx_z2.at(((i) * (ny - 1) + (j)) * (pml_size - 1) + (k))

#define psiHy_x1(i, j, k) p->PsiHy_x1.at(((i)*ny + (j)) * (nz) + (k))
#define psiHy_x2(i, j, k) p->PsiHy_x2.at(((i)*ny + (j)) * (nz) + (k))
#define psiHy_z1(i, j, k) p->PsiHy_z1.at(((i)*ny + (j)) * (pml_size - 1) + (k))
#define psiHy_z2(i, j, k) p->PsiHy_z2.at(((i)*ny + (j)) * (pml_size - 1) + (k))

#define psiHz_x1(i, j, k)                                                      \
  p->PsiHz_x1.at(((i) * (ny - 1) + (j)) * (nz - 1) + (k))
#define psiHz_x2(i, j, k)                                                      \
  p->PsiHz_x2.at(((i) * (ny - 1) + (j)) * (nz - 1) + (k))
#define psiHz_y1(i, j, k)                                                      \
  p->PsiHz_y1.at(((i) * (pml_size - 1) + (j)) * ((nz)-1) + (k))
#define psiHz_y2(i, j, k)                                                      \
  p->PsiHz_y2.at(((i) * (pml_size - 1) + (j)) * ((nz)-1) + (k))

#define psiEx_y1(i, j, k) p->PsiEx_y1.at(((i)*pml_size + (j)) * (nz - 1) + (k))
#define psiEx_y2(i, j, k) p->PsiEx_y2.at(((i)*pml_size + (j)) * (nz - 1) + (k))
#define psiEx_z1(i, j, k) p->PsiEx_z1.at(((i) * (ny) + (j)) * (pml_size) + (k))
#define psiEx_z2(i, j, k) p->PsiEx_z2.at(((i) * (ny) + (j)) * (pml_size) + (k))

#define psiEy_x1(i, j, k)                                                      \
  p->PsiEy_x1.at(((i) * (ny - 1) + (j)) * (nz - 1) + (k))
#define psiEy_x2(i, j, k)                                                      \
  p->PsiEy_x2.at(((i) * (ny - 1) + (j)) * (nz - 1) + (k))
#define psiEy_z1(i, j, k)                                                      \
  p->PsiEy_z1.at(((i) * (ny - 1) + (j)) * (pml_size) + (k))
#define psiEy_z2(i, j, k)                                                      \
  p->PsiEy_z2.at(((i) * (ny - 1) + (j)) * (pml_size) + (k))

#define psiEz_x1(i, j, k) p->PsiEz_x1.at(((i)*ny + (j)) * (nz) + (k))
#define psiEz_x2(i, j, k) p->PsiEz_x2.at(((i)*ny + (j)) * (nz) + (k))
#define psiEz_y1(i, j, k) p->PsiEz_y1.at(((i)*pml_size + (j)) * (nz) + (k))
#define psiEz_y2(i, j, k) p->PsiEz_y2.at(((i)*pml_size + (j)) * (nz) + (k))

#endif
