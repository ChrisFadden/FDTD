#ifndef FIELD_H
#define FIELD_H

#include "constants.h"
#include "grid.h"

struct field {
  array hx;
  array hy;
  array hz;

  array ex;
  array ey;
  array ez;

  array DenHx;
  array DenHy;
  array DenHz;

  array DenEx;
  array DenEy;
  array DenEz;

  array CurlA;
  array CurlB;

  real DA;
  real DB;
  real CA;
  real CB;
};

// Access Arrays
#define Hx(i, j, k) f->hx.at(((i) * ((ny)-1) + (j)) * ((nz)-1) + (k))
#define Hy(i, j, k) f->hy.at(((i) * (ny) + (j)) * ((nz)-1) + (k))
#define Hz(i, j, k) f->hz.at(((i) * ((ny)-1) + (j)) * (nz) + (k))

#define Ex(i, j, k) f->ex.at(((i) * (ny) + (j)) * (nz) + (k))
#define Ey(i, j, k) f->ey.at(((i) * ((ny)-1) + (j)) * (nz) + (k))
#define Ez(i, j, k) f->ez.at(((i) * (ny) + (j)) * ((nz - 1)) + (k))

#define denHx(i) f->DenHx.at(i)
#define denHy(j) f->DenHy.at(j)
#define denHz(k) f->DenHz.at(k)

#define denEx(i) f->DenEx.at(i)
#define denEy(j) f->DenEy.at(j)
#define denEz(k) f->DenEz.at(k)

#define curlA(i, j, k) f->CurlA.at(((i) * (ny) + (j)) * (nz) + (k))
#define curlB(i, j, k) f->CurlB.at(((i) * (ny) + (j)) * (nz) + (k))

// Access variables
#define Da (f->DA)
#define Db (f->DB)
#define Ca (f->CA)
#define Cb (f->CB)

#endif
