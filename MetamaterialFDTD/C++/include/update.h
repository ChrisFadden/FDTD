#ifndef UPDATE_H
#define UPDATE_H
#include "field.h"
#include "grid.h"
#include "pml.h"
void updateE(field *f, grid *g, PML *p);
void calcCurlEx(field *f, grid *g);
void calcCurlEy(field *f, grid *g);
void calcCurlEz(field *f, grid *g);
void calcCurlHx(field *f, grid *g);
void calcCurlHy(field *f, grid *g);
void calcCurlHz(field *f, grid *g);
void updateH(field *f, grid *g, PML *p);
#endif
