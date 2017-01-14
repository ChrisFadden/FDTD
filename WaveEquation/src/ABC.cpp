#include <iostream>
#include "Field.h"
#include "Grid.h"
#include "ABC.h"

ABC::ABC(Grid* grid, Field* field) {
  g = grid;
  f = field;

  ABC_C1 = (g->cc * g->dt - g->dx) / (g->cc * g->dt + g->dx);
  ABC_C2 = 2 * g->dx / (g->cc * g->dt + g->dx);
  ABC_C3 =
      (g->cc * g->dt) * (g->cc * g->dt) / (2 * g->dx * (g->cc * g->dt + g->dx));
}

void ABC::ApplyBC() {
  /****************
   * i = 0
   ***************/

  for (int j = 1; j < g->SizeY - 1; j++) {
    f->UNew(0, j) = -1 * f->UOld(1, j) +
                    ABC_C1 * (f->UNew(1, j) + f->UOld(0, j)) +
                    ABC_C2 * (f->U(0, j) + f->U(1, j)) +
                    ABC_C3 * (f->U(0, j + 1) - 2 * f->U(0, j) + f->U(0, j - 1) +
                              f->U(1, j + 1) - 2 * f->U(1, j) + f->U(1, j - 1));
  }

  /****************
   * i = SizeX
   ***************/
  int i = g->SizeX - 1;

  for (int j = 1; j < g->SizeY - 1; j++) {
    f->UNew(i, j) =
        -1 * f->UOld(i - 1, j) + ABC_C1 * (f->UNew(i - 1, j) + f->UOld(i, j)) +
        ABC_C2 * (f->U(i, j) + f->U(i - 1, j)) +
        ABC_C3 * (f->U(i, j + 1) - 2 * f->U(i, j) + f->U(i, j - 1) +
                  f->U(i - 1, j + 1) - 2 * f->U(i - 1, j) + f->U(i - 1, j - 1));
  }

  /***************
   * j = 0
   **************/

  for (i = 1; i < g->SizeX - 1; i++) {
    f->UNew(i, 0) = -1*f->UOld(i, 1) + ABC_C1 * (f->UNew(i, 1) + f->UOld(i, 0)) +
                    ABC_C2 * (f->U(i, 0) + f->U(i, 1)) +
                    ABC_C3 * (f->U(i + 1, 0) - 2 * f->U(i, 0) + f->U(i - 1, 0) +
                              f->U(i + 1, 1) - 2 * f->U(i, 1) + f->U(i - 1, 1));
  }

  /**************
   * j = SizeY
   *************/
  int j = g->SizeX - 1;

  for (i = 1; i < g->SizeX - 1; i++) {
    f->UNew(i, j) =
        -1*f->UOld(i, j - 1) + ABC_C1 * (f->UNew(i, j - 1) + f->UOld(i, j)) +
        ABC_C2 * (f->U(i, j) + f->U(i, j-1)) +
        ABC_C3 * (f->U(i + 1, j) - 2 * f->U(i, j) + f->U(i - 1, j) +
                  f->U(i + 1, j - 1) - 2 * f->U(i, j - 1) + f->U(i - 1, j - 1));
  }

  return;
}









