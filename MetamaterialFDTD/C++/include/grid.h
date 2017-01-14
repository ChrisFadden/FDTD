#ifndef GRID_H
#define GRID_H

#include "constants.h"
#include <cstddef>
struct grid {
  real Dx, Dy, Dz;
  real Dt;
  real Theta, Phi, Psi;
  std::size_t SizeX, SizeY, SizeZ;
  std::size_t I_src, J_src, K_src;
  std::size_t TFSF_x0, TFSF_y0, TFSF_z0, TFSF_x1, TFSF_y1, TFSF_z1;
  std::size_t MaxTime;
};

// Variable Access Macros
#define dx g->Dx
#define dy g->Dy
#define dz g->Dz
#define dt g->Dt
#define theta g->Theta
#define phi g->Phi
#define TFSFpsi g->Psi
#define sizeX g->SizeX
#define sizeY g->SizeY
#define sizeZ g->SizeZ
#define i_src g->I_src
#define j_src g->J_src
#define k_src g->K_src
#define tfsf_x0 g->TFSF_x0
#define tfsf_x1 g->TFSF_x1
#define tfsf_y0 g->TFSF_y0
#define tfsf_y1 g->TFSF_y1
#define tfsf_z0 g->TFSF_z0
#define tfsf_z1 g->TFSF_z1
#define maxTime g->MaxTime

#define nx (sizeX + 2 * (pml_size))
#define ny (sizeY + 2 * (pml_size))
#define nz (sizeZ + 2 * (pml_size))

#endif
