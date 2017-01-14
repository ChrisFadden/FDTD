/** @file preprocess.cpp
 *  @brief Validate input parameters, and allocate arrays
 */

#include "pml.h"
#include "process.h"
#include "source.h"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>

/**
 * @brief validate input parameters
 *
 * This function validates user inputted
 * parameters from a file.  If no file
 * is specified, it is assumed to be
 * parameters.txt.  Quantities derived
 * from user input parameters are set
 * in this function as well.
 */
void init_sim_parameters(field *f, grid *g, Source *&s, std::string file) {

  if (file.empty()) {
    file = "parameters.txt";
  }
  std::ifstream inputFile(file.c_str());
  std::string line;
  std::vector<std::string> parameters(18);

  // Read contents of file
  for (int i = 0; std::getline(inputFile, line); i++) {
    parameters[i] = line;
  }
  inputFile.close();

  // User Inputs
  real freq = atof(parameters[0].c_str()) * GHz;
  sizeX = atoi(parameters[1].c_str());
  sizeY = atoi(parameters[2].c_str());
  sizeZ = atoi(parameters[3].c_str());
  maxTime = atoi(parameters[4].c_str());

  i_src = atoi(parameters[6].c_str()) + pml_size;
  j_src = atoi(parameters[7].c_str()) + pml_size;
  k_src = atoi(parameters[8].c_str()) + pml_size;

  tfsf_x0 = atoi(parameters[9].c_str()) + pml_size;
  tfsf_y0 = atoi(parameters[10].c_str()) + pml_size;
  tfsf_z0 = atoi(parameters[11].c_str()) + pml_size;
  tfsf_x1 = atoi(parameters[12].c_str()) + pml_size;
  tfsf_y1 = atoi(parameters[13].c_str()) + pml_size;
  tfsf_z1 = atoi(parameters[14].c_str()) + pml_size;

  theta = atof(parameters[15].c_str()) * pi / 180;
  phi = atof(parameters[16].c_str()) * pi / 180;
  TFSFpsi = atof(parameters[17].c_str()) * pi / 180;

  // Derived Grid Quantities
  dx = cc / (freq * NLambda);
  dy = dx;
  dz = dx;
  dt = CFL * dx / cc;
  Da = 1.0;
  Db = dt / mu0;
  Ca = 1.0;
  Cb = dt / eps0;

  // Set Sources
  SourceFunc *srcFunc = NULL;
  if (!parameters[5].compare("Harmonic")) {
    srcFunc = new HarmonicSource(freq, g);
  } else if (!parameters[5].compare("Gaussian")) {
    srcFunc = new GaussianSource(freq, g);
  }

  if (tfsf_x0 != pml_size) {
    s = new TFSF_Source(f, g, srcFunc);
  } else {
    s = new Source(f, g, srcFunc);
  }

  return;
}
/**
 * @brief initialize all arrays to correct size
 */
void array_initialize(field *f, grid *g, PML *p, Source *s) {
  f->hx.resize(nx * (ny - 1) * (nz - 1));
  f->hy.resize((nx - 1) * ny * (nz - 1));
  f->hz.resize((nx - 1) * (ny - 1) * nz);

  f->ex.resize((nx - 1) * ny * nz);
  f->ey.resize(nx * (ny - 1) * nz);
  f->ez.resize(nx * ny * (nz - 1));

  f->DenHx.resize(((nx)-1));
  f->DenHy.resize(((ny)-1));
  f->DenHz.resize(((nz)-1));

  f->DenEx.resize(((nx)-1));
  f->DenEy.resize(((ny)-1));
  f->DenEz.resize(((nz)-1));

  f->CurlA.resize((nx) * (ny) * (nz));
  f->CurlB.resize((nx) * (ny) * (nz));

  p->initialize();
  s->initialize();
  return;
}
