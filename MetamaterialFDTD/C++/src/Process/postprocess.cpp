/** @file postprocess.cpp
 *  @brief Post processing of simulation results
 */

#include "process.h"
#include <fstream>
#include <iostream>
/**
 * @brief prints field to file
 *
 * This function prints the Ez
 * field to a given file, or output.csv
 * if no filename is given
 */

void printField(const field *f, const grid *g, std::string fname) {

  if (fname.empty()) {
    fname = "output.csv";
  }
  std::ofstream outFile(fname.c_str());

  /*
  for (loop j = pml_size; j != pml_size + sizeY; ++j) {
    for (loop i = pml_size; i != pml_size + sizeX; ++i) {
      outFile << Ez(i, j, k_src + 3) << ", ";
    }
    outFile << std::endl;
  }*/

  for (loop j = 0; j != ny; ++j) {
    for (loop i = 0; i != nx; ++i) {
      outFile << Ez(i, j, k_src + 3) << ", ";
    }
    outFile << std::endl;
  }

  outFile.close();
  return;
}
