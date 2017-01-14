#include "constants.h"
#include "error.h"
#include "field.h"
#include "grid.h"
#include "pml.h"
#include "process.h"
#include "source.h"
#include "update.h"
#include <cstddef>
#include <fstream>
#include <iostream>

int main() {
  field *f = new field;
  grid *g = new grid;
  PML *p = new PML(f, g);
  Source *s = NULL;
  init_sim_parameters(f, g, s);

  if (check_input_parameters(g)) {
    return 0; // error in input parameters
  }
  array_initialize(f, g, p, s);

  for (std::size_t n = 0; n != maxTime; n++) {
    updateH(f, g, p);
    s->applySource(n);
    updateE(f, g, p);
    std::cout << n + 1 << " / " << maxTime << std::endl;
  }

  printField(f, g);

  //  std::cout << "Hello World" << std::endl;

  // Clean Up Pointers
  delete f;
  delete g;
  delete p;
  delete s;

  return 0;
}
