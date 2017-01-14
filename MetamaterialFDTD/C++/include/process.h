#ifndef PROCESS_H
#define PROCESS_H
#include "field.h"
#include "grid.h"
#include "pml.h"
#include "source.h"
#include <string>

void init_sim_parameters(field *f, grid *g, Source *&s, std::string file = "");
void array_initialize(field *f, grid *g, PML *p, Source *s);
void printField(const field *f, const grid *g, std::string fname = "");

#endif
