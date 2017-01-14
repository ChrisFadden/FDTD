#ifndef ERROR_H
#define ERROR_H
#include "grid.h"
#include <string>
void error_msg(std::string);
void warning_msg(std::string);
bool check_input_parameters(grid *g);
#endif
