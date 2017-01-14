#include "Grid.h"
#include <cmath>
#include <fstream>
#include <string>

Grid::Grid() {
  std::string line;
  std::ifstream inputFile("grid.txt");
  std::vector<std::string> parameters(11);
  int i = 0;
  while (std::getline(inputFile, line)) {
    parameters[i] = line;
    i++;
  }

  int paramNum = 0;

  SizeX = std::stoi(parameters[paramNum]);
  paramNum++;
  SizeY = std::stoi(parameters[paramNum]);
  paramNum++;
  dx = std::stod(parameters[paramNum]) * pow(10, -3);
  paramNum++;
  CFL = std::stod(parameters[paramNum]);
  paramNum++;
  MaxTime = std::stoi(parameters[paramNum]);
  paramNum++;
  freq = std::stod(parameters[paramNum]) * pow(10, 9);
  paramNum++;
  isrc = SizeX / 2;
  jsrc = SizeY / 2;

  dt = CFL * dx / (sqrt(2) * cc);

  src = HARMONIC;

  tw = 0.5 / freq;

  t0 = 4.0 * tw;

  src = GAUSSIAN;

  if (std::stoi(parameters[paramNum]) == 0) src = HARMONIC;

  epsr.resize(SizeX * SizeY, 1);
  mur.resize(SizeX * SizeY, 1);
}

int Grid::getMaxTime() { return MaxTime; }

double Grid::HarmonicSource(int t) {
  return 45 * sin(8 * atan(1) * freq * dt * t);
}
double Grid::GaussianSource(int t) {
  return (-2.0 * ((t * dt - t0) / tw) *
          exp(-((t * dt - t0) / tw) * ((t * dt - t0) / tw)));
}

double Grid::Source(int t) {
  switch (src) {
    case HARMONIC:
      return HarmonicSource(t);
      break;
    case GAUSSIAN:
      return GaussianSource(t);
      break;
  }  // end switch
}
