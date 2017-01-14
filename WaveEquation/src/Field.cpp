#include <iostream>
#include <fstream>
#include "Field.h"
#include "Grid.h"

Field::Field(Grid* grid) {
  g = grid;

  uOld.resize(SizeX() * SizeY(), 0);
  u.resize(SizeX() * SizeY(), 0);
  uNew.resize(SizeX() * SizeY(), 0);

  ca.resize(SizeX() * SizeY(), (dt() * cc() / dx()) * (dt() * cc() / dx()));

  for (unsigned i = 0; i != ca.size(); i++) ca[i] = ca[i] / sqrt(Epsr(i) * Mur(i));

  MinValue = 0;
  MaxValue = 0;
  
}

void Field::Update(int t) {

  for (int i = 1; i < SizeX() - 1; i++)
    for (int j = 1; j < SizeY() - 1; j++) {
      UNew(i, j) = 2 * U(i, j) - UOld(i, j) +
                    Ca(i, j) * (U(i + 1, j) + U(i - 1, j) + U(i, j + 1) +
                                U(i, j - 1) - 4 * U(i, j));
    }

  UNew(isrc(), jsrc()) += g->Source(t);

  return;
}

void Field::ReplaceFields()
{
  uOld = u;
  u = uNew;
}

void Field::Print(std::string s) {
  std::ofstream Output;
  Output.open("output/" + s + ".csv");

  for (int i = 0; i < SizeX(); i++)
  {
    for(int j = 0; j < SizeY(); j++)
    {
      Output << U(i,j) << ", ";
      if(U(i,j) > MaxValue) MaxValue = U(i,j);
      if(U(i,j) < MinValue) MinValue = U(i,j);
    }
      Output << std::endl;
  }

  return;
}

void Field::PrintLimits()
{
  std::ofstream Output;
  Output.open("output/Limits.txt");
  Output << MinValue << std::endl;
  Output << MaxValue << std::endl;
}

inline double Field::Dim2(int i, int j) { return (SizeX() * i + j); }
inline double& Field::UOld(int i, int j) { return uOld.at(Dim2(i, j)); }
inline double& Field::UNew(int i, int j) { return uNew.at(Dim2(i, j)); }
inline double& Field::U(int i, int j) { return u.at(Dim2(i, j)); }
inline double& Field::Ca(int i, int j) { return ca.at(Dim2(i, j)); }

inline int Field::SizeX() { return g->SizeX; }
inline int Field::SizeY() { return g->SizeY; }
inline int Field::isrc() { return g->isrc; }
inline int Field::jsrc() { return g->jsrc; }

inline double Field::dt() { return g->dt; }
inline double Field::cc() { return g->cc; }
inline double Field::dx() { return g->dx; }
inline double Field::Epsr(int i) { return g->epsr[i]; }
inline double Field::Mur(int i) { return g->mur[i]; }
