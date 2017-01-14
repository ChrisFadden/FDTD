/** @file softSource.cpp
 *  @brief implementation of a soft source
 *
 *  This implements a soft source, as opposed
 *  to a total field / scattered field source.
 *  This means that it does not inject plane
 *  waves, but rather spherical waves from each
 *  point applied.  Fields reflected from a
 *  scatterer will not be affected by the source.
 */

#include "constants.h"
#include "source.h"
#include <iostream>
#include <math.h>

GaussianSource::GaussianSource(real Freq, grid *Grid) {
  g = Grid;
  tw = 0.5 / Freq;
  t0 = 4 * tw;
}

HarmonicSource::HarmonicSource(real Freq, grid *Grid) {
  freq = Freq;
  g = Grid;
}

/**
 * @brief returns a sinusoid source
 *
 * \f[
 *    f(x) = sin(2 \pi f t + \phi)
 * \f]
 */
real HarmonicSource::calcSource(const std::size_t t, const real offset) {
  return sin(2 * pi * freq * (dt * t - offset));
}

/**
 * @brief returns a gaussian source
 *
 * \f[
 *  f(x) = e^{\frac{(t - t_0 - \phi)}{t_w}^2}
 * \f]
 */
real GaussianSource::calcSource(const std::size_t t, const real offset) {
  return exp(-pow(((t * dt - t0 - offset) / tw), 2));
}

/**
 * @brief applies source to specified grid points
 */
void Source::applySource(const std::size_t t) {
  // Change this to hard source i.e. = instead of += for Gaussian
  Ez(i_src, j_src, k_src) += Cb * source->calcSource(t, 0);
}

void Source::initialize() { return; }
