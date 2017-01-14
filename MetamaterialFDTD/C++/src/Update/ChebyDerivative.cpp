#include <cstdlib>
#include <fftw3.h>
#include <iostream>
#include <math.h>
int main() {

  fftw_plan p;
  fftw_plan p2;
  std::size_t N = 20;
  double *in = (double *)malloc(sizeof(double) * N);
  double *out = (double *)malloc(sizeof(double) * N);
  double *beta = (double *)malloc(sizeof(double) * N);
  double *test = (double *)malloc(sizeof(double) * N);

  for (std::size_t i = 0; i != N; ++i) {
    in[i] = cos(-cos(M_PI * (i + 0.5) / N));
    test[i] = sin(-cos(M_PI * (i + 0.5) / N));
  }

  p = fftw_plan_r2r_1d(N, in, out, FFTW_REDFT10, FFTW_ESTIMATE);
  p2 = fftw_plan_r2r_1d(N, beta, out, FFTW_REDFT01, FFTW_ESTIMATE);

  fftw_execute(p);
  out[0] = out[0] / 2;
  for (std::size_t i = 0; i != N; ++i) {
    out[i] = out[i] / N; // * 2.0 / N;
  }

  beta[N - 1] = 0;
  beta[N - 2] = 2.0 * (N - 2) * out[N - 1];

  for (std::size_t i = N - 2; i != 0; --i) {
    beta[i - 1] = beta[i + 1] + 2.0 * (i)*out[i];
  }
  fftw_execute(p2);
  for (std::size_t i = 0; i != N; ++i) {
    out[i] = out[i] / 2;
  }

  for (std::size_t i = 0; i != N; ++i) {
    std::cout << test[i] - out[i] << std::endl;
  }

  fftw_destroy_plan(p);
  fftw_destroy_plan(p2);
  free(in);
  free(out);

  return 0;
}
