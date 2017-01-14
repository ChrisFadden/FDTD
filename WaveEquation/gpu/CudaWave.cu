// CUDA Implementation of the 2D wave equation

#include <stdlib.h>
#include <cuda_runtime.h>
#include <vector>
#include <math.h>

__device__ float src(int t, float dt) {
  float freqMag = 1000000000;  // 10^9 = GHz

  return 45 * sin(2 * acosf(-1.0) * 10 * freqMag * dt * t);
}

__global__ void Update(float* UOld, float* UNew, float* U, float* Ca, int SizeX,
                       int SizeY, int isrc, int jsrc, int t, float dt) {
  int r = blockDim.x * blockIdx.x + threadIdx.x;

  int i = r % SizeX;
  int j = (r - i) / SizeX;  // verify this

  int Idx = SizeX * i + j;
  int IdxLeft = SizeX * (i - 1) + j;
  int IdxRight = SizeX * (i + 1) + j;
  int IdxUp = SizeX * i + (j + 1);
  int IdxDown = SizeX * i + (j - 1);

  // Make sure we're in the grid loop
  if (i > 0 && i < (SizeX - 1) && j > 0 && j < (SizeY - 1))
    UNew[Idx] = 2 * U[Idx] - UOld[Idx] +
                Ca[Idx] * (U[IdxRight] + U[IdxLeft] + U[IdxUp] + U[IdxDown] -
                           4 * U[Idx]);

  if (i == isrc && j == isrc) UNew[Idx] = UNew[Idx] + src(t, dt);

  return;
}

// Mur Boundary Condition Kernel
__global__ void ApplyBC(float* UOld, float* UNew, float* U, int SizeX,
                        int SizeY, float dt, float dx) {
  int r = blockDim.x * blockIdx.x + threadIdx.x;

  int i = r % SizeX;
  int j = (r - i) / SizeX;  // verify this

  float cc = 299792458;

  float ABC_C1 = (cc * dt - dx) / (cc * dt + dx);
  float ABC_C2 = 2 * dx / (cc * dt + dx);
  float ABC_C3 = (cc * dt) * (cc * dt) / (2 * dx * (cc * dt + dx));

  int Idx = SizeX * i + j;
  int IdxLeft = SizeX * (i - 1) + j;
  int IdxRight = SizeX * (i + 1) + j;
  int IdxUp = SizeX * i + (j + 1);
  int IdxDown = SizeX * i + (j - 1);

  /**************
   * i == 0
   *************/
  if (i == 0 && j > 0 && j < (SizeY - 1)) {
    UNew[Idx] = -1 * UOld[IdxRight] + ABC_C1 * (UNew[IdxRight] + UOld[Idx]) +
                ABC_C2 * (U[Idx] + U[IdxRight]) +
                ABC_C3 * (U[IdxUp] - 2 * U[Idx] + U[IdxDown] + U[IdxRight + 1] -
                          2 * U[IdxRight] + U[IdxRight - 1]);
  }

  /**************
   * i == SizeX
   *************/
  if (i == SizeX - 1 && j > 0 && j < (SizeY - 1)) {
    UNew[Idx] = -1 * UOld[IdxLeft] + ABC_C1 * (UNew[IdxLeft] + UNew[Idx]) +
                ABC_C2 * (U[Idx] + U[IdxLeft]) +
                ABC_C3 * (U[IdxUp] - 2 * U[Idx] + U[IdxDown] + U[IdxLeft + 1] -
                          2 * U[IdxLeft] + U[IdxLeft - 1]);
  }
  
  /************
   * j == 0
   ***********/
  if(j == 0 && i > 0 && i < SizeX - 1)
  {
    UNew[Idx] = -1 * UOld[IdxUp] + ABC_C1 * (UNew[IdxUp] + UNew[Idx]) +
                ABC_C2 * (U[Idx] + U[IdxUp]) +
                ABC_C3 * (U[IdxRight] - 2*U[Idx] + U[IdxLeft] +
                          U[IdxRight+1] - 2*U[IdxUp] + U[IdxLeft + 1]);

  }

  /**************
   *  j == SizeY
   *************/
  if(j == SizeY-1 && i > 0 && i < SizeX - 1)
  {
    UNew[Idx] = -1 * UOld[IdxDown] + ABC_C1 * (UNew[IdxDown] + UNew[Idx]) +
                ABC_C2 * (U[Idx] + U[IdxDown]) + 
                ABC_C3 * (U[IdxRight] - 2*U[Idx] + U[IdxLeft] + 
                          U[IdxRight-1] - U[IdxDown] + U[IdxLeft - 1]);
  }
  return;
}  // End of ApplyBC function

int main()
{

  /****************
   * Initialization
   ***************/
  int SizeX = 100;
  int SizeY = 100;
  int MaxTime = 1000; 

  int isrc = 50;
  int jsrc = 50;

  float dx = 0.001;
  float cc = 299792458.0;
  float dt = 0.99 / (sqrt(2) * cc);
  
  float caInit = dt*cc / dx;

  std::vector<float> h_OldU(SizeX*SizeY,0.0);
  std::vector<float> h_NewU(SizeX*SizeY,0.0);
  std::vector<float> h_U(SizeX*SizeY,0.0);
  std::vector<float> h_Ca(SizeX*SizeY,caInit*caInit);
  
  float* d_OldU;
  float* d_NewU;
  float* d_U;
  float* d_Ca;
  
  /******************
   *  Allocate Memory
   *****************/
  cudaMalloc((void**)&d_OldU, h_OldU.size());
  cudaMalloc((void**)&d_NewU, h_NewU.size());
  cudaMalloc((void**)&d_U, h_U.size());
  cudaMalloc((void**)&d_Ca, h_Ca.size());

  cudaMemcpy(d_OldU, h_OldU.data(), h_OldU.size(), cudaMemcpyHostToDevice);
  cudaMemcpy(d_NewU, h_NewU.data(), h_NewU.size(), cudaMemcpyHostToDevice);
  cudaMemcpy(d_U, h_U.data(), h_U.size(), cudaMemcpyHostToDevice);
  cudaMemcpy(d_Ca, h_Ca.data(), h_Ca.size(), cudaMemcpyHostToDevice);
  
  int BlockSize = 32;
  int NumBlocks = (SizeX*SizeY - 1) / BlockSize + 1;

  for(int t = 0; t < MaxTime; t++)
  {
    Update<<<NumBlocks,BlockSize>>>(d_OldU,d_NewU,d_U,d_Ca,SizeX,SizeY,isrc,jsrc,t,dt);
    ApplyBC<<<NumBlocks,BlockSize>>>(d_OldU,d_NewU,d_U,SizeX,SizeY,dt,dx); 
    d_OldU = d_U;
    d_U = d_NewU;
  }

  /***************
   *  Free Memory
   **************/
  cudaFree(d_OldU);
  cudaFree(d_NewU);
  cudaFree(d_U);
  cudaFree(d_Ca);

  return 0;
}








