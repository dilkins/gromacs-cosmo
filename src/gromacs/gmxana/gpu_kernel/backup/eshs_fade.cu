#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "sysstuff.h"
#include "typedefs.h"
#include "macros.h"
#include "vec.h"
#include "pbc.h"
#include "xvgr.h"
#include "gromacs/fileio/futil.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trxio.h"
#include "physics.h"
#include "index.h"
#include "gromacs/utility/smalloc.h"
#include "calcgrid.h"
#include "nrnb.h"
#include "coulomb.h"
#include "../gstat.h"
#include "gromacs/fileio/matio.h"
#include "../gmx_ana.h"
#include "../hyperpol.h"
#include "names.h"
#include "cuda.h"
#include "gromacs/legacyheaders/gmx_fatal.h"
#include "gromacs/gmxlib/cuda_tools/cudautils.cuh"


#define block_size 1024  // number of threads per block for kernel launching

__global__ 
void kernel_fade(t_pbc *pbc,real beta_i,real *beta,rvec x_i,rvec *x,real rmax2,rvec qvec,real *result,int n,int N,
                 real fade,real inv_width) {

  __shared__ real buffer[block_size];
  int j,k;
  rvec dx;
  real r,r2,temp=0;

  // Obtain the thread ID
  j=threadIdx.x+blockIdx.x*blockDim.x;

  // Prevent the unnecessary threads
  if (j>=N)  return;

  // Initialize the shared memory in each Block
  //buffer[threadIdx.x]=0.0;

  // Calculate vector from the ith to jth molecule
  dx[XX]=(x[j][XX]-x_i[XX]);
  dx[YY]=(x[j][YY]-x_i[YY]);
  dx[ZZ]=(x[j][ZZ]-x_i[ZZ]);

  // Apply PBC condition
  dx[XX]-=rintf(dx[XX]/pbc->box[XX][XX])*pbc->box[XX][XX];
  dx[YY]-=rintf(dx[YY]/pbc->box[YY][YY])*pbc->box[YY][YY];
  dx[ZZ]-=rintf(dx[ZZ]/pbc->box[ZZ][ZZ])*pbc->box[ZZ][ZZ];

  // Distance sqaure between ith and jth
  r2=(dx[XX]*dx[XX]+dx[YY]*dx[YY]+dx[ZZ]*dx[ZZ]);
  
  if (r2>0 && r2<=rmax2) {
    r=sqrt(r2);

    if (r <= fade) {
      temp=beta_i*beta[j]*cos(qvec[XX]*dx[XX]+qvec[YY]*dx[YY]+qvec[ZZ]*dx[ZZ]);
    }
    else {
       temp=beta_i*beta[j]*cos(qvec[XX]*dx[XX]+qvec[YY]*dx[YY]+qvec[ZZ]*dx[ZZ])*sqr(cos((r-fade)*inv_width));
    }
  }

  // Save each element to the shared memory
  buffer[threadIdx.x]=temp; 

  // Synchronize all the threads
  __syncthreads(); 

  // Summation Reduction
  k=blockDim.x*0.5;

  while (k!=0) {
    if (threadIdx.x<k) {
      buffer[threadIdx.x]+=buffer[threadIdx.x+k];
    }
    __syncthreads();
    k*=0.5;
  }

  // Return the first element which contains the sum of the block
  if (threadIdx.x==0) {
    result[blockIdx.x]=buffer[0];
    //printf("result[%d][%d] = %f\n",blockIdx.x,blockIdx.y,buffer[0]);
  }
}


extern void double_sum_fade(t_pbc *pbc,real *beta,rvec *x,int n,real rmax2,rvec *qvec,int nbinq,real fade,real inv_width,real *temp_method) {

  int i,j,q,N=0,grid_size=0;
  t_pbc *pbc_d=NULL;
  rvec *x_d=NULL,*qvec_d=NULL;
  real *beta_d=NULL,*result=NULL,*result_d=NULL,*Sum=NULL;
  //cudaError_t status;

  // Determine grid_size of the kernel
  grid_size=(int)(n+block_size-1)/block_size;

  // Determine total number of threads
  N=block_size*grid_size;

  // Allocate Host Memory
  result=(real *)calloc(grid_size,sizeof(real));
  Sum=(real *)calloc(nbinq,sizeof(real));  

  // Allocate Device Memory
  cudaMalloc(&pbc_d,sizeof(t_pbc));
  cudaMalloc(&x_d,sizeof(rvec)*n);
  cudaMalloc(&qvec_d,sizeof(rvec)*nbinq);
  cudaMalloc(&beta_d,sizeof(real)*n);
  cudaMalloc(&result_d,sizeof(real)*grid_size);

  // Copy memory from Host to Device
  cudaMemcpy(pbc_d,pbc,sizeof(t_pbc),cudaMemcpyHostToDevice);
  cudaMemcpy(x_d,x,sizeof(rvec)*n,cudaMemcpyHostToDevice);
  cudaMemcpy(qvec_d,qvec,sizeof(rvec)*nbinq,cudaMemcpyHostToDevice);
  cudaMemcpy(beta_d,beta,sizeof(real)*n,cudaMemcpyHostToDevice);

  // Calculate individual water in series
  for (i=0;i<n;i++) {  
    for (q=0;q<nbinq;q++) {
      // Launch kernel with N threads
      kernel_fade<<<grid_size,block_size>>>(pbc_d,beta[i],beta_d,x_d[i],x_d,rmax2,qvec_d[q],result_d,n,N,fade,inv_width);
      //status=cudaGetLastError();
      //printf("kernel: %s\n",cudaGetErrorString(status));

      // Copy result back from Device to Host 
      //status=
      cudaMemcpy(result,result_d,sizeof(real)*grid_size,cudaMemcpyDeviceToHost);
      //printf("copy D2H %s\n",cudaGetErrorString(status)); 

      // Sum over the individual sum in each block
      for (j=0;j<grid_size;j++) {      
        Sum[q]+=result[j];
        //printf("block = %d, q = %d, result = %f\n",j,q,result[j+grid_size*q]);
      }
    }
  }

  // Send the number back
  for (q=0;q<nbinq;q++) {
    temp_method[q]=Sum[q]*0.5;
  }

  // Free Device Memory
  cudaFree(pbc_d);
  cudaFree(x_d);   
  cudaFree(beta_d); 
  cudaFree(result_d);

  // Free Host Memory
  free(result);

  return;
}


