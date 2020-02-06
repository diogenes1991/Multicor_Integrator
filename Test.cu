#include <iostream>
#include <sched.h>
#include <pthread.h>
#include <thread>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <random> 
#include "Mersenne_Twister.h"

using namespace std;

thread_local mt19937 rng(5);

class CUDA_Mersenne_Twister{
// Create a length n array to store the state of the generator
    public:
        const int w = 32;
        const  int n = 624;
        const  int m = 397;
        const  int r = 31;
        const unsigned long int a = 0x990bb0df16;
        const  int u = 11;
        const unsigned long  int d = 0xffffffff16;
        const  int s = 7;
        const unsigned long int b = 0x9D2C568016;
        const  int t = 15;
        const unsigned long int c = 0xefc6000016;
        const  int l = 18;
        const unsigned long int f = 1812433253;
        
        unsigned long long int State[MT_NBITS];
        __device__ void Seed(int se);
        __device__ void Twist();
        __device__ unsigned long long int Extract();
    
    int index  = MT_NBITS+1;
    const int lower_mask = (1 << r) - 1; // That is, the binary number of r 1's
    const int upper_mask = (~lower_mask) & ((1<<(w-1))-1);  //lowest w bits of (not lower_mask)
 

};

// __device__ CUDA_Mersenne_Twister(){
//    CUDA_Mersenne_Twister::Seed(threadIdx) ;
// }

__device__ void CUDA_Mersenne_Twister::Seed(int se){
        index  = MT_NBITS;
        State[0]  = se;
        for (int i=1;i<MT_NBITS;i++){
            State[i] = (f * (State[i-1] xor (State[i-1] >> (w-2))) + i) & ((1<<(w-1))-1);
        }
}

__device__ unsigned long long int CUDA_Mersenne_Twister::Extract(){
        if(index >= MT_NBITS){
            if(index > n){
//             std::cout << "Warning: Mersenne Twister Generator was never seeded."<<std::endl;
//             std::cout << "Using default seed of 5489..." <<std::endl;
            Seed(5489);
            }
            Twist();
        }

    int y  = State[index];
    y  = y xor ((y >> u) & d);
    y  = y xor ((y << s) & b);
    y  = y xor ((y << t) & c);
    y  = y xor (y >> l);
 
    index  = index + 1;
    return y & ((1<<(w-1))-1);
}
 
__device__ void CUDA_Mersenne_Twister::Twist(){
        for(int i=0;i<MT_NBITS;i++){
            int x  = (State[i] & upper_mask)
                   + (State[(i+1) %MT_NBITS] & lower_mask);
            int xA  = x >> 1;
         if ((x%2) != 0)xA  = xA xor a;
         State[i]  = State[(i + m)%MT_NBITS] xor xA;
     }
     index  = 0;
}

__global__ void CUDA_RAND(){
    CUDA_Mersenne_Twister MT1;
    
}

// class Dummy_Mersenne_Twister{
// // Create a length n array to store the state of the generator
//     public:
//         __device__ int dummy_extract();
//  
// 
// };
// 
// __device__ int Dummy_Mersenne_Twister::dummy_extract(){
//     return 7;
// }

double aleatorio(double a, double b){
    return uniform_real_distribution<double>{a, b}(rng);
}



// Kernel function to add the elements of two arrays
__global__
void add(int n, float *x, float *y){
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;
  for (int i = index; i < n; i += stride)
    y[i] = x[i] + y[i];
}


////////////////////////////////////////////////////////////////
///
///             GPU Compiled MC Integrator: 
///  The basic idea is to have a function that performs the 
///  calls to the underlying integrand using the GPU.
///  for simplicity we will start with a MC Integrator that does a 
///  single MC Integration with 1000 evaluations.
///
///
////////////////////////////////////////////////////////////////

__device__ double funct(double*a){
    return 1+pow(a[0],5);
}

// __device__ int CUDAMT(){
//     Mersenne_Twister MT1;
//     return MT1.Extract();
//     
// }


__global__ void integrand(double *x,double *y,int *c){
    double fx = funct(x);
    if ( fx > y[threadIdx.x] && y[threadIdx.x] > 0 ) c[threadIdx.x]+=1;
    if ( fx < y[threadIdx.x] && y[threadIdx.x] < 0 ) c[threadIdx.x]-=1;
}

__global__ void Dummy_Init_Mersenne(int N, int* x, int* y, bool* retval){
    retval[threadIdx.x] = true;
    for (int i=0;i<N;i++){
    if(x[threadIdx.x]<y[i]) break;
    if((x[threadIdx.x]%y[i])==0) retval[threadIdx.x] = false;
    break;
    }
    
}

__global__ void CUDA_Extract(CUDA_Mersenne_Twister* Gen, int* retval){
//   return 0;
}


int main(void)
{
  int N = 1<<10;
  int *x,*y;
  bool* r;
  cudaMallocManaged(&x,N*sizeof(int));
  cudaMallocManaged(&y,N*sizeof(int));
  cudaMallocManaged(&r,N*sizeof(bool));
  
  for(int i=0;i<N;i++){
      x[i] = 2*i+3;
      y[i] = 2*i+3;
}

  int t = time(NULL);
  Dummy_Init_Mersenne<<<1,N>>>(N,x,y,r);
  cudaDeviceSynchronize();
  for(int i=0;i<N;i++){
      if (r[i]){
      std::cout << "We have that " << x[i] << ":"
      << (r[i] ? "  is":"  is not") << " prime." 
      <<std::endl;
      }
  }

  cout << "This took " << time(NULL)-t << " secs to finish" <<endl;
  
  cudaFree(r);
  cudaFree(x);
  cudaFree(y);
  
  
  
  CUDA_Mersenne_Twister *Gen;
  
  cudaMallocManaged(&Gen,N*sizeof(CUDA_Mersenne_Twister));  
  cout << "Size of our CUDA MT = " << sizeof(CUDA_Mersenne_Twister) << endl;
  Gen = new CUDA_Mersenne_Twister[N];
  
  
  
  
  
  
  /*
  int *c;
  double *x,*y;
  cudaMallocManaged(&x,N*sizeof(double));
  cudaMallocManaged(&y,N*sizeof(double));
  cudaMallocManaged(&c,N*sizeof(int));
  
  
  
  
  int BatchSize = 10000;
  int NIterations = 10;
  
  cout.precision(16);
  double integral[2]={0,0};
  int t = time(NULL);
  int GLOBAL_COUNTER = 0;
    for(int k=1;k<=NIterations;k++){
        int count = 0;
        for(int i=0;i<N;i++)c[i]=0;
        for(int j=0;j<BatchSize;j++){
            for(int i=0;i<N;i++){x[i]=aleatorio(0,1);y[i]=aleatorio(0,2);}
            integrand<<<1,N>>>(x,y,c);
            cudaDeviceSynchronize();
        }
        for(int i=0;i<N;i++){
            count += c[i];
//             cout << "Kernel #"<<i<<" has local counter of: "<<c[i]<<endl;
        }
        GLOBAL_COUNTER += count;
//         cout << "Accumulated counter of : " << count <<endl;
        cout << "Iteration["<<k<<"]: " << k*BatchSize*N;
        cout << " integrand evaluations so far" <<endl;
        
        
        double delta = double(count)/(N*BatchSize);
        integral[0]  += delta;
        integral[1]  += pow(delta,2);
        cout << "Integral = ";
        cout << integral[0]/k << " +/- "; 
        cout << sqrt((integral[1]/k - pow(integral[0]/k,2))) << endl;
        cout << endl;
    }
    
    cout << "The global counter is = " << GLOBAL_COUNTER <<endl;
    cout << "Predicting  = " << double(GLOBAL_COUNTER) / (N*BatchSize*NIterations)<<endl;
    cout << "This took " << time(NULL)-t << " secs to finish" <<endl;
//     cout << "Total counenter  = " << count
    
  cudaFree(c);
 */ 
  return 0;
  
}
