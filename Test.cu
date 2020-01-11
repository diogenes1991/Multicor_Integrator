#include <iostream>
#include <sched.h>
#include <pthread.h>
#include <thread>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <random> 

using namespace std;

thread_local mt19937 rng(5);

double aleatorio(double a, double b){
    return uniform_real_distribution<double>{a, b}(rng);
}



// Kernel function to add the elements of two arrays
__global__
void add(int n, float *x, float *y)
{
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
    /// We wish to integrate x^5 from 0 to 1
    return 1+pow(a[0],5);
}

__global__ void integrand(double *x,double *y,int *c){
    double fx = funct(x);
    if ( fx > y[threadIdx.x] && y[threadIdx.x] > 0 ) c[threadIdx.x]+=1;
    if ( fx < y[threadIdx.x] && y[threadIdx.x] < 0 ) c[threadIdx.x]-=1;
}

int main(void)
{
  int N = 1<<10;
  int *c;
  double *x,*y;
  cudaMallocManaged(&x,N*sizeof(double));
  cudaMallocManaged(&y,N*sizeof(double));
  cudaMallocManaged(&c,N*sizeof(double));
  
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
  
  return 0;
  
}
