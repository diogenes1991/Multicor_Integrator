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

thread_local mt19937 rng(0);
#define NTHREADS 8

double RANDOM(double a, double b){
    return uniform_real_distribution<double>{a, b}(rng);
}

double funct(double* a){
    return pow(a[0],5);
}

void EstimateBounds(int ndim, double (*f)(double*), double* bounds){
    double x[ndim];
    for(int i=1;i<=1000;i++){
        for(int j=0;j<ndim;j++) x[j] = RANDOM(0,1);
        double fx = f(x);
        if ( fx > bounds[1]) bounds[1] = fx;
        if ( fx < bounds[0]) bounds[0] = fx;
    }
}

void Integrate(double (*f)(double*), int ndim, double* integral, 
               int batchsize, int maxeval, int verbose, int seed){
    
    double x[ndim];
    rng.seed(seed);
    
    
    
    /// Algorithm to estimate the maxima and minima ///
    for(int j=0;j<ndim;j++) x[j] = 0.5;
    double bounds[2] = {f(x),f(x)};
    EstimateBounds(ndim,f,bounds);
    
    /// Integral initialization ///
    int niter = int(maxeval/batchsize);
    
    
    for(int k=1;k<=niter;k++)
    {
        double loc_min = bounds[0];
        double loc_max = bounds[1];
        
        int count = 0;
        for (int i=1; i<=batchsize; i++)
        {
            for(int j=0;j<ndim;j++) x[j] = RANDOM(0,1);
            double fx = f(x);
            double y = RANDOM(bounds[0],bounds[1]);
            if ( fx > loc_max )   loc_max = fx;
            if ( fx < loc_min )   loc_min = fx;
            if ( fx > y && y > 0 ) count++;
            if ( fx < y && y < 0 ) count--;

        }
        
        double delta;
        if (bounds[0]*bounds[1] > 0){
             delta = bounds[0]+((bounds[1]-bounds[0])*double(count)/batchsize);}
        else delta = ((bounds[1]-bounds[0])*double(count)/batchsize);
        
        integral[0]  += delta;
        integral[1] += pow(delta,2);
        bounds[0] = loc_min;
        bounds[1] = loc_max;
        
        if(verbose>0){
        cout << "Iteration["<<k<<"]: " << k*batchsize;
        cout << " integrand evaluations so far" <<endl;
        if(verbose>1){
            cout << "The bounds for this iteration were = ["<<bounds[0]<<","<<bounds[1]<<"]"<<endl;}
        cout << "Integral = ";
        cout << integral[0]/k << " +- "; 
        cout << sqrt((integral[1]/k - pow(integral[0]/k,2))) << endl;
        cout << endl;
            
        }
        
    }
        
    integral[0] /= niter;
    integral[1] = sqrt((integral[1]/niter - pow(integral[0],2)));
    
}

struct IntegratorArguments{

    double (*Integrand)(double*);
    int SizeofBatches;
    int IntegrandEvaluations;
    int NumberOfVariables;
    double* Integral;
    int VerboseLevel;
    int Seed;
    
};

void LayeredIntegrate(IntegratorArguments IA){
    Integrate(IA.Integrand,IA.NumberOfVariables,IA.Integral,
              IA.SizeofBatches,IA.IntegrandEvaluations,IA.VerboseLevel,IA.Seed);
}

void* ThreadIntegrate(void * IntArgs){
    IntegratorArguments *IA = (IntegratorArguments*)IntArgs;
    LayeredIntegrate(*IA);
    return NULL;
}


void MultithreadIntegrate(IntegratorArguments IA){
    
    //////////////////////////////////////////////////////////////////////////
    ///
    ///   We will have (for now) hard coded batch sizes of 10.000 integrand 
    ///   evaluatios, each core will run 100 batches then a master thread 
    ///   will compile and show the accumulated results for the 100xNTHREADS
    ///   observations so far.
    ///
    //////////////////////////////////////////////////////////////////////////
    
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    
    pthread_t threads[NTHREADS];
    double integral_value[NTHREADS][2];
    IntegratorArguments IntArgs[NTHREADS];
    int rc[NTHREADS];
    
    int ver = IA.VerboseLevel;
    
    for(int i=0;i<NTHREADS;i++){
        IntArgs[i].Integrand = IA.Integrand;
        IntArgs[i].NumberOfVariables = IA.NumberOfVariables;
        IntArgs[i].SizeofBatches = IA.SizeofBatches;
        IntArgs[i].IntegrandEvaluations = IA.IntegrandEvaluations;
        IntArgs[i].VerboseLevel = 0;
        IntArgs[i].Integral = integral_value[i];
        CPU_SET(i, &cpuset);
    }
    
    int NIter = IA.IntegrandEvaluations/(NTHREADS*IA.SizeofBatches);
    IA.Integral[0]=0;
    IA.Integral[1]=0;
    
    for(int j=1;j<=NIter;j++){
        for(int i=0;i<NTHREADS;i++){
            integral_value[i][0]=0;
            integral_value[i][1]=0;
            IntArgs[i].Seed = i*j;
            if (ver)cout << "Now Attempting to create thread "<<i<<endl;
            rc[i] = pthread_create(&threads[i], NULL, ThreadIntegrate,&IntArgs[i]);
            if (rc[i]) {
                cout << "Error:unable to create thread," << rc[i] << endl;
                exit(-1);
            }
            else if(ver) cout << "Thread "<<i<<" has been succesfuly created" << endl;
            
            int set_result = pthread_setaffinity_np(threads[i], sizeof(cpu_set_t), &cpuset);
            if(set_result) cout << "Error: Thread "<<i<<" could not be commited to a new core"<<endl;
            else if (ver) cout << "Thread reassignment succesful" << endl;
        }
        
        for(int i=0;i<NTHREADS;i++){
            pthread_join(threads[i],NULL);
            if(ver){
                cout << "Thread "<<i<<" had the results:" <<endl;
                cout << "Integral = ";
                cout << integral_value[i][0] << " +/- ";
                cout << integral_value[i][1] << endl<<endl;
            }
            IA.Integral[0] += integral_value[i][0];
            IA.Integral[1] += pow(integral_value[i][0],2);
        }
        
        cout << "Iteration["<<j<<"]: " << j*NTHREADS*IA.SizeofBatches;
        cout << " integrand evaluations so far" <<endl;
        cout << "Integral = ";
        cout << IA.Integral[0]/(j*NTHREADS) << " +- "; 
        cout << sqrt(IA.Integral[1]/(NTHREADS*NIter)-pow(IA.Integral[0]/(j*NTHREADS),2)) << endl;
        cout << endl;
        
    }
    
    IA.Integral[0] /= (NTHREADS*NIter);
    IA.Integral[1] = sqrt(IA.Integral[1]/(NTHREADS*NIter)-pow(IA.Integral[0],2));
    
    
}



int main(void)
{
  
   cout.precision(16);
   bool execute_single_core = false;
   bool execute_multi_core = false;
   bool execute_multi_core_2 = false;
   bool execute_multi_core_3 = false;
   bool execute_mersenne_twister = true;
   
   if(execute_mersenne_twister){
       Mersenne_Twister MT1;
       bool ver = true;
       if (ver){
       cout << "Parameters : " << endl;
       cout << "w = " << MT1.w << endl;
       cout << "n = " << MT1.n << endl;
       cout << "m = " << MT1.m << endl;
       cout << "r = " << MT1.r << endl;
       cout << "a = " << MT1.a << endl;
       cout << "u = " << MT1.u << endl;
       cout << "d = " << MT1.d << endl;
       cout << "s = " << MT1.s << endl;
       cout << "b = " << MT1.b << endl;
       cout << "t = " << MT1.t << endl;
       cout << "c = " << MT1.c << endl;
       cout << "l = " << MT1.l << endl;
       }
//        MT1.Seed(time(NULL));
       cout << "MAX BIT DEPTH = " << (1<<(MT1.w-1))-1 << endl;
       int NI = 1000;
       double AVG1 = 0;
       double AVG2 = 0;
       double AVG3 = 0;
       double AVG4 = 0;
       for(int i=0;i<NI;i++){
           unsigned long int aux = MT1.Extract();
           AVG1 += MT1.Extract();
           AVG2 += (MT1.Extract()+MT1.Extract())/2;
           AVG3 += (MT1.Extract()+MT1.Extract()+MT1.Extract())/3;
           AVG4 += (MT1.Extract()+MT1.Extract()+MT1.Extract()+MT1.Extract())/4;
           
       }
       AVG1/= NI;
       AVG2/= NI;
       AVG3/= NI;
       AVG4/= NI;
       cout << "<x^1> = " << AVG1 << endl;
       cout << "<x^2> = " << AVG2 << endl;
       cout << "<x^3> = " << AVG3 << endl;
       cout << "<x^4> = " << AVG4 << endl;
          
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
}
   
   ///////////////////////////////////////////////////////////////////////////
   ///
   ///          Single Thread Execution 
   ///   This is a dummy example of the creation of a single thread 
   ///   This piece is equivalent to the simpler call for 
   ///   LayeredIntegrate(IntArg0);
   ///   The reason to do this more sophisticated call is to show the pthread
   ///   call in a way that is generalizable to more threads.
   ///
   ///////////////////////////////////////////////////////////////////////////
   
    if(execute_single_core){
    pthread_t thr0;
    double integral_value0[2] = {0,0};
    IntegratorArguments IntArg0;
    IntArg0.Integrand = funct;
    IntArg0.IntegrandEvaluations = 10000000;
    IntArg0.SizeofBatches = 500000;
    IntArg0.NumberOfVariables = 1;
    IntArg0.VerboseLevel = 2;
    IntArg0.Seed = 5;
   
    IntArg0.Integral = integral_value0;
    int t = time(NULL);
    cout << "Now Attempting to create thread "<<0<<endl;
    int rc0 = 0;
    rc0 = pthread_create(&thr0, NULL, ThreadIntegrate,&IntArg0);   
    if (rc0) {
        cout << "Error:unable to create thread," << rc0 << endl;
        exit(-1);
    }
    else cout << "Thread "<<0<<" has been succesfuly created" << endl;
    pthread_join(thr0,NULL);
    cout << "Thread 0 has finished, it took " << time(NULL)-t <<" secs to finish"  << endl;
    cout << "Integral Value = "<< integral_value0[0] << " +/- " << integral_value0[1] <<endl;
    }
    
    
    ////////////////////////////////////////////////////////////////////////////////
    ///
    ///             Multiple Threads Creation 
    ///
    ///   This implementation works, however it creates all the cores on the same CPU.
    ///   If you top while executing you'll see a CPU has a load over 100% if NTHREADS>1.
    ///   Since this CPU has to handle more load than it can it ends up taking longer 
    ///   than sequantial execution of NTHREADS processes.
    ///   This is obviously not what we want, we want concurrency 
    ///   all the processes to fully occupy different CPUS.
    ///   
    /////////////////////////////////////////////////////////////////////////////// 
    
    if(execute_multi_core){
        
    pthread_t threads[NTHREADS];
    double integral_value[NTHREADS][2];
    IntegratorArguments IntArgs[NTHREADS];
    int rc[NTHREADS];
    for(int i=0;i<NTHREADS;i++){
        integral_value[i][0]=0;
        integral_value[i][1]=0;
        IntArgs[i].Integrand = funct;
        IntArgs[i].IntegrandEvaluations = 10000000;
        IntArgs[i].SizeofBatches = 50000;
        IntArgs[i].NumberOfVariables = 2;
        IntArgs[i].VerboseLevel = 0;
        IntArgs[i].Seed = i;
        IntArgs[i].Integral = integral_value[i];        
    }
        
    int t = time(NULL);
    for(int i=0;i<NTHREADS;i++){
        cout << "Now Attempting to create thread "<<i<<endl;
        rc[i] = pthread_create(&threads[i], NULL, ThreadIntegrate,&IntArgs[i]);
        if (rc[i]) {
            cout << "Error:unable to create thread," << rc[i] << endl;
            exit(-1);
        }
        else cout << "Thread "<<i<<" has been succesfuly created" << endl;
    }
    /// Thread Waiting Phase ///
    for(int i=0;i<NTHREADS;i++) pthread_join(threads[i],NULL);
    cout << "All threads have now finished" <<endl;
    cout << "This took " << time(NULL)-t << " secs to finish" <<endl;
    cout << "Or " << (time(NULL)-t)/NTHREADS << " secs per core" <<endl;
    for(int i = 0; i < NTHREADS; i++ ) {
        cout << "Thread " << i << " has as the value for the integral" << endl;
        cout << "Integral = ";
        cout << integral_value[i][0] << " +- "; 
        cout << integral_value[i][1] << endl;
    }
      
    }
    
    ////////////////////////////////////////////////////////////////////////
    ///
    ///             Multiple Cores Execution 
    ///   Here we attempt to commit different threads to different cores,
    ///   we use pthread_setaffinity_np to do this.
    ///
    ///////////////////////////////////////////////////////////////////////
    

    if(execute_multi_core_2){
        
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    
    pthread_t threads[NTHREADS];
    double integral_value[NTHREADS][2];
    IntegratorArguments IntArgs[NTHREADS];
    int rc[NTHREADS];
    for(int i=0;i<NTHREADS;i++){
        integral_value[i][0]=0;
        integral_value[i][1]=0;
        IntArgs[i].Integrand = funct;
        IntArgs[i].IntegrandEvaluations = 100000000;
        IntArgs[i].SizeofBatches = 5000;
        IntArgs[i].NumberOfVariables = 1;
        IntArgs[i].VerboseLevel = 0;
        IntArgs[i].Seed = i;
        IntArgs[i].Integral = integral_value[i];        
    }
        
    int t = time(NULL);
    for(int i=0;i<NTHREADS;i++){
        cout << "Now Attempting to create thread "<<i<<endl;
        rc[i] = pthread_create(&threads[i], NULL, ThreadIntegrate,&IntArgs[i]);
        if (rc[i]) {
            cout << "Error:unable to create thread," << rc[i] << endl;
            exit(-1);
        }
        else cout << "Thread "<<i<<" has been succesfuly created" << endl;
        CPU_SET(i, &cpuset);
    }
    
    cout << "Now attempting to commit different threads to different cores" << endl;
    for(int i=0;i<NTHREADS;i++){
        const int set_result = pthread_setaffinity_np(threads[i], sizeof(cpu_set_t), &cpuset);
        if(set_result) cout << "Error: Thread "<<i<<" could not be commited to a new core"<<endl;
        else cout << "Thread reassignment succesful" << endl;
    }
    
    /// Thread Waiting Phase ///
    for(int i=0;i<NTHREADS;i++) pthread_join(threads[i],NULL);
    cout << "All threads have now finished" <<endl;
    cout << "This took " << time(NULL)-t << " secs to finish" <<endl;
    cout << "Or " << (time(NULL)-t)/NTHREADS << " secs per core" <<endl;
    for(int i = 0; i < NTHREADS; i++ ) {
        cout << "Thread " << i << " has as the value for the integral" << endl;
        cout << "Integral = ";
        cout << integral_value[i][0] << " +- "; 
        cout << integral_value[i][1] << endl;
    }
      
    }
    
    ////////////////////////////////////////////////////////////////////////////
    ///               
    ///                NEW API No-Common pointers will be used now
    ///
    ///////////////////////////////////////////////////////////////////////////
    
    if(execute_multi_core_3){
        
        
    double integral_value0[2] = {0,0};
    IntegratorArguments IntArg0;
    IntArg0.Integral = integral_value0;
    IntArg0.Integrand = funct;
    IntArg0.IntegrandEvaluations = 10000000;
    IntArg0.SizeofBatches = 100000;
    IntArg0.NumberOfVariables = 1;
    IntArg0.VerboseLevel = 0;
    IntArg0.Seed = 5;
   
    int t = time(NULL);
    MultithreadIntegrate(IntArg0);
    cout << "Integration has finished, it took " << time(NULL)-t <<" secs to finish"  << endl;
    cout << "Integral Value = "<< integral_value0[0] << " +/- " << integral_value0[1] <<endl;
    
    
    
    }
    
    
    
pthread_exit(NULL);
    
}
