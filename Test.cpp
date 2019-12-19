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

thread_local mt19937 rng(0);

double aleatorio(double a, double b){
    return uniform_real_distribution<double>{a, b}(rng);
}

double funct(double* a){
    return pow(a[0],6);
}

void EstimateBounds(int ndim, double (*f)(double*), double* bounds){
    double x[ndim];
    for(int i=1;i<=1000;i++){
        for(int j=0;j<ndim;j++) x[j] = aleatorio(0,1);
        double fx = f(x);
        if ( fx > bounds[1]) bounds[1] = fx;
        if ( fx < bounds[0]) bounds[0] = fx;
    }
}

void Integrate(double (*f)(double*), int ndim, double* integral, int verbose, int seed){
    
    int nbatch = 5000000;
    const int maxeval = 25*nbatch;
    double x[ndim];
    rng.seed(seed);
    
    
    
    /// Algorithm to estimate the maxima and minima ///
    /// Plug and Play ///
    for(int j=0;j<ndim;j++) x[j] = 0.5;
    double bounds[2] = {f(x),f(x)};
    EstimateBounds(ndim,f,bounds);
    
    /// Integral initialization ///
    int niter = int(maxeval/nbatch);
    
    
    for(int k=1;k<=niter;k++)
    {
        
        
        double loc_min = bounds[0];
        double loc_max = bounds[1];
        
        int count = 0;
        for (int i=1; i<=nbatch; i++)
        {
            for(int j=0;j<ndim;j++) x[j] = aleatorio(0,1);
            double fx = f(x);
            double y = aleatorio(bounds[0],bounds[1]);
            if ( fx > loc_max )   loc_max = fx;
            if ( fx < loc_min )   loc_min = fx;
            if ( fx > y && y > 0 ) count++;
            if ( fx < y && y < 0 ) count--;

        }
        
        double delta = (bounds[1]-bounds[0])*double(count)/nbatch;
        integral[0]  += delta;
        integral[1] += pow(delta,2);
        bounds[0] = loc_min;
        bounds[1] = loc_max;
        
        if(verbose>0){
        cout << "Iteration["<<k<<"]: " << k*nbatch;
        cout << " integrand evaluations so far" <<endl;
        if(verbose>1){
            cout << "The bounds for this iteration were = ["<<bounds[0]<<","<<bounds[1]<<"]"<<endl;}
        cout << "Integral = ";
        cout << integral[0]/k << " +- "; 
        cout << sqrt((integral[1]/k - pow(integral[0]/k,2)))/(k) << endl;
        cout << endl;
        }
        
    }
    integral[0] /= niter;
    integral[1] = sqrt((integral[1]/niter - pow(integral[0],2)))/niter;
    
}

struct IntegratorArguments{

    double (*Integrand)(double*);
    int NumberOfVariables;
    double* Integral;
    int VerboseLevel;
    int Seed;
    
};

void LayeredIntegrate(IntegratorArguments IA){
    Integrate(IA.Integrand,IA.NumberOfVariables,IA.Integral,IA.VerboseLevel,IA.Seed);
}

void * ThreadIntegrate(void * IntArgs){
    IntegratorArguments *IA = (IntegratorArguments*)IntArgs;
    LayeredIntegrate(*IA);
    return NULL;
    pthread_exit(NULL);
    
}

#define NTHREADS 4

int main(void)
{
  
   cout.precision(16);
   bool execute_single_core = true;
   bool execute_multi_core = true;
   bool execute_multi_core_2 = true;
   bool execute_multi_core_3 = false;
   
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
    IntArg0.NumberOfVariables = 2;
    IntArg0.VerboseLevel = 0;
    IntArg0.Seed = 1;
   
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
    cout << "Integral Value = "<< integral_value0[0] << "+/-" << integral_value0[1] <<endl;
    }
    
    
    ////////////////////////////////////////////////////////////////////////////////
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
        
    cout << "Manual Thread creation and core commitment" <<endl;
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    
    pthread_t thr1,thr2;
    double integral_value1[2] = {0,0};
    double integral_value2[2] = {0,0};
    IntegratorArguments IntArgs1,IntArgs2;
    int rc1,rc2;
    
    IntArgs1.Integrand = funct;
    IntArgs2.Integrand = funct;
    IntArgs1.NumberOfVariables = 2;
    IntArgs2.NumberOfVariables = 2;
    IntArgs1.VerboseLevel = 0;
    IntArgs2.VerboseLevel = 0;
    IntArgs1.Seed = 1;
    IntArgs2.Seed = 2;
    IntArgs1.Integral = integral_value1;        
    IntArgs2.Integral = integral_value2;        
        
    int t = time(NULL);
    
    cout << "Now Attempting to create thread "<<1<<endl;
    rc1 = pthread_create(&thr1, NULL, ThreadIntegrate,&IntArgs1);
    if (rc1) {
        cout << "Error:unable to create thread," << rc1 << endl;
        exit(-1);
    }
    else cout << "Thread "<<1<<" has been succesfuly created" << endl;
    CPU_SET(1, &cpuset);
    
    cout << "Now Attempting to create thread "<<2<<endl;
    rc2 = pthread_create(&thr2, NULL, ThreadIntegrate,&IntArgs2);
    if (rc2) {
        cout << "Error:unable to create thread," << rc2 << endl;
        exit(-1);
    }
    else cout << "Thread "<<2<<" has been succesfuly created" << endl;
    CPU_SET(2, &cpuset);
    
    cout << "Now attempting to commit different threads to different cores" << endl;
    
    const int set_result1 = pthread_setaffinity_np(thr1, sizeof(cpu_set_t), &cpuset);
    if(set_result1) cout << "Error: Thread "<<1<<" could not be commited to a new core"<<endl;
    else cout << "Thread reassignment succesful" << endl;
    
    const int set_result2 = pthread_setaffinity_np(thr2, sizeof(cpu_set_t), &cpuset);
    if(set_result2) cout << "Error: Thread "<<2<<" could not be commited to a new core"<<endl;
    else cout << "Thread reassignment succesful" << endl;
    
    /// Thread Waiting Phase ///
    pthread_join(thr1,NULL);
    pthread_join(thr2,NULL);
    cout << "All threads have now finished" <<endl;
    cout << "This took " << time(NULL)-t << " secs to finish" <<endl;
    
    cout << "Or " << (time(NULL)-t)/NTHREADS << " secs per core" <<endl;
    cout << "Thread " << 1 << " has as the value for the integral" << endl;
    cout << "Integral = ";
    cout << integral_value1[0] << " +- "; 
    cout << integral_value1[1] << endl;
    
    cout << "Thread " << 2 << " has as the value for the integral" << endl;
    cout << "Integral = ";
    cout << integral_value2[0] << " +- "; 
    cout << integral_value2[1] << endl;
      
    }
    
    
    
pthread_exit(NULL);
    
}
