#ifndef _Mersenne_Twister_H
#define _Mersenne_Twister_H

#include <stdio.h>

#define MT_NBITS 32 


class Mersenne_Twister{
// Create a length n array to store the state of the generator
    public:
        const unsigned long int w = 32;
        const unsigned long int n = 624;
        const unsigned long int m = 397;
        const unsigned long int r = 31;
        const unsigned long int a = 0x990bb0df16;
        const unsigned long int u = 11;
        const unsigned long int d = 0xffffffff16;
        const unsigned long int s = 7;
        const unsigned long int b = 0x9D2C568016;
        const unsigned long int t = 15;
        const unsigned long int c = 0xefc6000016;
        const unsigned long int l = 18;
        const unsigned long int f = 1812433253;
        
        unsigned long int State[MT_NBITS];
        void Seed(int se);
        void Twist();
        int Extract();
        
    int index  = MT_NBITS+1;
    const int lower_mask = (1 << r) - 1; // That is, the binary number of r 1's
    const int upper_mask = (~lower_mask) & 0xffffff;  //lowest w bits of (not lower_mask)
 

};

void Mersenne_Twister::Seed(int se){
        int index  = MT_NBITS;
        State[0]  = se;
        for (int i=1;i<MT_NBITS;i++){
            State[i] = (f * (State[i-1] xor (State[i-1] >> (w-2))) + i) & 0xffffff;
        }
}

int Mersenne_Twister::Extract(){
        if(index >= MT_NBITS){
            if(index > n){
            std::cout << "Warning: Mersenne Twister Generator was never seeded."<<std::endl;
            std::cout << "Using default seed of 5489..." <<std::endl;
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
    return y & 0xffffff;
}
 
void Mersenne_Twister::Twist(){
        for(int i=0;i<MT_NBITS;i++){
            int x  = (State[i] & upper_mask)
                   + (State[(i+1) %MT_NBITS] & lower_mask);
            int xA  = x >> 1;
         if ((x%2) != 0)xA  = xA xor a;
         State[i]  = State[(i + m)%MT_NBITS] xor xA;
     }
     index  = 0;
}





// #define RANDOM_BITS 53
// 
// 
// static inline state_t Twist(state_t a, state_t b)
// {
//   state_t mixbits = (a & 0x80000000) | (b & 0x7fffffff);
//   state_t matrixA = (-(b & 1)) & 0x9908b0df;
//   return (mixbits >> 1) ^ matrixA;
// }
// 
// 
// static inline void MersenneReload(state_t *state)
// {
//   state_t *s = state;
//   int j;
// 
//   for( j = MERSENNE_N - MERSENNE_M + 1; --j; ++s )
//     *s = s[MERSENNE_M] ^ Twist(s[0], s[1]);
//   for( j = MERSENNE_M; --j; ++s )
//     *s = s[MERSENNE_M - MERSENNE_N] ^ Twist(s[0], s[1]);
//   *s = s[MERSENNE_M - MERSENNE_N] ^ Twist(s[0], state[0]);
// }
// 
// 
// static inline state_t MersenneInt(state_t s)
// {
//   s ^= s >> 11;
//   s ^= (s << 7) & 0x9d2c5680;
//   s ^= (s << 15) & 0xefc60000;
//   return s ^ (s >> 18);
// }
// 
// 
// static void MersenneGet(This *t, real *x)
// {
//   count next = t->rng.mersenne.next, dim;
// 
//   for( dim = 0; dim < t->ndim; ++dim ) {
// #if RANDOM_BITS == 53
//     state_t a, b;
// #endif
// 
//     if( next >= MERSENNE_N ) {
//       MersenneReload(t->rng.mersenne.state);
//       next = 0;
//     }
// 
// #if RANDOM_BITS == 53
//     a = MersenneInt(t->rng.mersenne.state[next++]) >> 5;
//     b = MersenneInt(t->rng.mersenne.state[next++]) >> 6;
//     x[dim] = (67108864.*a + b)/9007199254740992.;
// #else
//     x[dim] = MersenneInt(t->rng.mersenne.state[next++])/4294967296.;
// #endif
//   }
// 
//   t->rng.mersenne.next = next;
// }
// 
// 
// static void MersenneSkip(This *t, number n)
// {
// #if RANDOM_BITS == 53
//   n = 2*n*t->ndim + t->rng.mersenne.next;
// #else
//   n = n*t->ndim + t->rng.mersenne.next;
// #endif
//   t->rng.mersenne.next = n % MERSENNE_N;
//   n /= MERSENNE_N;
//   while( n-- ) MersenneReload(t->rng.mersenne.state);
// }
// 
// 
// static inline void MersenneIni(This *t)
// {
//   state_t seed = t->seed;
//   state_t *next = t->rng.mersenne.state;
//   count j;
// 
//   for( j = 1; j <= MERSENNE_N; ++j ) {
//     *next++ = seed;
//     seed = 0x6c078965*(seed ^ (seed >> 30)) + j;
//     /* see Knuth TAOCP Vol 2, 3rd Ed, p. 106 for multiplier */
//   }
// 
//   MersenneReload(t->rng.mersenne.state);
//   t->rng.mersenne.next = 0;
// 
//   t->rng.getrandom = MersenneGet;
//   t->rng.skiprandom = MersenneSkip;
// }


#endif
