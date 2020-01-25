#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "pcg_basic.h"

#define D 100
#define UP 1
#define DOWN -1
#define S(x,y) spin[D*x + y]


/*RNG selection macros taken from mat_gen.c by CDHW.*/
#define NO 0
#define YES 1
#define USE_PCG YES
#define SEED_TIME YES

#if ( USE_PCG ) /* Bad random number generator */
#  define VAR_TYPE pcg32_random_t
#  define RANDOM_MAX 2 /* Maximum value from random(3) man page */
#  define RANDOM_BOUND_0_1(A,B) pcg32_boundedrand_r(A,B)
#  define SRANDOM_R(X,Y,Z) pcg32_srandom_r(X,Y,Z)
#else /* Better random number generator */
#  define VAR_TYPE double
#  define RANDOM_MAX RAND_MAX
#  define RANDOM_BOUND_0_1(X,Y) random_bound_0_1()
#  define SRANDOM_R(X,Y,Z) srand(Y)
#endif

#if ( SEED_TIME)
#  define SEED time(NULL)
#else
#  define SEED 42u
#endif

/*Program to model time evolution of the Ising Model of Ferromagnetism





*/




/*function to create and initialise a random 2D array of spins*/
static void initialize_spin_array(){
    /*map 2D array onto 1D array as its faster*/
    /*elements accessed as element_ij = spin_array(D*i+j)*/
    int *spin_array = malloc(D*D*sizeof(int));



}

/*wrapping function for default rand() rng
any rng without inbuilt capability to generate integers in a specific
region must have a similar wrapper written and combined with
appropriate macros.  */
static int random_bound_0_1(){
    return rand()% 2;
}


int main()
{

    int rounds = 5;
    int bound = 2;
    int i = 0;
    VAR_TYPE rng1;
    SRANDOM_R(&rng1, SEED, (intptr_t)&rng1);

    printf("%d\n",rng1);

        for(i = 0; i< rounds; i++){
            printf("bounded rn %d is %ld \n", i, RANDOM_BOUND_0_1(&rng1, bound));
        }


    printf("Hello world!\n");
    int x = 99;
    int r = (D+(x))%D;
    printf("%d\n",r);
    printf("%d\n", RANDOM_MAX);
    return 0;
}
