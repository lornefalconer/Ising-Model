#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "pcg_basic.h"

#define MAXLENGTH 256
#define D 100
#define UP 1
#define DOWN -1
#define S(x,y) spin_array[D*((D+(x))%D) + ((D+(y))%D)]
#define J 1.0e-21
#define MU0 12.57e-7
#define B 0.01

#define MEM_ERROR 3


#define GNUPLOT "gnuplot"
#define SCRIPT "spin_array_gnuplot.script"

/*RNG selection macros taken from mat_gen.c by CDHW.*/
#define NO 0
#define YES 1
#define USE_PCG YES
#define SEED_TIME YES
#define BOUND 2
#define TESTING 0

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


/*function to open any file and perform error check*/
/*taken directly from my n_body.c program*/
void err_check_file(FILE * fptr){
    if (fptr == NULL){
        printf("file cannot be found\n");
        exit(MEM_ERROR);
        }
}


/*Program to model time evolution of the Ising Model of Ferromagnetism





*/



/*function to create and initialise a random 2D array of spins*/
void initialize_spin_array(int* spin_array){
    /*map 2D array onto 1D array as its faster*/
    /*elements accessed as element_ij = spin_array(D*i+j)*/
    VAR_TYPE rng1;
    SRANDOM_R(&rng1, SEED, (intptr_t)&rng1);




    for(int x = 0; x<D; x++){
        for(int y = 0; y<D; y++){
            int result = RANDOM_BOUND_0_1(&rng1, BOUND);

            if(result == 0){
                result = -1;
            }
            spin_array[D*x +y] = result;

            if(TESTING){/*Test initialization.*/
                printf("TEST spin element_%d%d is: %d\n",x,y,spin_array[D*x +y]);

            }
        }


    }



}

/*function to step through the apin array and write to file.
Wrtten as a separate function to the current state of the spin can be analysed during
MC iteration.*/
static void spin_array_write_to_file(int* spin_array){
    FILE *fspin = fopen("spin_array.txt", "w+");/*opens file to write the ith body data to*/
    err_check_file(fspin);
    fprintf(fspin,"#x\ty\tS\n");

    for(int x = 0; x<D; x++){
        for(int y = 0; y<D; y++){
            fprintf(fspin,"%d\t%d\t%d\n",x,y,spin_array[D*x +y]);

        }
    }
    fclose(fspin);
}

/*wrapping function for default rand() rng
any rng without inbuilt capability to generate integers in a specific
region must have a similar wrapper written and combined with
appropriate macros.  */
int random_bound_0_1(){
    return rand()% 2;
}

/*Calculates total energy by summing contribution from successive nearest neighbours in the
+ve direction in the x-y plane. Method taken from project brief by CDHW.*/
static double evaluate_total_energy(int* spin_array){
    double total_E = 0.0;/*ensure to include dcp in double defined integer re:feedback.*/

    for(int x = 0; x < D; x++){
        for(int y = 0; y < D; y++){
            total_E -= S(x,y)*(J*(S(x+1,y) + S(x,y+1)) + MU0*B );

        }
    }

    return total_E;

}



/*Plots spin array. Cuurently hardcoded into random initializastion.
Next step is to make this an independent function, callable at any point.*/
static void plot_spin_array(){

    char command[MAXLENGTH];

    snprintf(command, sizeof(command), "%s %s", GNUPLOT, SCRIPT );
    system( command );

}

int main()
{



    int rounds = 5;
    int bound = 2;
    int i = 0;
    VAR_TYPE rng1;
    SRANDOM_R(&rng1, SEED, (intptr_t)&rng1);/*Seeds RND, Make sure RNG is seeded first.*/

    printf("%d\n",rng1);



    int *spin_array = malloc(D*D*sizeof(int));/*Allocates memory for spin array.*/
    initialize_spin_array(spin_array);/*Initializes spin array to random arrangement of +/-1.*/

    if(TESTING){/*Tests if spins are read correctly outside function they are assigned in.*/
        for(int x = 0; x<D;x++){
            for(int y = 0; y<D;y++){
                printf("spin element OUTSIDE %d %d is: %d\n",x,y,spin_array[D*x+y]);
            }
        }
    }

    if(TESTING){/*Tests if arrays are correctly initialized.*/
            for(i = 0; i< rounds; i++){
                printf("bounded rn %d is %ld \n", i, RANDOM_BOUND_0_1(&rng1, BOUND));
            }


        int L = 5;
        //test of the wrapper method. WORKS.
        int y=0;
        for(int x = -1; x< L +1; x++){/*Test of periodic boundary condition element access.*/
                int r = S(x,y);
                printf("x value %d is %d\n",x, r);/*Appears to work.*/


            }

    }
    double result = evaluate_total_energy(spin_array);

    printf("Total energy is %.15lg\n",result);
    spin_array_write_to_file(spin_array);
    plot_spin_array();
    free(spin_array);
    return 0;
}
