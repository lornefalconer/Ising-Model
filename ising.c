#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include "pcg_basic.h"

#define ITERATIONS 1000000
#define MAXLENGTH 256
#define D 100/*100 seems the upper limit to reasonable processing time for serial computing.*/
#define UP 1
#define DOWN -1
#define S(x,y) spin_array[D*((D+(x))%D) + ((D+(y))%D)]
#define N 5000

#define MU0 12.57e-7
#define B 0.01
#define K_B 1.38064852e-23
#define T 500.0
#define J (K_B*T)


#define MEM_ERROR 3


#define GNUPLOT "gnuplot"
#define SCRIPT_SPIN "spin_array_gnuplot.script"
#define SCRIPT_EN "total_E_gnuplot.script"

/*RNG selection macros taken from mat_gen.c by CDHW.*/
#define FAIL 0.0
#define NO 0
#define YES 1
#define USE_PCG YES
#define USE_RAND YES
#define SEED_TIME YES
#define SPIN_BOUND 2
#define TESTING NO

#if ( USE_PCG ) /* Bad random number generator */
#  define VAR_TYPE pcg32_random_t
#  define RANDOM_MAX 2 /* Maximum value from random(3) man page */
#  define RANDOM_BOUND(A,B) pcg32_boundedrand_r(A,B)
#  define SRANDOM_R(X,Y,Z) pcg32_srandom_r(X,Y,Z)
#  define RANDOM_DOUBLE_0_1(P,Q) ldexp(P,Q)
#else if( USE_RAND )/* Better random number generator */
#  define VAR_TYPE double
#  define RANDOM_MAX RAND_MAX
#  define RANDOM_BOUND(X,Y) random_bound(Y)
#  define SRANDOM_R(X,Y,Z) srand(Y)
#  define RANDOM_DOUBLE_0_1(P,Q) random_double_0_1()
#endif

#if ( SEED_TIME )
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



static void plot(char* prog_name, char* script){

    char command[MAXLENGTH];

    snprintf(command, sizeof(command), "%s %s", prog_name, script );
    system( command );

}


/*Program to model time evolution of the Ising Model of Ferromagnetism





*/

/*C library random number generator function to generate double [0,1).*/
static double random_double_0_1(){
    return (double)rand()/(double)RANDOM_MAX;
}


/*function to create and initialise a random 2D array of spins*/
void initialize_spin_array(int* spin_array){
    /*map 2D array onto 1D array as its faster*/
    /*elements accessed as element_ij = spin_array(D*i+j)*/
    VAR_TYPE rng1;
    SRANDOM_R(&rng1, SEED, (intptr_t)&rng1);




    for(int x = 0; x<D; x++){
        for(int y = 0; y<D; y++){
            int result = RANDOM_BOUND(&rng1, SPIN_BOUND);

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

static void write_total_E_to_file(double* total_energy_array){
    FILE *fen = fopen("total_energy.txt", "w+");/*opens file to write the ith body data to*/
    err_check_file(fen);
    fprintf(fen,"#E\tIteration\n");

    for(int i = 0; i<ITERATIONS; i++){
            fprintf(fen,"%d\t%lg\n",i,total_energy_array[i]);
        }
    fclose(fen);
}




/*wrapping function for default rand() rng
any rng without inbuilt capability to generate integers in a specific
region must have a similar wrapper written and combined with
appropriate macros.  */
int random_bound(int bound){
    return rand()% bound;
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

/*Calculates the energy required to flip the spin.*/
static double calculate_delta_E(int* spin_array,int x,int y){

    double delta_E = S(x,y)*( 2.0*J*(S(x-1,y)+S(x,y-1)+S(x+1,y)+S(x,y+1)) + MU0*B );

    return delta_E;
}




/*Function to write current delta_E of the ith iteration to an array using the modulo operator.
once all array elements have been filled an average is calculated and stored in an appropriate variable.
the array is then overwritten from the beginning. When full of new values an average is calculated again and compared to
the previous average. If they differ, the system is not in equilibrium and the old average is overwritten by the new.
If they are the same (within a suitable uncertainty window oweing to floating point arithmetic) the return value of the function is
modified to break the outer loop.*/

static double store_consecutive_total_energy_calc_average(double* array_of_consecutive_total_energy, double total_energy, int iteration){

    int modulo_element = (N+iteration)%N;

    //printf("%d\n",modulo_element);

    array_of_consecutive_total_energy[modulo_element] = total_energy;

    if(N-1 == modulo_element){/*If the array has filled up*/
       double total_energy_sum = 0.0;/*Initializes total_energy to floating point 0.0*/
       for(int i = 0; i<N;i++){
            total_energy_sum += array_of_consecutive_total_energy[i];
       }
       double total_energy_av = total_energy_sum/N;/*store average value*/
       //printf("test total energy ac %g\n",total_energy_av);
       return total_energy_av;
    }



    return FAIL;
}





/*Function to evolve system with randomly selected spins.
Ideally would have had a single "flip random spin" function executed in a loop in main,
however I had issues with passing pointers to seeds of the pcg rng.
It was losing sequence information with each iteration and starting from the beginning of the stream each time.
This method produces varying rns.*/
static void evolve_system_random_spin(int* spin_array){




    VAR_TYPE rng_x_coord, rng_y_coord, rng_Metropolis_U;
    SRANDOM_R(&rng_x_coord, SEED, (intptr_t)&rng_x_coord);/*Seeds RNG stream for x_coord of spin.*/
    SRANDOM_R(&rng_y_coord, SEED, (intptr_t)&rng_y_coord);/*Seeds RNG stream for y coord of spin.*/
    SRANDOM_R(&rng_Metropolis_U, SEED, (intptr_t)&rng_Metropolis_U);/*Seeds RNG for Metropolis rn.*/


    double *total_energy_array = malloc(ITERATIONS*sizeof(double));
    double *array_of_consecutive_total_E = malloc(N*sizeof(double));/*Array to store consecutive calculated total_E values.*/

    double average_total_energy_1 = 0.0;
    double average_total_energy_2 = NAN;
    double relative_error = 100.0;
    int count = 0;
    for(int i = 0; i<ITERATIONS; i++){
        //printf("%d\n",i);
        double total_energy = evaluate_total_energy(spin_array);
        total_energy_array[i] = total_energy;

        double result = store_consecutive_total_energy_calc_average(array_of_consecutive_total_E,total_energy,i);
        if(result != 0){
                printf("result %g\n",result);
        }
        if((result != FAIL) && (count == 0)){
            average_total_energy_1 = result;
            count = 1;
        }
        else if(result!= FAIL){
            average_total_energy_2 = result;
            count = 0;
            //printf("average_1 is %g average_2 is %g\n", average_total_energy_1,average_total_energy_2);

        }

        if(fabs((average_total_energy_1 - average_total_energy_2)/average_total_energy_2) < 0.00001){
            break;
        }

        int spin_coord_x = RANDOM_BOUND(&rng_x_coord, (D+1));/*Select spin coords from uniform distribution 0-D*/
        int spin_coord_y = RANDOM_BOUND(&rng_y_coord, (D+1));

        //printf("Spin coord x: %d \nSpin coord y: %d\n",spin_coord_x,spin_coord_y);

        double delta_E = calculate_delta_E(spin_array,spin_coord_x,spin_coord_y);/*Call to func to calculate delta_E.*/




        /*The below method for determining whether to flip spin taken from CDHW and [1].*/
        if(delta_E <= 0.0){/*If delta_E < 0, flipped spin is the more stable state so flips.*/
            S(spin_coord_x,spin_coord_y) = -S(spin_coord_x,spin_coord_y);
        }
        else{
            /*Generate random number U between 0.0 and 1.0 re:CDHWs brief.*/
            /*Input here is for use with PCG RNG. If default rand() selected arguments are
            ignored. Use of PCG RNG to generate doubles in the range [0,1) taken from [2].*/
            double U = RANDOM_DOUBLE_0_1(pcg32_random_r(&rng_Metropolis_U), -32);

            //printf("U is: %.15lg\n",U);
            if(U < exp(-delta_E/(K_B * T))){
                S(spin_coord_x,spin_coord_y) = -S(spin_coord_x,spin_coord_y);
            }
        }


    }

    spin_array_write_to_file(spin_array);
    write_total_E_to_file(total_energy_array);
    free(total_energy_array);
    plot(GNUPLOT,SCRIPT_EN);
}



/*write a function to calculate an updating average of total_E
If the average is unchanging for a given number of iterations the system has reached equilibrium
and we can terminate the simulation.*/





/*Plots spin array. Currently hardcoded into random initializastion.
Next step is to make this an independent function, callable at any point.*/

int main()
{



    int rounds = 5;
    int bound = 2;
    int i = 0;
    VAR_TYPE rng_initialize_spin_array;
    SRANDOM_R(&rng_initialize_spin_array, SEED, (intptr_t)&rng_initialize_spin_array);/*Seeds RND, Make sure RNG is seeded first.*/

    printf("%d\n",rng_initialize_spin_array);



    int *spin_array = malloc(D*D*sizeof(int));/*Allocates memory for spin array.*/
    initialize_spin_array(spin_array);/*Initializes spin array to random arrangement of +/-1.*/
    spin_array_write_to_file(spin_array);
    plot(GNUPLOT,SCRIPT_SPIN);

    if(TESTING){/*Tests if spins are read correctly outside function they are assigned in.*/
        for(int x = 0; x<D;x++){
            for(int y = 0; y<D;y++){
                printf("spin element OUTSIDE %d %d is: %d\n",x,y,spin_array[D*x+y]);
            }
        }
    }

    if(TESTING){/*Tests if arrays are correctly initialized.*/
            for(i = 0; i< rounds; i++){
                printf("bounded rn %d is %ld \n", i, RANDOM_BOUND(&rng_initialize_spin_array, SPIN_BOUND));
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
/*
    for(int x = 0; x < D ; x++){
        for(int y = 0; y < D; y++){
            double delta_E = calculate_delta_E(spin_array,x,y);
            printf("Flipping Energy is: %.15lg\n",delta_E);
        }
    }
*/

    evolve_system_random_spin(spin_array);
    //spin_array_write_to_file(spin_array);
    //plot_spin_array();



    free(spin_array);
    plot(GNUPLOT,SCRIPT_SPIN);


    return 0;
}
