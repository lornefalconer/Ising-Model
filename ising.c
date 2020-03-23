#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include "pcg_basic.h"

#define CHUNKSIZE 20

#define MU0 12.57e-7
#define B 0.0
#define K_B 1.38064852e-23

#define J 1.0e-22
#define T_0 (J/K_B)

#define ITERATIONS 100000
#define TEMP_START 1.0
#define TEMP_MAX (5*T_0)
#define TEMP_INTERVAL (TEMP_MAX/500)
#define MAXLENGTH 256
#define D 20/*100 seems the upper limit to reasonable processing time for serial computing.*/
#define UP 1
#define DOWN -1
#define S(x,y) spin_array[D*((D+(x))%D) + ((D+(y))%D)]
#define STDEV_ARRAY_LENGTH 1000

#define INVESTIGATE_C_VS_D YES
#define CUTOFF_TESTING YES

#define MEM_ERROR 3


#define GNUPLOT "gnuplot"
#define SCRIPT_SPIN "spin_array_gnuplot.script"
#define SCRIPT_EN "total_E_gnuplot.script"
#define SCRIPT_C "specific_heat_with_temp_plot.script"
#define SCRIPT_M "gnuplot_magnetism.script"
#define SCRIPT_S "gnuplot_susceptability.script"

/*RNG selection macros taken from mat_gen.c by CDHW.*/
#define FAIL 0.0
#define NO 0
#define YES 1
#define USE_PCG YES
#define USE_RAND YES
#define SEED_TIME YES
#define SPIN_BOUND 2
#define TESTING NO
#define PLOTTING YES
#define PLOTTING_DEBUG NO

#define RANDOM_SPIN YES

#define SERIAL YES
#define PARALLEL NO

#if ( USE_PCG ) /* Bad random number generator */
#  define VAR_TYPE pcg32_random_t
#  define RANDOM_MAX 2 /* Maximum value from random(3) man page */
#  define RANDOM_BOUND(A,B) pcg32_boundedrand_r(A,B)
#  define SRANDOM_R(X,Y,Z) pcg32_srandom_r(X,Y,Z)
#  define RANDOM_DOUBLE_0_1(P,Q) ldexp(P,Q)
#else/* Better random number generator */
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

static void write_specific_heat_to_file(double* specific_heat_array, double* temp_array){
    FILE *fc = fopen("specific_heat.txt", "w+");/*opens file to write the ith body data to*/
    err_check_file(fc);
    fprintf(fc,"#C\tTemp\n");

    for(int i = 0; i<(TEMP_MAX/TEMP_INTERVAL); i++){
            fprintf(fc,"%lg\t%lg\n",T_0*fabs(temp_array[i]/T_0 - 2.27),specific_heat_array[i]/(D*D*K_B));
        }
    fclose(fc);
}

static void write_susceptability_to_file(double* susceptability, double* temp_array){
    FILE *fc = fopen("susceptability.txt", "w+");/*opens file to write the ith body data to*/
    err_check_file(fc);
    fprintf(fc,"#Chi\tTemp\n");

    for(int i = 0; i<(TEMP_MAX/TEMP_INTERVAL); i++){
            fprintf(fc,"%lg\t%lg\n",T_0*fabs(temp_array[i]/T_0 - 2.27),susceptability[i]/(D*D*K_B));
        }
    fclose(fc);
}

static void write_magnetism_to_file(double* magnetism_array, double* temp_array){
    FILE *fc = fopen("magnetism.txt", "w+");/*opens file to write the ith body data to*/
    err_check_file(fc);
    fprintf(fc,"#M\tTemp\n");

    for(int i = 0; i<(TEMP_MAX/TEMP_INTERVAL); i++){
            fprintf(fc,"%lg\t%lg\n",T_0*fabs(temp_array[i]/T_0 - 2.27),magnetism_array[i]/(MU0*D*D));
        }
    fclose(fc);
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

static double store_consecutive_quantity_calc_average(double* array_of_consecutive_quantity, double quantity, int iteration){

    int modulo_element = (STDEV_ARRAY_LENGTH+iteration)%STDEV_ARRAY_LENGTH;

    //printf("%d\n",modulo_element);

    array_of_consecutive_quantity[modulo_element] = quantity;

    if(STDEV_ARRAY_LENGTH-1 == modulo_element){/*If the array has filled up*/
       double total_quantity_sum = 0.0;/*Initializes total_energy to floating point 0.0*/
       for(int i = 0; i<STDEV_ARRAY_LENGTH;i++){
            total_quantity_sum += array_of_consecutive_quantity[i];
       }
       double total_quantity_av = total_quantity_sum/STDEV_ARRAY_LENGTH;/*store average value*/
       //printf("test total energy ac %g\n",total_energy_av);
       return total_quantity_av;
    }



    return FAIL;/*If array has not filled up, return a value signifying failure.*/
}




static double calculate_magnetisation(int* spin_array){
    int total_spin = 0;
    for(int i = 0; i< D*D; i++){
        total_spin += spin_array[i];

    }

    //printf("TEST: TOTAL SPIN: %d\n",total_spin);
    double Magnetisation = fabs((double)total_spin/(D*D));
    //printf("MAGE TEST: %lf\n",Magnetisation);


    return Magnetisation;
}




/*Function to evolve system with randomly selected spins.
Ideally would have had a single "flip random spin" function executed in a loop in main,
however I had issues with passing pointers to seeds of the pcg rng.
It was losing sequence information with each iteration and starting from the beginning of the stream each time.
This method produces varying rns.*/


static double calculate_specific_heat(double* energy_squared_array, double* array_of_consecutive_total_E, double average_total_energy,double T){
    double total_energy_squared = 0.0;
    for(int j = 0; j<STDEV_ARRAY_LENGTH;j++){
                    energy_squared_array[j] = array_of_consecutive_total_E[j]*array_of_consecutive_total_E[j];
                    total_energy_squared += energy_squared_array[j];

            }
    double average_total_energy_squared = total_energy_squared/STDEV_ARRAY_LENGTH;
    double specific_heat_cap = (average_total_energy_squared - pow(average_total_energy,2))/(K_B*pow(T,2));
    printf("specific heat: %.25e\n",specific_heat_cap);
    return specific_heat_cap;

}

static double calculate_susceptability(double* magnetism_squared_array, double* array_of_consecutive_magnetism, double average_magnetism,double T){
    double total_magnetism_squared = 0.0;
    for(int j = 0; j<STDEV_ARRAY_LENGTH;j++){
                    magnetism_squared_array[j] = array_of_consecutive_magnetism[j]*array_of_consecutive_magnetism[j];
                    total_magnetism_squared += magnetism_squared_array[j];

            }
    double average_magnetism_squared = total_magnetism_squared/STDEV_ARRAY_LENGTH;
    double susceptability = (average_magnetism_squared - pow(average_magnetism,2))/(K_B*T);
    printf("susceptability: %.25e\n",susceptability);
    return susceptability;

}



///Need to fiddle with this so it re-initializes the spin array before the next temperature iteration///

static double* evolve_system_random_spin_beta(int* spin_array, VAR_TYPE rng_x_coord, VAR_TYPE rng_y_coord, VAR_TYPE rng_Metropolis_U, double Temp){


    double *total_energy_array = malloc(ITERATIONS*sizeof(double));
    double *array_of_consecutive_total_E = malloc(STDEV_ARRAY_LENGTH*sizeof(double));/*Array to store consecutive calculated total_E values.*/
    double *energy_squared_array = malloc(STDEV_ARRAY_LENGTH*sizeof(double));
    double *magnetism_squared_array = malloc(STDEV_ARRAY_LENGTH*sizeof(double));
    double *array_of_consecutive_magnetism = malloc(STDEV_ARRAY_LENGTH*sizeof(double));

    double average_total_energy_1 = 0.0;
    double average_total_energy_2 = NAN;
    double average_total_magnetism = 0.0;
    double relative_error = 100.0;
    int count = 0;

    int spin_coord_x = 0;
    int spin_coord_y = 0;

        for(int i = 0; i<ITERATIONS; i++){
            //printf("MAGNETSATION(IN): %f\n", calculate_magnetisation(spin_array));
            //printf("%d\n",i);
            double total_energy = evaluate_total_energy(spin_array);/*evaluate total energy of spin array*/
            total_energy_array[i] = total_energy;/*Store the energy in an array for later plotting and eqm analysis.*/




            double magnetism = calculate_magnetisation(spin_array);

            if(CUTOFF_TESTING){
                double result_E = store_consecutive_quantity_calc_average(array_of_consecutive_total_E,total_energy,i);
                double result_M = store_consecutive_quantity_calc_average(array_of_consecutive_magnetism,magnetism,i);
                if(result_E != 0){
                    //printf("result %g\n",result);
                }
                if((result_E != FAIL) && (count == 0)){
                    average_total_energy_1 = result_E;
                    count = 1;
                }
                else if(result_E!= FAIL){
                    average_total_energy_2 = result_E;
                    average_total_magnetism = result_M;
                    count = 0;
                //printf("average_1 is %g average_2 is %g\n", average_total_energy_1,average_total_energy_2);

                }

                if(fabs((average_total_energy_1 - average_total_energy_2)/average_total_energy_2) < 0.0001){

                    break;
                }
            }

            if( RANDOM_SPIN ){
                spin_coord_x = RANDOM_BOUND(&rng_x_coord, (D+1));/*Select spin coords from uniform distribution 0-D*/
                spin_coord_y = RANDOM_BOUND(&rng_y_coord, (D+1));

                //printf("Spin coord x: %d \nSpin coord y: %d\n",spin_coord_x,spin_coord_y);
            }
            else{
                spin_coord_x = i % D;

                if(spin_coord_x == (D-1)){
                    spin_coord_y += 1;

                }


            }
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
                    if(U < exp(-delta_E/(K_B * Temp))){
                        S(spin_coord_x,spin_coord_y) = -S(spin_coord_x,spin_coord_y);
                    }
                }



        }

    double specific_heat = calculate_specific_heat(energy_squared_array, array_of_consecutive_total_E, average_total_energy_2,Temp);
    double susceptability = calculate_susceptability(magnetism_squared_array,array_of_consecutive_magnetism,average_total_magnetism,Temp);
    double magnetism = calculate_magnetisation(spin_array);

    double *ret_vals = (double*)malloc(3*sizeof(double));
    ret_vals[0] = specific_heat;
    ret_vals[1] = magnetism;
    ret_vals[2] = susceptability;
    //printf("MAGNETSATION: %.15g\n", calculate_magnetisation(spin_array));

    spin_array_write_to_file(spin_array);
    write_total_E_to_file(total_energy_array);

    free(total_energy_array);
    //free(spin_array);
    if(PLOTTING_DEBUG){
        plot(GNUPLOT,SCRIPT_EN);
        //plot(GNUPLOT,SCRIPT_C);
    }


    return ret_vals;
}


/*write a function to calculate an updating average of total_E
If the average is unchanging for a given number of iterations the system has reached equilibrium
and we can terminate the simulation.*/


/**Use the program to investigate how eqm magnetisation depends on T.


Magnetisation is the sum of spins over the number of spins.



**/



static double repeat_with_iterating_temp(){

    VAR_TYPE rng_initialize_spin_array;/*declare variable to initialize pcg rng*/
    SRANDOM_R(&rng_initialize_spin_array, SEED, (intptr_t)&rng_initialize_spin_array);/*Seeds RND, Make sure RNG is seeded first.*/

    printf("%d\n",rng_initialize_spin_array);



    int *spin_array = malloc(D*D*sizeof(int));/*Allocates memory for spin array.*/

    initialize_spin_array(spin_array);/*Initializes spin array to random arrangement of +/-1.*/




    spin_array_write_to_file(spin_array);/*writes the initial spin configuration to a file for debug plotting.*/

    //if( PLOTTING_DEBUG ){
      //  plot(GNUPLOT,SCRIPT_SPIN);
    //}



    printf("MAGNETSATION: %.25e\n", calculate_magnetisation(spin_array));

/*
    for(int x = 0; x < D ; x++){
        for(int y = 0; y < D; y++){
            double delta_E = calculate_delta_E(spin_array,x,y);
            printf("Flipping Energy is: %.15lg\n",delta_E);
        }
    }
*/

    VAR_TYPE rng_x_coord, rng_y_coord, rng_Metropolis_U;/*decalre variables to initialize 3 other rngs*/
    SRANDOM_R(&rng_x_coord, SEED, (intptr_t)&rng_x_coord);/*Seeds RNG stream for x_coord of spin.*/
    SRANDOM_R(&rng_y_coord, SEED, (intptr_t)&rng_y_coord);/*Seeds RNG stream for y coord of spin.*/
    SRANDOM_R(&rng_Metropolis_U, SEED, (intptr_t)&rng_Metropolis_U);/*Seeds RNG for Metropolis rn.*/


    double Temp = TEMP_START;/*Declare starting temperature for interation.*/

    /*allocate memory for plotting of C,M and Chi*/
    double* specific_heat_array = malloc((TEMP_MAX/TEMP_INTERVAL)*sizeof(double));
    double* magnetism_array = malloc((TEMP_MAX/TEMP_INTERVAL)*sizeof(double));
    double* susceptability_array = malloc((TEMP_MAX/TEMP_INTERVAL)*sizeof(double));

    /*Allocate memory for plotting of temperature againt above quantities*/
    double* temp_array = malloc((TEMP_MAX/TEMP_INTERVAL)*sizeof(double));
    /*Count variable will ensure allocation of memory to correct indices in memory*/
    int count = 0;

    double max_val = 0.0;

    while(Temp < TEMP_MAX){/*While the temp is less tan some upper limit*/
        printf("\n\nTEMP:%f\n\n",Temp);

        initialize_spin_array(spin_array);/*Initialize spin array to random spins.*/
                if(PLOTTING_DEBUG && Temp == 1.0){
                    spin_array_write_to_file(spin_array);/*Write to file in case you want to debug the plot*/
                }
        /*Evolvoe the system with time and the current Temperature*/

        double *ret_vals = evolve_system_random_spin_beta(spin_array, rng_x_coord, rng_y_coord, rng_Metropolis_U, Temp);
        double specific_heat = ret_vals[0];

        if(INVESTIGATE_C_VS_D){/*If you are investigating the relationship between C and D*/
            if(Temp == TEMP_START){
                max_val = specific_heat;

            }
            else if(specific_heat > max_val){
                    max_val = specific_heat;

            }
        }

        double magnetism = ret_vals[1];
        double susceptability = ret_vals[2];
        free(ret_vals);


                if(PLOTTING_DEBUG){
                        spin_array_write_to_file(spin_array);
                        plot(GNUPLOT,SCRIPT_SPIN);
                }

        specific_heat_array[count] = specific_heat;
        magnetism_array[count] = magnetism;
        susceptability_array[count] = susceptability;
        temp_array[count] = Temp;


        Temp += TEMP_INTERVAL;
        count += 1;

        printf("Specific heat: %.25e\n", specific_heat);
        printf("MAGNETSATION: %f\n", magnetism);

    }
    printf("MAX C %.15e\n",max_val);
    write_specific_heat_to_file(specific_heat_array, temp_array);
    write_magnetism_to_file(magnetism_array,temp_array);
    write_susceptability_to_file(susceptability_array,temp_array);

    free(susceptability_array);
    free(specific_heat_array);
    free(magnetism_array);
    free(temp_array);
    free(spin_array);


    if(PLOTTING){
        //plot(GNUPLOT,SCRIPT_SPIN);
        plot(GNUPLOT,SCRIPT_C);
        plot(GNUPLOT,SCRIPT_M);
        plot(GNUPLOT,SCRIPT_S);

    }
}


/*Plots spin array. Currently hardcoded into random initializastion.
Next step is to make this an independent function, callable at any point.*/

int main()
{

if(SERIAL){/**EXECUTES IN ABOUT 22s FOR 100000 ITERATIONS**/

    int rounds = 5;
    int bound = 2;
    int i = 0;
    VAR_TYPE rng_initialize_spin_array;/*declare variable to initialize pcg rng*/
    SRANDOM_R(&rng_initialize_spin_array, SEED, (intptr_t)&rng_initialize_spin_array);/*Seeds RND, Make sure RNG is seeded first.*/

    printf("%d\n",rng_initialize_spin_array);



    int *spin_array = malloc(D*D*sizeof(int));/*Allocates memory for spin array.*/

    initialize_spin_array(spin_array);/*Initializes spin array to random arrangement of +/-1.*/




    spin_array_write_to_file(spin_array);/*writes the initial spin configuration to a file for debug plotting.*/

    //if( PLOTTING_DEBUG ){
      //  plot(GNUPLOT,SCRIPT_SPIN);
    //}



    printf("MAGNETSATION: %.25e\n", calculate_magnetisation(spin_array));

/*
    for(int x = 0; x < D ; x++){
        for(int y = 0; y < D; y++){
            double delta_E = calculate_delta_E(spin_array,x,y);
            printf("Flipping Energy is: %.15lg\n",delta_E);
        }
    }
*/

    VAR_TYPE rng_x_coord, rng_y_coord, rng_Metropolis_U;/*decalre variables to initialize 3 other rngs*/
    SRANDOM_R(&rng_x_coord, SEED, (intptr_t)&rng_x_coord);/*Seeds RNG stream for x_coord of spin.*/
    SRANDOM_R(&rng_y_coord, SEED, (intptr_t)&rng_y_coord);/*Seeds RNG stream for y coord of spin.*/
    SRANDOM_R(&rng_Metropolis_U, SEED, (intptr_t)&rng_Metropolis_U);/*Seeds RNG for Metropolis rn.*/


    double Temp = TEMP_START;/*Declare starting temperature for interation.*/

    /*allocate memory for plotting of C,M and Chi*/
    double* specific_heat_array = malloc((TEMP_MAX/TEMP_INTERVAL)*sizeof(double));
    double* magnetism_array = malloc((TEMP_MAX/TEMP_INTERVAL)*sizeof(double));
    double* susceptability_array = malloc((TEMP_MAX/TEMP_INTERVAL)*sizeof(double));

    /*Allocate memory for plotting of temperature againt above quantities*/
    double* temp_array = malloc((TEMP_MAX/TEMP_INTERVAL)*sizeof(double));
    /*Count variable will ensure allocation of memory to correct indices in memory*/
    int count = 0;

    double max_val = 0.0;

    while(Temp < TEMP_MAX){/*While the temp is less tan some upper limit*/
        printf("\n\nTEMP:%f\n\n",Temp);

        initialize_spin_array(spin_array);/*Initialize spin array to random spins.*/
                if(PLOTTING_DEBUG && Temp == 1.0){
                    spin_array_write_to_file(spin_array);/*Write to file in case you want to debug the plot*/
                }
        /*Evolvoe the system with time and the current Temperature*/

        double *ret_vals = evolve_system_random_spin_beta(spin_array, rng_x_coord, rng_y_coord, rng_Metropolis_U, Temp);
        double specific_heat = ret_vals[0];

        if(INVESTIGATE_C_VS_D){/*If you are investigating the relationship between C and D*/
            if(Temp == TEMP_START){
                max_val = specific_heat;

            }
            else if(specific_heat > max_val){
                    max_val = specific_heat;

            }
        }

        double magnetism = ret_vals[1];
        double susceptability = ret_vals[2];
        free(ret_vals);


                if(PLOTTING_DEBUG){
                        spin_array_write_to_file(spin_array);
                        plot(GNUPLOT,SCRIPT_SPIN);
                }

        specific_heat_array[count] = specific_heat;
        magnetism_array[count] = magnetism;
        susceptability_array[count] = susceptability;
        temp_array[count] = Temp;


        Temp += TEMP_INTERVAL;
        count += 1;

        printf("Specific heat: %.25e\n", specific_heat);
        printf("MAGNETSATION: %f\n", magnetism);

    }
    printf("MAX C %.15e\n",max_val);
    write_specific_heat_to_file(specific_heat_array, temp_array);
    write_magnetism_to_file(magnetism_array,temp_array);
    write_susceptability_to_file(susceptability_array,temp_array);

    free(susceptability_array);
    free(specific_heat_array);
    free(magnetism_array);
    free(temp_array);
    free(spin_array);


    if(PLOTTING){
        //plot(GNUPLOT,SCRIPT_SPIN);
        plot(GNUPLOT,SCRIPT_C);
        plot(GNUPLOT,SCRIPT_M);
        plot(GNUPLOT,SCRIPT_S);

    }

}










else if(PARALLEL){/**EXECUTES IN 14s WITH 100000 ITERATIONS AND STATIC ALLOCATION**/

    printf("Hello (Parallel) world\n");
    /*Initialise vairables and seed rng outside of parallel regions*/
    int rounds = 5;
    int bound = 2;


    double *total_energy_array = malloc(ITERATIONS*sizeof(double));
    double *array_of_consecutive_total_E = malloc(STDEV_ARRAY_LENGTH*sizeof(double));/*Array to store consecutive calculated total_E values.*/


    VAR_TYPE rng_initialize_spin_array;
    SRANDOM_R(&rng_initialize_spin_array, SEED, (intptr_t)&rng_initialize_spin_array);/*Seeds RND, Make sure RNG is seeded first.*/

    printf("%d\n",rng_initialize_spin_array);

    int *spin_array = malloc(D*D*sizeof(int));/*Allocates memory for spin array.*/
    initialize_spin_array(spin_array);/*Initializes spin array to random arrangement of +/-1.*/
    spin_array_write_to_file(spin_array);
    plot(GNUPLOT,SCRIPT_SPIN);

    /*Will need to re-engineer original function to pass coord rng seed to function.*/



    VAR_TYPE rng_x_coord, rng_y_coord,rng_Metropolis_U1,rng_Metropolis_U2,rng_Metropolis_U3,rng_Metropolis_U4;
    SRANDOM_R(&rng_x_coord, SEED, (intptr_t)&rng_x_coord);/*Seeds RNG stream for x_coord of spin.*/
    SRANDOM_R(&rng_y_coord, SEED, (intptr_t)&rng_y_coord);/*Seeds RNG stream for y coord of spin.*/
    SRANDOM_R(&rng_Metropolis_U1, SEED, (intptr_t)&rng_Metropolis_U1);/*Seeds RNG for Metropolis rn.*/
    SRANDOM_R(&rng_Metropolis_U2, SEED, (intptr_t)&rng_Metropolis_U2);/*Seeds RNG for Metropolis rn.*/
    SRANDOM_R(&rng_Metropolis_U3, SEED, (intptr_t)&rng_Metropolis_U3);/*Seeds RNG for Metropolis rn.*/
    SRANDOM_R(&rng_Metropolis_U4, SEED, (intptr_t)&rng_Metropolis_U4);/*Seeds RNG for Metropolis rn.*/

    VAR_TYPE Metropolis_array[4] = {rng_Metropolis_U1,rng_Metropolis_U2,rng_Metropolis_U3,rng_Metropolis_U4};

    //evolve_system_random_spin_beta(spin_array, rng_x_coord, rng_y_coord, rng_Metropolis_U);


    int chunk = CHUNKSIZE;

    int tid, nthreads, strips;
    int i = 0;
    /*declare parallel environment and variables which are shared and private to each thread.*/
    #pragma omp parallel shared(chunk,spin_array,total_energy_array,rng_x_coord,rng_y_coord,nthreads,Metropolis_array) private(tid,i)
    {
        tid = omp_get_thread_num();
        printf("Hello World from thread = %d\n", tid);

        if (tid == 0)
        {
        nthreads = omp_get_num_threads();
        strips = nthreads;
        printf("Number of threads = %d\n", nthreads);
        /*allocates an array to hold delta_E of each strip/core*/
        }
        printf("Outside test %d\n", strips);
        int strip_width = D/strips;
        int *numb_strips_x_coords = malloc(strips*sizeof(int));/*dynamically allocates array of strips, one for each core.*/
        double* delta_E_strip_array = malloc(strips*sizeof(int));

        //VAR_TYPE* Metropolis_seed_array = malloc(strips*sizeof(int));



    /**need to re-engineer this to to have iteration outside the loop**/

    #pragma omp for schedule(static,chunk)

        for (int i=0; i<ITERATIONS; i++)
        {
            double T = 0.0;
            double total_energy = evaluate_total_energy(spin_array);/*evaluate total energy of spin array*/
            total_energy_array[i] = total_energy;/*Store the energy in an array for later plotting and eqm analysis.*/

            //printf("ITERATION: %d\n",i);
            int spin_coord_x_strip = RANDOM_BOUND(&rng_x_coord, (strip_width+1));/*Select x spin coord from uniform distribution 0-D/nthreads*/
            int spin_coord_y = RANDOM_BOUND(&rng_y_coord, (D+1));/*selects y spin coord along whole length of strip.*/
            //printf("TEST: %d\n", i);
            /*dynamically generate x coord for each strip (symmetric to avoid deadlocks)*/
            for(int j = 0; j<nthreads; j++){/*updates array of x coords for each strip.*/
                numb_strips_x_coords[j] = spin_coord_x_strip + i*strip_width;
                delta_E_strip_array[j] = calculate_delta_E(spin_array,numb_strips_x_coords[j],spin_coord_y);

                if(delta_E_strip_array[j] <= 0.0){/*If delta_E < 0, flipped spin is the more stable state so flips.*/
                    S(numb_strips_x_coords[j],spin_coord_y) = -S(numb_strips_x_coords[j],spin_coord_y);
                }
                else{
                    double U = RANDOM_DOUBLE_0_1(pcg32_random_r(&Metropolis_array[j]), -32);

                    if(U < exp(-delta_E_strip_array[j]/(K_B * T))){
                        S(numb_strips_x_coords[j],spin_coord_y) = -S(numb_strips_x_coords[j],spin_coord_y);
                    }
                }

            }


        }
        //evolve_system_random_spin_beta(spin_array, rng_x_coord, rng_y_co\ord, rng_Metropolis_U);



    printf("TEST");
    }

    spin_array_write_to_file(spin_array);
    write_total_E_to_file(total_energy_array);
    free(total_energy_array);
    plot(GNUPLOT,SCRIPT_EN);
    plot(GNUPLOT,SCRIPT_SPIN);


    free(spin_array);
    printf("END");

}



    return 0;
}
