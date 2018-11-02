// HEADERS
#include <iostream>
#include <fstream>
#include <string>
#include<time.h>
#include<cstdlib>
#include <math.h>
#include<ctime>
#include <Windows.h>
using namespace std;

//CONSTANTS
//How many times to iterate for the test cases
#define GENERATION_LIMIT 1900
//Number of organisms in Population
#define POPULATION_LIMIT 190
//Propability for crossover to occur
#define CROSSOVER_PROBABILITY 50
//Dependency factor for the mutation equation
#define DEPENDENCY_FACTOR 2
//Lower bound
#define LB -10.0
//Upper bound
#define UB 10.0
//Name of the input file
#define INPUT_FILE_NAME "test_cases/input.txt"

/*Incase we wanted to implement the mutation improvment technique
*Giving each gene a fixed chance to improve fitness value on mutation
*If the value it produced makes the fitness of the organism worse
*If enabled it will take longer to run the algorithm but the results are better
*/
#define ENABLE_MUTATION_IMPROVMENT_TECHNIQUE true
#if ENABLE_MUTATION_IMPROVMENT_TECHNIQUE
#define FIXED_IMPROVMENT_MUTATION_CHANCE 6
//Name of output file
#define OUTPUT_FILE_NAME "result_with_mutation_enhancment.txt"
#endif
#if not  ENABLE_MUTATION_IMPROVMENT_TECHNIQUE
#define OUTPUT_FILE_NAME "result_without_mutation_enhancment.txt"
#endif


//Struct for the points of the polynomial equation
struct Coordinates
{   //X coordinate
	float x_coordinate;
	//Y coordinate
	float y_coordinate;

};


