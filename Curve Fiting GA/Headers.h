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
#define FILE_NAME "input.txt"
/*Incase we wanted to implement the mutation improvment technique
Giving each gene a fixed chance to improve fitness value on mutation
If the value it produced makes the fitness of the organism worse
*/
#define ENABLE_MUTATION_IMPROVMENT_TECHNIQUE true
#if ENABLE_MUTATION_IMPROVMENT_TECHNIQUE
#define FIXED_IMPROVMENT_MUTATION_CHANCE 6
#endif

//Struct for the points of the polynomial equation
struct Coordinates
{   //X coordinate
	float x_coordinate;
	//Y coordinate
	float y_coordinate;

};

//Struct for the organism of our Genetic Algorithm
struct Organism
{   //Floating point Array for our organism (GenoType) -> array to hold  each gene value
	float * chromosome;
	//Error value of the organism
	float fitness;
	//Function to copy 1 organism into another
	void copyOrganism(Organism x)
	{   //Copy fitness value 
		fitness = x.fitness;
		//Iterate over each gene and copy it
		for (int i = 0; i<degree_of_polynomial[current_test_case_number] + 1; i++)
		{
			chromosome[i] = x.chromosome[i];


		}

	}
	//Function to calculate the value of y coordinate given the cofficients and the x coordinate (ycalc)
	float calculateYCoordinate(float x_coordinate, int number_of_coefficients)
	{   //ex: y = a0+a1*x+a2*x^2......
		float total_sum = 0.0, current_sum;
		for (int i = 0; i<number_of_coefficients; i++)
		{
			current_sum = chromosome[i];
			for (int j = 0; j<i; j++)
			{
				current_sum *= x_coordinate;
			}
			total_sum += current_sum;
		}
		return total_sum;
	}
	//Function to calculate the fitness value of the organism using the given sigma equaion 1/N*(sigma((ycalc-yactual)^2))
	void calculateFitness()
	{
		fitness = 0.0;
		float current_x_coordinate, current_actual_y_coordinate, current_calculated_y_coordinate;
		for (int i = 0; i<number_of_points[current_test_case_number]; i++)
		{   //Get the each x coordinate of the current test case
			current_x_coordinate = list_of_coordinates[current_test_case_number][i].x_coordinate;
			//Get each actual y coordinate of the current test case
			current_actual_y_coordinate = list_of_coordinates[current_test_case_number][i].y_coordinate;
			//Calculate the y coordinate
			current_calculated_y_coordinate = calculateYCoordinate(current_x_coordinate, degree_of_polynomial[current_test_case_number] + 1);
			//Applying the sigma equation for each ycalc and yactual
			fitness = fitness + ((current_calculated_y_coordinate - current_actual_y_coordinate)*(current_calculated_y_coordinate - current_actual_y_coordinate));
		}
		// 1/N * sigma result
		fitness = fitness / (float)number_of_points[current_test_case_number];
	}


};