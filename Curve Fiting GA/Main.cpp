//Include Header file
#include "Headers.h"

//Global Variables

//Number of test cases to enter to the program
int number_of_test_cases;
//Array for the number of point for each test case
int * number_of_points;
//Array for the degree of polynomial for each test case
int *degree_of_polynomial;
//2D array for the points of each test case
Coordinates ** list_of_coordinates;
//current test case
int current_test_case_number;

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
		float total_sum = 0.0, current_sum=0.0;
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

//Class to hold the Population of the GA and apply selection , crossover , mutation and replacment on 
class Population
{   //Array for the Initial Population and New Population at the end of each generation
	Organism my_population[POPULATION_LIMIT];

	//Array for the offspring that are generated during crossover
	Organism offspring_population[POPULATION_LIMIT];
public:
	//Get the best organism in the population
	Organism getBest()
	{
		return my_population[0];
	}
	//Free memory used during the end of each test case to be used for the next test case
	void freeResourcesForNextTestCase()
	{
		for (int i = 0; i<POPULATION_LIMIT; i++)
			delete my_population[i].chromosome;
		for (int i = 0; i<POPULATION_LIMIT; i++)
			delete offspring_population[i].chromosome;



	}
	//Biased crossover technique -> if the organism has a good fitness value then increase chance of crossover 
	int evaluateChanceOfCrossover(int index)
	{
		if (index <= 10)
			return 20;
		else if (index <= 20)
			return 15;
		else if (index <= 30)
			return 10;
		else if (index <= 40)
			return 5;
		else return 0;

	}
	//Function to print the population
	void print()
	{
		cout << "start of population" << endl << endl;
		for (int i = 0; i<POPULATION_LIMIT; i++)
		{
			for (int j = 0; j<degree_of_polynomial[current_test_case_number] + 1; j++)
			{
				cout << my_population[i].chromosome[j] << " ";
			}
			cout << endl << my_population[i].fitness << endl;
		}
		cout << "end of population" << endl << endl;
	}
	//Bubble Sort Ascending (min fitness value is best)
	void sortPopulation(bool flag)
	{
		int i, j;
		for (i = 0; i < POPULATION_LIMIT; i++)
		{


			// Last i elements are already in place
			for (j = 0; j < POPULATION_LIMIT - i - 1; j++)
				if (my_population[j].fitness > my_population[j + 1].fitness)
				{
					swapOrganisms(my_population[j], my_population[j + 1]);
				}


		}
		if (flag)
			return;
		for (i = 0; i < POPULATION_LIMIT; i++)
		{


			// Last i elements are already in place
			for (j = 0; j < POPULATION_LIMIT - i - 1; j++)
				if (offspring_population[j].fitness > offspring_population[j + 1].fitness)
				{
					swapOrganisms(offspring_population[j], offspring_population[j + 1]);
				}


		}

	}
	//Function to make deep copy used in sortPopulation()
	void swapOrganisms(Organism &first, Organism &second)
	{
		Organism temp;
		temp.fitness = first.fitness;
		temp.chromosome = new float[degree_of_polynomial[current_test_case_number] + 1];
		for (int i = 0; i<degree_of_polynomial[current_test_case_number] + 1; i++)
		{
			temp.chromosome[i] = first.chromosome[i];
		}
		first.fitness = second.fitness;
		for (int i = 0; i<degree_of_polynomial[current_test_case_number] + 1; i++)
		{
			first.chromosome[i] = second.chromosome[i];
		}
		second.fitness = temp.fitness;
		for (int i = 0; i<degree_of_polynomial[current_test_case_number] + 1; i++)
		{
			second.chromosome[i] = temp.chromosome[i];
		}
		delete temp.chromosome;
	}
	//Selection Function (Tournment)
	void selection()
	{   //Temp array to store current  organisms before selection 
		Organism selection_population[POPULATION_LIMIT];
		//Initialize the selection_population array
		for (int i = 0; i<POPULATION_LIMIT; i++)
			selection_population[i].chromosome = new float[degree_of_polynomial[current_test_case_number] + 1];
		//First and second index of organisms to be compared in tournment selection chosen at random
		int first_index, second_index;
		//Copy current population to selection_population array
		for (int i = 0; i<POPULATION_LIMIT; i++)
			selection_population[i].copyOrganism(my_population[i]);
		
		//Selection process
		//Best 2 organisms are always taken
		for (int i = 2; i<POPULATION_LIMIT; i++)
		{//choose 2 random index from 0 till POPULATION_LIMIT - 1
			first_index = rand() % POPULATION_LIMIT;
			second_index = rand() % POPULATION_LIMIT;
			//choose the organism with the better fitness (Lower)
			if (selection_population[first_index].fitness >= selection_population[second_index].fitness)
			{
				my_population[i].copyOrganism(selection_population[second_index]);
			}
			else
			{
				my_population[i].copyOrganism(selection_population[first_index]);

			}
		}
		//free memory used 
		for (int i = 0; i<POPULATION_LIMIT; i++)
			delete selection_population[i].chromosome;

	}
	//Initialize the population with random floats ranging from -10 till 10
	void initPopulation()
	{   //both variables will be used in generating a random float for the gene
		float random_gene;
		int random_number;
		//Propability that the gene will be negative
		int probability_of_negative_number;

		for (int i = 0; i<POPULATION_LIMIT; i++)
		{   //Initialize the population
			my_population[i].chromosome = new float[degree_of_polynomial[current_test_case_number] + 1];
			for (int j = 0; j<degree_of_polynomial[current_test_case_number] + 1; j++)
			{
				//Choose a random number from 0 till 10
				random_number = rand() % 11;
				/*To generate a float divide by a random number from 1-1000
				*Notice we start from 1 to avoid dividing by zero
				*/
				random_gene = random_number / ((rand() % 1001) + 1.0);
				//Get a random number to be used in determining negative or not
				probability_of_negative_number = rand() % 100;
				//If the random number generated is less than or equal 50 then the gene will be negative
				if (probability_of_negative_number <= 50)
				{
					random_gene *= -1.0;

				}
				//Add the randomly generated float to the chromosome
				my_population[i].chromosome[j] = random_gene;


			}
			//Calculate the fitness value of the organism
			my_population[i].calculateFitness();
		}
	}
	//Crossover Function (randomly generated single point)
	void  crossover()
	{   //Initialize the offspring array
		for (int i = 0; i<POPULATION_LIMIT; i++)
			offspring_population[i].chromosome = new float[degree_of_polynomial[current_test_case_number] + 1];
		/*First and second index will be chosen next to each other
		*Crossover_counter to get exactly same number of population from offspring as current population
		*Crossover_point randomly generate the point where the 2 parents will cross at to produce new offspringss
		*/
		int first_index, second_index, crossover_counter = 0, crossover_point;
		int i = 0;
		//Stop when we have offsprings equal to POPULATION_LIMIT
		while (crossover_counter<POPULATION_LIMIT)
		{   //get 2 indices next to each other to preform crossover on

			first_index = i % POPULATION_LIMIT;
			second_index = (i + 1) % POPULATION_LIMIT;
			//Evaluate the chance of crossover (the better the fitness the better the chance of crossover)  
			if (rand() % 100<CROSSOVER_PROBABILITY + evaluateChanceOfCrossover(first_index) + evaluateChanceOfCrossover(second_index))
			{
				//Randomly generate the crossover point
				crossover_point = rand() % (degree_of_polynomial[current_test_case_number] + 1);
				//Incase the crossover_point is = 0
				if (crossover_point == 0)
				{
					crossover_point++;

				}
				//Produce new offspring based on previous calculations
				produceOffspring(first_index, second_index, crossover_point, crossover_counter);
				//Increase crossover_counter by 2 as 2 new offsprings were produced
				crossover_counter += 2;

			}
			i+=2;
		}


	}
	//Generate 2 new offsprings based on values calculated from Crossover() function
	void produceOffspring(int index_1, int index_2, int crossover_point, int crossover_counter)
	{   //First offspring
		//From 0 till the crossover_point -1 (first parent)
		for (int i = 0; i<crossover_point; i++)
		{
			offspring_population[crossover_counter].chromosome[i] = my_population[index_1].chromosome[i];
		}
		//From crossover_point till the end of the chromosome (second parent)
		for (int i = crossover_point; i<degree_of_polynomial[current_test_case_number] + 1; i++)
		{
			offspring_population[crossover_counter].chromosome[i] = my_population[index_2].chromosome[i];
		}
		//Second offsrpings
		//From 0 till the crossover_point -1 (second parent)
		for (int i = 0; i<crossover_point; i++)
		{
			offspring_population[crossover_counter + 1].chromosome[i] = my_population[index_2].chromosome[i];
		}
		//From crossover_point till the end of the chromosome (first parent)
		for (int i = crossover_point; i<degree_of_polynomial[current_test_case_number] + 1; i++)
		{
			offspring_population[crossover_counter + 1].chromosome[i] = my_population[index_1].chromosome[i];
		}
		//Caluclate fitness value for the new offspringss
		offspring_population[crossover_counter].calculateFitness();
		offspring_population[crossover_counter + 1].calculateFitness();


	}
	//Mutatuion function (non uniform)
	void mutation(int current_itr)
	{
		/*
		Change_effect checks which range of mutation occured  delta_change = my_population[i].chromosome[j] - LB (true) or
		delta_change = delta_change = UB - my_population[i].chromosome[j] (false)
		increase or decrease value of gene during mutation
		*/
		bool  change_effect;
		/*New_gene_value calculated from the equation given y=(1-randomnumber^(currentgeneration/GENERATION_LIMIT)^DEPENDENCY_FACTOR)
		*Equation_variable will be used for the equation
		*/
		float new_gene_value, delta_change, equation_variable;

		float random_number;
		//Storing old values incase we implement Mutation improvment technique
		//Mutation_chance_counter counts how many times mutation produced a bad gene 
#if ENABLE_MUTATION_IMPROVMENT_TECHNIQUE
		float old_fitness_value, old_gene_value;
		int mutation_chance_counter;
#endif

		for (int i = 0; i < POPULATION_LIMIT; i++)
		{   
#if ENABLE_MUTATION_IMPROVMENT_TECHNIQUE
			
			//Set or Reset the mutation_chance_counter
			mutation_chance_counter = 0;
#endif

			for (int j = 0; j < degree_of_polynomial[current_test_case_number] + 1; j++)
			{
				
#if ENABLE_MUTATION_IMPROVMENT_TECHNIQUE
				//Save the old fitness value
				old_fitness_value = my_population[i].fitness;
				//Save the old gene value
				if (mutation_chance_counter < FIXED_IMPROVMENT_MUTATION_CHANCE)
				{
					old_gene_value = my_population[i].chromosome[j];
				}
				else
				{
					mutation_chance_counter = 0;
					continue;
				}

#endif

				//We will Randomly generate a number that will tell if we are gonna increase or decrease the value of the gene during mutation
				if (rand() % 100 <= 50)
				{

					delta_change = my_population[i].chromosome[j] - LB;
					change_effect = true;
				}
				else
				{
					delta_change = UB - my_population[i].chromosome[j];
					change_effect = false;
				}
				//Equation : y=(1-randomnumber^(currentgeneration/GENERATION_LIMIT)^DEPENDENCY_FACTOR)
				equation_variable = 1.0 - ((float)current_itr / GENERATION_LIMIT);
				equation_variable = pow(equation_variable, DEPENDENCY_FACTOR);
				/*Range of values are float 0 -> 1
				*This line of code has imporved the results if we inc the number of zeroes
				*/
				random_number = (rand() % 10000) / 10000.0;
				equation_variable = 1.0 - (pow(random_number, equation_variable));
				equation_variable *= delta_change;

				if (change_effect)
					my_population[i].chromosome[j] = my_population[i].chromosome[j] - equation_variable;


				else
					my_population[i].chromosome[j] = my_population[i].chromosome[j] + equation_variable;

#if ENABLE_MUTATION_IMPROVMENT_TECHNIQUE
				my_population[i].calculateFitness();

				if (my_population[i].fitness > old_fitness_value)
				{
					my_population[i].chromosome[j] = old_gene_value;
					my_population[i].fitness = old_fitness_value;
					mutation_chance_counter++;
					j--;
				}

				else if (my_population[i].fitness < old_fitness_value)
				{
					mutation_chance_counter = 0;
					continue;
				}



#endif
			}
#if (not ENABLE_MUTATION_IMPROVMENT_TECHNIQUE)
			my_population[i].calculateFitness();
#endif


		}
		//Applying the same algorithm for the offspring array
		for (int i = 0; i<POPULATION_LIMIT; i++)
		{
#if ENABLE_MUTATION_IMPROVMENT_TECHNIQUE
			
			//Set or Reset the mutation_chance_counter
			mutation_chance_counter = 0;
#endif


			for (int j = 0; j<degree_of_polynomial[current_test_case_number] + 1; j++)
			{
#if ENABLE_MUTATION_IMPROVMENT_TECHNIQUE
				old_fitness_value = offspring_population[i].fitness;
				if (mutation_chance_counter < FIXED_IMPROVMENT_MUTATION_CHANCE)
				{
					old_gene_value = offspring_population[i].chromosome[j];
				}
				else
				{
					mutation_chance_counter = 0;
					continue;
				}
#endif


				if (rand() % 100 <= 50)
				{

					delta_change = offspring_population[i].chromosome[j] - LB;
					change_effect = true;
				}
				else
				{
					delta_change = UB - offspring_population[i].chromosome[j];
					change_effect = false;
				}

				equation_variable = 1.0 - ((float)current_itr / GENERATION_LIMIT);
				equation_variable = pow(equation_variable, DEPENDENCY_FACTOR);
				random_number = (rand() % 10000) / 10000.0;
				equation_variable = 1.0 - (pow(random_number, equation_variable));
				equation_variable *= delta_change;
				if (change_effect)
					offspring_population[i].chromosome[j] = offspring_population[i].chromosome[j] - equation_variable;


				else
					offspring_population[i].chromosome[j] = offspring_population[i].chromosome[j] + equation_variable;

#if ENABLE_MUTATION_IMPROVMENT_TECHNIQUE
				offspring_population[i].calculateFitness();
				if (offspring_population[i].fitness > old_fitness_value)
				{
					offspring_population[i].chromosome[j] = old_gene_value;
					offspring_population[i].fitness = old_fitness_value;
					mutation_chance_counter++;
					j--;
				}

				else if (offspring_population[i].fitness <= old_fitness_value)
				{
					mutation_chance_counter = 0;
					continue;
				}
#endif
			}
#if (not ENABLE_MUTATION_IMPROVMENT_TECHNIQUE)
			offspring_population[i].calculateFitness();
#endif



		}




	}
	//Replacment function (Elitism)
	void replacment()
	{
		Organism replacment_population[POPULATION_LIMIT];
		for (int i = 0; i<POPULATION_LIMIT; i++)
			replacment_population[i].chromosome = new float[degree_of_polynomial[current_test_case_number] + 1];
		int first_index, second_index;
		for (int i = 0; i<POPULATION_LIMIT; i++)
			replacment_population[i].copyOrganism(my_population[i]);
		int new_generation_counter = 0, my_population_counter = 0, replacment_population_counter = 0;
		while (my_population_counter<POPULATION_LIMIT)
		{
			if (replacment_population[replacment_population_counter].fitness<offspring_population[new_generation_counter].fitness)
			{
				my_population[my_population_counter++].copyOrganism(replacment_population[replacment_population_counter++]);
			}
			else
				my_population[my_population_counter++].copyOrganism(offspring_population[new_generation_counter++]);
		}
		for (int i = 0; i < POPULATION_LIMIT; i++)
		{
			delete replacment_population[i].chromosome;
			
		}
	}
	//Dealocate offspring dynamic data for next iteration of genetic algorithm 
	void freeOffspring()
	{ 
		
		for (int i = 0; i < POPULATION_LIMIT; i++)
		{
			
			
			delete offspring_population[i].chromosome;
			
		}
	}

};


//Function Prototypes
void dealocateResources(Organism * &solutions);
void getInputFromFile();
void outputToFile(Organism * &solutions);
int main()
{   //To calculate how many seconds the genetic algorithm took for each test case 
	clock_t begin, end;
	//Array of the best solutions in our program
	Organism * solutions;
	//Seeding the rand function
	srand(time(0));
	//Change color of console text = black bgcolor = lightgrey
	system("color 70");
	//Get input from file
	getInputFromFile();
	//Float to calculate percentage done by the algorithm
	float percentage_of_current_test_case_done;
	//Population for the test cases
	Population current_test_case_population;
	//Initalize the solutions array with the number of test cases
	solutions = new Organism[number_of_test_cases];
	for (int i = 0; i<number_of_test_cases; i++)
		solutions[i].chromosome = new float[degree_of_polynomial[i] + 1];
	for (int i = 0; i<number_of_test_cases; i++)
	{   //Initialize and sort the current population
		cout << "***Starting test case "<<i+1 <<"***"<< endl;
		current_test_case_number = i;
		current_test_case_population.initPopulation();
		current_test_case_population.sortPopulation(true);
		begin = clock();
		percentage_of_current_test_case_done = 10.0;
		 //Genetic Algorithm	
		for (int j = 0; j<GENERATION_LIMIT; j++)
		{
		    //Indication for percentage done by the algorithm	
			if (percentage_of_current_test_case_done <= (float)((j + 1.0) / GENERATION_LIMIT)*100.0)
			{
				cout << "***" << percentage_of_current_test_case_done << "% completed***" << endl;
				percentage_of_current_test_case_done += 10.0;
			}
			current_test_case_population.selection();
			current_test_case_population.sortPopulation(true);
			current_test_case_population.crossover();
			current_test_case_population.mutation(j);
			current_test_case_population.sortPopulation(false);
			current_test_case_population.replacment();
			if (j < GENERATION_LIMIT-1)
			{
				current_test_case_population.freeOffspring();
			}
		}
		//Save the best organism for each test case at the end of the genetic algorithm
		solutions[i].copyOrganism(current_test_case_population.getBest());
		//Dealocate memory to utilize memory usage
		current_test_case_population.freeResourcesForNextTestCase();
		end = clock();
		cout << "***Test case " << i + 1 << " done.  time taken for this test : ";
		cout<< double(end - begin) / CLOCKS_PER_SEC << " seconds***"<< endl;
	}
	//Output solutions for each test case to the file named "Result.txt"
	outputToFile(solutions);
	//Free all used dynamic memory to prevent memory leakage
	dealocateResources(solutions);
	system("pause");
	return 0;
}
//Function to output results to file
void outputToFile(Organism * &solutions)
{
	ofstream output_file;
	output_file.open(OUTPUT_FILE_NAME);
	for (int j = 0; j<number_of_test_cases; j++)
	{
		output_file << "test case " << j + 1 << " : " << endl;
		output_file << "coefficients: ";
		for (int i = 0; i<degree_of_polynomial[j] + 1; i++)
		{   
			
			output_file << solutions[j].chromosome[i];
			if (i == 1)
				output_file << "x";
			else if (i > 1)
				output_file << "x^" << i;
			if (i != degree_of_polynomial[j])
			{
				output_file << " + ";
			}
		}
		output_file << endl << "error = " << solutions[j].fitness << endl;
	}




}
//Function to  Read Input from file
void getInputFromFile()
{
	ifstream input_file;
	input_file.open(INPUT_FILE_NAME);
	input_file >> number_of_test_cases;
	number_of_points = new int[number_of_test_cases];
	degree_of_polynomial = new int[number_of_test_cases];
	list_of_coordinates = new Coordinates*[number_of_test_cases];
	for (int i = 0; i<number_of_test_cases; i++)
	{
		input_file >> number_of_points[i];
		input_file >> degree_of_polynomial[i];

		list_of_coordinates[i] = new Coordinates[number_of_points[i]];
		for (int j = 0; j<number_of_points[i]; j++)
		{

			input_file >> list_of_coordinates[i][j].x_coordinate;
			input_file >> list_of_coordinates[i][j].y_coordinate;

		}

	}

	input_file.close();

}
//Function to prevent memory leakage after program finishes
void dealocateResources(Organism * &solutions)
{
	delete number_of_points;
	delete degree_of_polynomial;
	for (int i = 0; i<number_of_test_cases; i++)
		delete list_of_coordinates[i];
	delete[]list_of_coordinates;
	for (int i = 0; i<number_of_test_cases; i++)
	{
		delete solutions[i].chromosome;

	}
	delete solutions;

}
