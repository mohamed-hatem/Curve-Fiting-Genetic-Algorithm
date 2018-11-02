//Include Header file
#include "Headers.h"
//Class to hold the Population of the GA
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
		// first and second index of organisms to be compared in tournment selection chosen at random
		int first_index, second_index;
		//copy current population to selection_population array
		for (int i = 0; i<POPULATION_LIMIT; i++)
			selection_population[i].copyOrganism(my_population[i]);
		//selection process
		for (int i = 0; i<POPULATION_LIMIT; i++)
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
				/*To generate a float divide by a random number from 1-10
				*Nnotice we start from 1 to avoid dividing by zero
				*/
				random_gene = random_number / ((rand() % 11) + 1.0);
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
			second_index = i + 1 % POPULATION_LIMIT;
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
			i++;
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

		double random_number;
		//Storing old values incase we implement Mutation improvment technique
		//Mutation_chance_counter counts how many times mutation produced a bad gene 
#if ENABLE_MUTATION_IMPROVMENT_TECHNIQUE
		float old_fitness_value, old_gene_value;
		int mutation_chance_counter;
#endif

		for (int i = 0; i < POPULATION_LIMIT; i++)
		{   //Save the old fitness value
#if ENABLE_MUTATION_IMPROVMENT_TECHNIQUE
			old_fitness_value = my_population[i].fitness;
			//Set or Reset the mutation_chance_counter
			mutation_chance_counter = 0;
#endif

			for (int j = 0; j < degree_of_polynomial[current_test_case_number] + 1; j++)
			{
				//Save the old gene value
#if ENABLE_MUTATION_IMPROVMENT_TECHNIQUE
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
				random_number = (rand() % 100) / 100.0;
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
			old_fitness_value = offspring_population[i].fitness;
			//Set or Reset the mutation_chance_counter
			mutation_chance_counter = 0;
#endif


			for (int j = 0; j<degree_of_polynomial[current_test_case_number] + 1; j++)
			{
#if ENABLE_MUTATION_IMPROVMENT_TECHNIQUE
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
				random_number = (rand() % 100) / 100.0;
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
		for (int i = 0; i<POPULATION_LIMIT; i++)
			delete replacment_population[i].chromosome;
	}

};