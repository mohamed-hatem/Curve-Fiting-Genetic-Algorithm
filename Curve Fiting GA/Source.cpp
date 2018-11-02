#include "Population.cpp"

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
		 //Genetic Algorithm	
		for (int j = 0; j<GENERATION_LIMIT; j++)
		{
			if (j == GENERATION_LIMIT / 2)
				cout << "***Half way done***"<<endl;
			current_test_case_population.selection();
			current_test_case_population.sortPopulation(true);
			current_test_case_population.crossover();
			current_test_case_population.mutation(j);
			current_test_case_population.sortPopulation(false);
			current_test_case_population.replacment();
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
	output_file.open("Result.txt");
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
	input_file.open(FILE_NAME);
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
