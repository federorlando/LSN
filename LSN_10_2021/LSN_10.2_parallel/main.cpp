/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include "random.h"
#include "functions.h"
#include "TSP_position.h"
#include "mpi.h"

using namespace std;



int main (int argc, char *argv[]){

//////////////////////////////////MPI stuff

int size, rank;
MPI_Init(&argc,&argv);
MPI_Comm_size(MPI_COMM_WORLD, &size);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Status stat1, stat2;
MPI_Request req;

if(size%2 != 0)
{
	cout << "Please run on even number of nodes." << endl << endl;
	return 1;
}

vector <int> processors;	//used to pick randomly pairs of nodes to exchange info
for(int i=0; i<size; i++)
        processors.push_back(0);

int node_a, node_b;
int best_element_index;
double walker;
int itag1, itag2;

//////////////////////////////////

Random rnd;
int seed[4];
int p1, p2;

ifstream Primes("Primes");

if (Primes.is_open())
{
	for(int i=0; i<rank; i++)
		Primes.ignore(10,'\n');		//each process must read a different couple here to be truly independent

	Primes >> p1 >> p2;
}

else cerr << "PROBLEM: Unable to open Primes" << endl;
Primes.close();
    
ifstream input("seed.in");
string property;
if (input.is_open())
{
        while ( !input.eof() )
	{
		input >> property;
		if( property == "RANDOMSEED" )
		{
	                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
	                rnd.SetRandom(seed,p1,p2);
		}
        }
        input.close();
}
else cerr << "PROBLEM: Unable to open seed.in" << endl;

////////////////PARAMETERS

int N = 50;		//Number of cities
int N_pop = 46;		//Number of chromosomes in a population
int steps = pow(10,4);	//Number of generations
int N_migr = 500;	//Number of generation between migrations j-th and (j+1)-th

vector <int> order;
for(int i=0; i<N; i++)
        order.push_back(i+1);

TSP_position pos[N];		//VECTOR OF N OBJECTS

vector <double> appo;
for(int i=0; i<N; i++)
        appo.push_back(0);





/////////////////////////////////////////////////////////    
/////////////////////////////////////////////////////////
//						       //
//                 AMERICAN CAPITALS                   //
//						       //
/////////////////////////////////////////////////////////
///////////////////////////////////////////////////////// 

double Latitude[N];
double Longitude[N];
string lat, lon;

ifstream load_capitals;
load_capitals.open("American_capitals.dat");

if (load_capitals.is_open())
{
	for(int i=0; i<N; i++)
	{
		load_capitals >> lat >> lon;
		Latitude[i] = stod(lat);
		Longitude[i] = stod(lon);
	}		
}

else cerr << "PROBLEM: Unable to open American_capitals.dat" << endl;

load_capitals.close();

for(int i=0; i<N; i++)
{
        pos[i].SetX(Latitude[i]);
        pos[i].SetY(Longitude[i]);
	//cout << "Lat " << Latitude[i] << "              long " << Longitude[i] << endl;
}
    
ofstream square_pos;
square_pos.open("path_rank" + to_string(rank));

for(int i=0; i<N; i++)
        square_pos << pos[i].GetX() << " " << pos[i].GetY() << endl;

//cout << "Check: " << TSP_check(order) << endl;
//cout << "Cost: " << TSP_cost(order, pos) << endl;
 
//////////////Randomly create first generation 
vector <vector <int>> pop;		//each vector is a chromosome

for(int i=0; i<N_pop; i++)		//careful: N_pop numer of chromosomes, N number of bits (cities)
{
        for(int j=0; j<N; j++)
		random_shuffle(order.begin()+1,order.end());	//+1 to left 1st city untouched

        pop.push_back(order);
}

//////////////Start procedure

vector <vector <int>> new_pop;
double Rand = 0, randnew = 0, randnewnew = 0, index1 = 0, index2 = 0;
    
vector <double> fitness;
for(int i=0; i<N_pop; i++)
        fitness.push_back(0.);
    
vector <double> average_cost;  		
for(int i=0; i<steps; i++)
        average_cost.push_back(0.);

////////////////////////////////////////////////////
//cout << "Initial vector: " << endl;
for(int j=0; j<size; j++)
{
   	processors[j] = j;				//vector is = [0,1,2....size] (ordered)
	cout << processors[j] << "  ";
}
cout << endl;

////////////////////////////////////////////////////
    
for(int i=0; i<steps; i++)
{       
        //Population fitness of the i-th generation
        for(int j=0; j<N_pop; j++)
		fitness[j] = TSP_cost(pop[j],pos);     
        
        //Creation of the new (i+1-th) generation
        while(new_pop.size() < N_pop)
	{
            
		//Select 1st & 2nd chromosome (pair of parent chr.) from old (i-th) generation
		index1 = selection(fitness,rnd.Rannyu());
		index2 = selection(fitness,rnd.Rannyu());
    
		//Accumulate them - neatly - in new (i+1-th) population       
            	new_pop.push_back(pop[index1]);
            	new_pop.push_back(pop[index2]);
            
		//Crossover: Pc > 0.5
		if(rnd.Rannyu()<0.7)
			crossover(new_pop[new_pop.size()-1], new_pop[new_pop.size()-2], (int)round(rnd.Rannyu(1,N)));
            
		//Mutation of the second chromosome
		Rand = rnd.Rannyu();
		if(Rand<0.05) 
			swap(new_pop[new_pop.size()-1],(int)round(rnd.Rannyu(1,N-1)),(int)round(rnd.Rannyu(1,N-1)));		//1, and not 0, to leave 1st city untouched
		if(Rand>0.25 && Rand<0.3)
			shift(new_pop[new_pop.size()-1],(int)round(rnd.Rannyu(1,N)));						//not touching 1st city implicit in function def
		if(Rand>0.5 && Rand<0.55)
			invert_order(new_pop[new_pop.size()-1],(int)round(rnd.Rannyu(1,N-1)),(int)round(rnd.Rannyu(1,N-1)));	//1, and not 0, to leave 1st city untouched
		if(Rand>0.75 && Rand<0.8)
			part_perm(new_pop[new_pop.size()-1], (int)round(rnd.Rannyu(1,9)), (int)round(rnd.Rannyu(10,17)), (int)round(rnd.Rannyu(18,26)), (int)round(rnd.Rannyu(27,34)));													//1, and not 0, to leave 1st city untouched
            
		//Mutation of the first chromosome
		Rand = rnd.Rannyu();
		if(Rand<0.05)
			swap(new_pop[new_pop.size()-2],(int)round(rnd.Rannyu(1,N-1)),(int)round(rnd.Rannyu(1,N-1)));
		if(Rand>0.25 && Rand<0.3)
			shift(new_pop[new_pop.size()-2],(int)round(rnd.Rannyu(1,N)));
		if(Rand>0.5 && Rand<0.55)
			invert_order(new_pop[new_pop.size()-2],(int)round(rnd.Rannyu(1,N-1)),(int)round(rnd.Rannyu(1,N-1)));
		if(Rand>0.75 && Rand<0.8)
			part_perm(new_pop[new_pop.size()-2], (int)round(rnd.Rannyu(1,9)), (int)round(rnd.Rannyu(10,17)), (int)round(rnd.Rannyu(18,26)), (int)round(rnd.Rannyu(27,34)));													//1, and not 0, to leave 1st city untouched
        }
           
        //Population average cost (based on the N_pop/2 best paths) of the i-th generation
        sort(fitness.begin(),fitness.end());

        for(int k=0; k<(N_pop/2); k++)
            average_cost[i] += fitness[k];
        
        average_cost[i]/=(N_pop/2);
        
        //New population becomes current one
        for(int k=0; k<N_pop; k++)
		pop[k] = new_pop[k];

        new_pop.clear();

	/////////////////////////////////////////////////////////
	//						       //
	//                      MIGRATION                      //
	//						       //
	/////////////////////////////////////////////////////////

	if(i%N_migr==0)
	//if(i==110)
	{
		cout << "i = " << i << "... Preparing for migration!" << endl;
		//random_shuffle(processors.begin(),processors.end()); //NO! only works if all nodes read the same in Primes!
		
		//Shuffling vector of nodes
		srand(5);
		randnew = rand();
		randnew = randnew/RAND_MAX;
		randnewnew = rand();
		randnewnew = randnewnew*(size-1)/RAND_MAX;
		//cout << "RANDNEW " << randnew << endl;
		if(randnew<0.3) 
			swap(processors,(int)round(randnewnew),(int)round(randnewnew));		
		if(randnew>0.35 && randnew<0.65)
			shift(processors,(int)round(randnewnew));						
		if(randnew>0.7 && randnew<1)
			invert_order(processors,(int)round(randnewnew),(int)round(randnewnew));	

		//cout << "Post shuffle: " << endl;
		//for(int j=0; j<size; j++)				
			//cout << processors[j] << "  ";
		//cout << endl;	

		for(int k=0; k<size; k=k+2)				
		{
			//choose pair of nodes 
			node_a = processors[k];
			node_b = processors[k+1];
			//cout << k << "           node_a " << node_a << "		node_b " << node_b << endl;
			
			//... for both nodes, get best-fitting chromosome out of sampled population 
			best_element_index = 0;
			walker = TSP_cost(pop[0],pos);

			for(int l=1; l<N_pop; l++)
			{
        			if(TSP_cost(pop[l],pos) < walker)
				best_element_index = l;
			}

			//copy chromosome in appo
			for(int t=0; t<N; t++)
				appo[t] = pop[best_element_index][t];
	
			//cout << "Before migration. Rank " << rank << "	best el index " << best_element_index << "	and component 33 of best el " << appo[33] << endl;

			//MIGRATE APPO'S
			itag1 = 1;
			itag2 = 2;
  
			if (rank == node_a)
			{
				MPI_Isend(appo.data(), N, MPI_DOUBLE, node_b, itag1, MPI_COMM_WORLD, &req);
				MPI_Recv(appo.data(), N, MPI_DOUBLE, node_b, itag2, MPI_COMM_WORLD, &stat2);
			}

			else if (rank == node_b)
			{
				MPI_Send(appo.data(), N, MPI_DOUBLE, node_a, itag2, MPI_COMM_WORLD);
				MPI_Recv(appo.data(), N, MPI_DOUBLE, node_a, itag1, MPI_COMM_WORLD, &stat1);
			}

			//cout << "After migr. a->b. Rank " << rank << "	best el index " << best_element_index << "	and component 33 of best el " << appo[33] << endl;

			//copy appo's back into chromosomes
			for(int t=0; t<N; t++)
				pop[best_element_index][t] = appo[t];

		}
	}

	/////////////////////////////////////////////////////////
	//						       //
	//                   MIGRATION ENDED                   //
	//						       //
	/////////////////////////////////////////////////////////
}
    
//////////////Find chromosome in last population which gives shortest path

int best_path_index = 0;
double min_fit = TSP_cost(pop[0],pos);

for(int i=1; i<N_pop; i++)
{
        if(TSP_cost(pop[i],pos) < min_fit)
		best_path_index = i;
}
    
square_pos <<  endl;
for(auto el: pop[best_path_index])
        square_pos << el << endl;

square_pos <<  endl;

double final_cost = TSP_cost(pop[best_path_index],pos);
cout << endl << "rank " << rank << "	FINAL COST = " << final_cost << endl << endl;

square_pos << "	FINAL COST = " << final_cost << endl;
square_pos.close();
  
//Data blocking
vector <double> sum_prog;
vector <double> err_prog;
    
for(int l=0; l<steps; l++)
{
        sum_prog.push_back(0);
        err_prog.push_back(0);
}
    
data_blocking(average_cost,sum_prog,err_prog);
    
ofstream square_output;
square_output.open("cost_rank"+to_string(rank)); 
for(int i=0; i<steps; i++)
	square_output << average_cost[i] << " " << sum_prog[i] << " " << err_prog[i] << endl;   
square_output.close();

rnd.SaveSeed();

MPI_Finalize();

return 0;
}



/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
