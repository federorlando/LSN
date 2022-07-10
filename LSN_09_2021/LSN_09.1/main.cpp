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

using namespace std;



int main (int argc, char *argv[]){

Random rnd;
int seed[4];
int p1, p2;

ifstream Primes("Primes");
if (Primes.is_open())
{
	Primes.ignore(10,'\n');
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

int N = 34;		//Number of cities
int N_pop = 46;		//Number of chromosomes in a population
int steps = pow(10,4);	//Number of generations

vector <int> order;
for(int i=0; i<N; i++)
        order.push_back(i+1);

TSP_position pos[N];		//VECTOR OF N OBJECTS





/////////////////////////////////////////////////////////    
/////////////////////////////////////////////////////////
//						       //
//                    UNIT CIRCLE                      //
//						       //
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////

double angle=0;

for(int i=0; i<N; i++)
{
        angle = 2*M_PI*rnd.Rannyu();
        pos[i].SetX(cos(angle));
        pos[i].SetY(sin(angle));
}
    
ofstream circ_pos;
circ_pos.open("circle_path.txt"); 
for(int i=0; i<N; i++)
        circ_pos << pos[i].GetX() << " " << pos[i].GetY() << endl;

//cout << "Check: " << TSP_check(order) << endl;
//cout << "Cost: " << TSP_cost(order, pos) << endl;
 
//Randomly create first generation 
vector <vector <int>> pop;		//each vector is a chromosome

for(int i=0; i<N_pop; i++)		//careful: N_pop numer of chromosomes, N number of bits (cities)
{
        for(int j=0; j<N; j++)
		random_shuffle(order.begin()+1,order.end());	//+1 to left 1st city untouched

        pop.push_back(order);
}
    
cout << "UNIT CIRCLE. First population: " << endl << endl;

for(int i=0; i<pop.size(); i++)		//de facto pop.size() = N_pop
{
        cout << "Chromosome "<< i+1 <<": ";

        for(auto el: pop[i])
		cout << el << " ";

        cout<<endl;
}

cout << endl << endl;

//Start procedure

vector <vector <int>> new_pop;
double rand = 0, index1 = 0, index2 = 0;
    
vector <double> fitness;
for(int i=0; i<N_pop; i++)
        fitness.push_back(0.);
    
vector <double> average_cost;  		
for(int i=0; i<steps; i++)
        average_cost.push_back(0.);
    
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
		rand = rnd.Rannyu();
		if(rand<0.05) 
			swap(new_pop[new_pop.size()-1],(int)round(rnd.Rannyu(1,N-1)),(int)round(rnd.Rannyu(1,N-1)));		//1, and not 0, to leave 1st city untouched
		if(rand>0.25 && rand<0.3)
			shift(new_pop[new_pop.size()-1],(int)round(rnd.Rannyu(1,N)));						//not touching 1st city implicit in function def
		if(rand>0.5 && rand<0.55)
			invert_order(new_pop[new_pop.size()-1],(int)round(rnd.Rannyu(1,N-1)),(int)round(rnd.Rannyu(1,N-1)));	//1, and not 0, to leave 1st city untouched
		if(rand>0.75 && rand<0.8)
			part_perm(new_pop[new_pop.size()-1], (int)round(rnd.Rannyu(1,9)), (int)round(rnd.Rannyu(10,17)), (int)round(rnd.Rannyu(18,26)), (int)round(rnd.Rannyu(27,34)));													//1, and not 0, to leave 1st city untouched
            
		//Mutation of the first chromosome
		rand = rnd.Rannyu();
		if(rand<0.05)
			swap(new_pop[new_pop.size()-2],(int)round(rnd.Rannyu(1,N-1)),(int)round(rnd.Rannyu(1,N-1)));
		if(rand>0.25 && rand<0.3)
			shift(new_pop[new_pop.size()-2],(int)round(rnd.Rannyu(1,N)));
		if(rand>0.5 && rand<0.55)
			invert_order(new_pop[new_pop.size()-2],(int)round(rnd.Rannyu(1,N-1)),(int)round(rnd.Rannyu(1,N-1)));
		if(rand>0.75 && rand<0.8)
			part_perm(new_pop[new_pop.size()-2], (int)round(rnd.Rannyu(1,9)), (int)round(rnd.Rannyu(10,17)), (int)round(rnd.Rannyu(18,26)), (int)round(rnd.Rannyu(27,34)));													//1, and not 0, to leave 1st city untouched
        }
           
        //Population average cost (based on the N_pop/2 best paths) of the i-th generation
        sort(fitness.begin(),fitness.end());

        for(int k=0; k<(N_pop/2); k++)
            average_cost[i] += fitness[k];
        
        average_cost[i]/=(N_pop/2);
        
        //New population becomes current one
        for(int k=0;k<N_pop;k++)
		pop[k]=new_pop[k];

        new_pop.clear();
}
    
//Find chromosome in last population which gives shortest path

cout << "Final population: " << endl << endl;

for(int i=0; i<pop.size(); i++)		
{
        cout << "Chromosome "<< i+1 <<": ";

        for(auto el: pop[i])
		cout << el << " ";

        cout<<endl;
}

int best_path_index = 0;
double min_fit = TSP_cost(pop[0],pos);

for(int i=1; i<N_pop; i++)
{
        if(TSP_cost(pop[i],pos) < min_fit)
		best_path_index = i;
}
    
circ_pos <<  endl;
for(auto el: pop[best_path_index])
        circ_pos << el << endl;

circ_pos.close();

double final_cost = TSP_cost(pop[best_path_index],pos);
cout << endl << "UNIT CIRCLE: FINAL COST = " << final_cost << endl << endl << endl;
  
//Data blocking
vector <double> sum_prog;
vector <double> err_prog;
    
for(int l=0; l<steps; l++)
{
        sum_prog.push_back(0);
        err_prog.push_back(0);
}
    
data_blocking(average_cost,sum_prog,err_prog);
    
ofstream circ_output;
circ_output.open("circle_output.txt"); 
for(int i=0; i<steps; i++)
	circ_output << average_cost[i] << " " << sum_prog[i] << " " << err_prog[i] << endl;   
circ_output.close();



/////////////////////////////////////////////////////////    
/////////////////////////////////////////////////////////
//						       //
//                SQUARE [0,1] x [0,1]                 //
//						       //
/////////////////////////////////////////////////////////
///////////////////////////////////////////////////////// 

for(int i=0;i<N;i++)
{
        rand = rnd.Rannyu();
        pos[i].SetX(rand);
        rand = rnd.Rannyu();
        pos[i].SetY(rand);
}
    
ofstream square_pos;
square_pos.open("square_path.txt");
    
for(int i=0; i<N; i++)
        square_pos << pos[i].GetX() << " " << pos[i].GetY() << endl;
    
for(int i=0; i<N_pop; i++)
{
        for(int j=0; j<N; j++)
		random_shuffle(order.begin()+1,order.end());
        
	pop[i]=order;
}

cout << "SQUARE. First population: " << endl << endl;

for(int i=0; i<pop.size(); i++)		//de facto pop.size() = N_pop
{
        cout << "Chromosome "<< i+1 <<": ";

        for(auto el: pop[i])
		cout << el << " ";

        cout<<endl;
}

cout << endl;
    
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
		rand = rnd.Rannyu();
		if(rand<0.05) 
			swap(new_pop[new_pop.size()-1],(int)round(rnd.Rannyu(1,N-1)),(int)round(rnd.Rannyu(1,N-1)));
		if(rand>0.25 && rand<0.3)
			shift(new_pop[new_pop.size()-1],(int)round(rnd.Rannyu(1,N)));
		if(rand>0.5 && rand<0.55)
			invert_order(new_pop[new_pop.size()-1],(int)round(rnd.Rannyu(1,N-1)),(int)round(rnd.Rannyu(1,N-1)));
		if(rand>0.75 && rand<0.8)
			part_perm(new_pop[new_pop.size()-1], (int)round(rnd.Rannyu(1,9)), (int)round(rnd.Rannyu(10,17)), (int)round(rnd.Rannyu(18,26)), (int)round(rnd.Rannyu(27,34)));	
            
		//Mutation of the first chromosome
		rand = rnd.Rannyu();
		if(rand<0.05)
			swap(new_pop[new_pop.size()-2],(int)round(rnd.Rannyu(1,N-1)),(int)round(rnd.Rannyu(1,N-1)));
		if(rand>0.25 && rand<0.3)
			shift(new_pop[new_pop.size()-2],(int)round(rnd.Rannyu(1,N)));
		if(rand>0.5 && rand<0.55)
			invert_order(new_pop[new_pop.size()-2],(int)round(rnd.Rannyu(1,N-1)),(int)round(rnd.Rannyu(1,N-1)));
		if(rand>0.75 && rand<0.8)
			part_perm(new_pop[new_pop.size()-2], (int)round(rnd.Rannyu(1,9)), (int)round(rnd.Rannyu(10,17)), (int)round(rnd.Rannyu(18,26)), (int)round(rnd.Rannyu(27,34)));	
        }
           
        //Population average cost (based on the N_pop/2 best paths) of the i-th generation
        sort(fitness.begin(),fitness.end());

        for(int k=0; k<(N_pop/2); k++)
            average_cost[i] += fitness[k];
        
        average_cost[i]/=(N_pop/2);
        
        //New population becomes current one
        for(int k=0;k<N_pop;k++)
		pop[k]=new_pop[k];

        new_pop.clear();
}
    
//Find chromosome in last population which gives shortest path 

cout << "Final population: " << endl << endl;

for(int i=0; i<pop.size(); i++)		
{
        cout << "Chromosome "<< i+1 <<": ";

        for(auto el: pop[i])
		cout << el << " ";

        cout<<endl;
}

best_path_index = 0;
min_fit = TSP_cost(pop[0],pos);

for(int i=1; i<N_pop; i++)
{
        if(TSP_cost(pop[i],pos) < min_fit)
		best_path_index = i;
}
    
square_pos <<  endl;
for(auto el: pop[best_path_index])
        square_pos << el << endl;
    
square_pos.close();

final_cost = TSP_cost(pop[best_path_index],pos);
cout << endl << "SQUARE: FINAL COST = " << final_cost << endl << endl;
    
//Data blocking
    
data_blocking(average_cost,sum_prog,err_prog);
    
ofstream square_output;
square_output.open("square_output.txt"); 
for(int i=0; i<steps; i++)
	square_output << average_cost[i] << " " << sum_prog[i] << " " << err_prog[i] << endl;   
square_output.close();
    
rnd.SaveSeed();
    
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
