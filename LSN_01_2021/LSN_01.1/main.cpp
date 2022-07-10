#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include "random.h"
#include "functions.h"

using namespace std;

int main (int argc, char *argv[]){

////////////////Prepare object rnd

Random rnd;
int seed[4];
int p1, p2;

ifstream Primes("Primes");
if (Primes.is_open())
	Primes >> p1 >> p2;
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

//////////////////////////////////////////////////////////////PART 1

fstream outfile;

int M=100000;           //Total number of throws
int N=100;              //Number of blocks
int L=M/N;            	//Number of throws in each block, please use for M a multiple of N

vector <double> r;
vector <double> ave;
vector <double> sum_prog;
vector <double> err_prog;

////////////////Fill r with M casual numbers and set all other vectors to 0 

for(int i=0; i<M; i++)			
	r.push_back(rnd.Rannyu());

for(int i=0; i<N; i++)
{							
        sum_prog.push_back(0);
        err_prog.push_back(0);
}

////////////////Compute average for each of the N blocks and apply function data_blocking

double sum=0;
int k=0;

for(int i=0; i<N; i++)	
{
	sum = 0;
	for(int j=0; j<L; j++)
	{
		k = i*L+j;
		sum += r[k];
	}
        ave.push_back(sum/L);
}

data_blocking(ave, sum_prog, err_prog);

////////////////Write to file

outfile.open("output_1.txt",ios::out);

if(!outfile)
{
	cerr << "Error in opening output file" << endl;
	return -1;
}

for(int i=0; i<N; i++)
	outfile << i << "   "  << sum_prog[i]-0.5 << "   "  << err_prog[i] << endl;

outfile.close();






//////////////////////////////////////////////////////////////PART 2

////////////////Compute average for each of the N blocks and apply function data_blocking

k=0;

for(int i=0; i<N; i++)	
{
	sum = 0;
	for(int j=0; j<L; j++)
	{
		k = i*L+j;
		sum += pow(r[k]-0.5,2);
	}
	ave[i] = sum/(1.*L);
}

data_blocking (ave, sum_prog, err_prog);

////////////////Write to file

outfile.open("output_2.txt",ios::out);

if(!outfile)
{
	cerr << "Error in opening output file" << endl;
	return -1;
}

for(int i=0; i<N; i++)
	outfile << i << "   "  << sum_prog[i]-1./12 << "   "  << err_prog[i] << endl;

outfile.close();






//////////////////////////////////////////////////////////////PART 3

int m = 100;
vector <int> n;     				//vector counting numbers in [0,1/m]
vector <double> chi;

for(int k=0; k<m; k++)
{
	n.push_back(0);
	chi.push_back(0);       		//since N=m=100
}
 
double rand = 0;

for(int i=0; i<N; i++)
{
        for(int l=0; l<m; l++)
		n[l]=0;

        for(int j=0; j<10000; j++)
	{
            	rand = rnd.Rannyu();
            	for(int k=0; k<m; k++)
                	if( rand>=(double)k/m & rand<(double)(k+1)/m ) n[k]++;
        }

        for(int s=0; s<m; s++)
		chi[i]=chi[i]+pow(n[s]-10000/m,2)/(10000/m); 
}
   
data_blocking(chi, sum_prog, err_prog);

outfile.open("output_3.txt",ios::out);

if(!outfile)
{
	cerr << "Error in opening output file" << endl;
	return -1;
}

    
for(int k=0; k<N; k++)
	outfile << k << " " << sum_prog[k] << " " << err_prog[k] << endl;

outfile.close();

rnd.SaveSeed();

return 0;
}

