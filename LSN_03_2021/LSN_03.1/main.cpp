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

////////////////Parameters

const double S_0 = 100;
const double T = 1;
const double K = 100;
const double r = 0.1;
const double sigma = 0.25;

////////////////

fstream outfile;

int M=100000;           //Total number of throws
int N=100;              //Number of blocks
int L=M/N;            	//Number of throws in each block, please use for M a multiple of N

vector <double> call;
vector <double> put;
vector <double> ave_call;
vector <double> ave_put;
vector <double> sum_prog;
vector <double> err_prog;




////////////////////////////////////////////////DIRECT SAMPLING

for(int i=0; i<N; i++)
{							
        sum_prog.push_back(0);
        err_prog.push_back(0);
        call.push_back(0);
        put.push_back(0);
}

double S = 0;

for(int i=0; i<N; i++)	
{
	for(int j=0; j<L; j++)
	{
		S = S_0*exp((r-0.5*pow(sigma,2))*T+sigma*rnd.Gauss(0.,T)*sqrt(T));
		call[i] += exp(-r*T)*max(0.,S-K);
		put[i] += exp(-r*T)*max(0.,K-S);
	}
	ave_call.push_back(call[i]/L);
        ave_put.push_back(put[i]/L);
}

////////////////Data blocking & Write to file (CALL)

data_blocking(ave_call, sum_prog, err_prog);

outfile.open("output_1_call.txt",ios::out);

if(!outfile)
{
	cerr << "Error in opening output file" << endl;
	return -1;
}

for(int i=0; i<N; i++)
	outfile << i << "   "  << sum_prog[i] << "   "  << err_prog[i] << endl;

outfile.close();

////////////////Data blocking & Write to file (PUT)

data_blocking(ave_put, sum_prog, err_prog);

outfile.open("output_1_put.txt",ios::out);

if(!outfile)
{
	cerr << "Error in opening output file" << endl;
	return -1;
}

for(int i=0; i<N; i++)
	outfile << i << "   "  << sum_prog[i] << "   "  << err_prog[i] << endl;

outfile.close();



////////////////////////////////////////////////DISCRETE SAMPLING

vector <double> S_ds;

for(int i=0; i<N; i++)
{
	S_ds.push_back(0);							
        call[i] = 0;
        put[i] = 0;
        ave_call[i] = 0;
        ave_put[i] = 0;
}

S_ds[0] = S_0;
const int Nstep = 100;	//number of time intervals in which [0,T] is divided
					
for(int i=0; i<N; i++)	
{
	for(int j=0; j<L; j++)
	{
		for(int k=1; k<Nstep; k++)	//start from k=1 because we have fixed S_ds[0] = S_0
		{
			S_ds[k] = S_ds[k-1]*exp(	(r-0.5*pow(sigma,2))*(T/Nstep)	+ 	sigma*sqrt(T/Nstep)*rnd.Gauss(0,1)	);
		}
		call[i] += exp(-r*T)*max(0.,S_ds[Nstep-1]-K);
		put[i] += exp(-r*T)*max(0.,K-S_ds[Nstep-1]);
	}
	ave_call[i] = call[i]/L;
	ave_put[i] = put[i]/L;
}

////////////////Data blocking & Write to file (CALL)

data_blocking(ave_call, sum_prog, err_prog);

outfile.open("output_2_call.txt",ios::out);

if(!outfile)
{
	cerr << "Error in opening output file" << endl;
	return -1;
}

for(int i=0; i<N; i++)
	outfile << i << "   "  << sum_prog[i] << "   "  << err_prog[i] << endl;

outfile.close();

////////////////Data blocking & Write to file (PUT)

data_blocking(ave_put, sum_prog, err_prog);

outfile.open("output_2_put.txt",ios::out);

if(!outfile)
{
	cerr << "Error in opening output file" << endl;
	return -1;
}

for(int i=0; i<N; i++)
	outfile << i << "   "  << sum_prog[i] << "   "  << err_prog[i] << endl;

outfile.close();




rnd.SaveSeed();

return 0;
}

