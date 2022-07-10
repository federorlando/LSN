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
		sum += (M_PI/2)*cos(M_PI*r[k]/2);
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
	outfile << i << "   "  << sum_prog[i]-1 << "   "  << err_prog[i] << endl;

outfile.close();






//////////////////////////////////////////////////////////////PART 2 

vector <double> s;

for(int i=0; i<M; i++)			
	s.push_back(	(1./100)*(101-sqrt(10201-10200*(r[i])))	   );
	//s.push_back(-2/M_PI*log(1-(1-exp(-M_PI/2))*r[i]));

////////////////Compute average for each of the N blocks and apply function data_blocking

k=0;

for(int i=0; i<N; i++)	
{
	sum = 0;
	for(int j=0; j<L; j++)
	{
		k = i*L+j;
		sum += ((M_PI/2)*cos(M_PI*s[k]/2))/((100*1./51)*((101*1./100)-s[k]));
		//sum=sum+(M_PI/2)*cos(M_PI*s[k]/2)/(M_PI/2*exp(-M_PI/2*s[k])/(1-exp(-M_PI/2)) );
	}
        ave[i] = (sum/L);			//NON USARE PUSHBACK QUI!!!!!!!!!!!!!!!!!!!!!!!
}

data_blocking(ave, sum_prog, err_prog);

////////////////Write to file

outfile.open("output_2.txt",ios::out);

if(!outfile)
{
	cerr << "Error in opening output file" << endl;
	return -1;
}

for(int i=0; i<N; i++)
	outfile << i << "   "  << sum_prog[i]-1 << "   "  << err_prog[i] << endl;

outfile.close();








/*

//////////////////////////////////////////////////////////////PART 1

///////////////////////////////Declare variables, prepare object rnd

fstream outfile;

int M=100000;           //Total number of throws
int N=100;              //Number of blocks
int L=M/N;            	//Number of throws in each block, please use for M a multiple of N

double r[M];
double x[N];
double *ave = new double[N];
double *ave2 = new double[N];
double *sum_prog = new double[N];
double *sum2_prog = new double[N];
double *err_prog = new double[N];

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

///////////////////////////////Fill r with M casual numbers, let x=(0,1...N-1) and set all other vectors to 0 

for(int i=0; i<M; i++)			
	r[i] = rnd.Rannyu();
rnd.SaveSeed();

for(int i=0; i<N; i++)
{			
	x[i] = i;		
	ave[i] = 0;		
	ave2[i] = 0;			
	sum_prog[i] = 0;			
	sum2_prog[i] = 0;		
	err_prog[i] = 0;
}

///////////////////////////////Compute average for each of the N blocks and apply function data_blocking

for(int i=0; i<N; i++)	
{
	double sum = 0;
	for(int j=0; j<L; j++)
	{
		int k = i*L+j;
		sum += (M_PI/2)*cos(M_PI*r[k]/2);
	}
	ave[i] = sum/(1.*L);
	ave2[i] = (ave[i])*(ave[i]);
}

data_blocking (ave, ave2, sum_prog, sum2_prog, err_prog, N);

///////////////////////////////Write to file

outfile.open("output_1.txt",ios::out);

if(!outfile)
{
	cerr << "Error in opening output file" << endl;
	return -1;
}

for(int i=0; i<N; i++)
	outfile << x[i] << "   "  << sum_prog[i]-1 << "   "  << err_prog[i] << endl;

outfile.close();



//////////////////////////////////////////////////////////////PART 2

double s[M];

for(int i=0; i<N; i++)
{					
	ave[i] = 0;		
	ave2[i] = 0;			
	sum_prog[i] = 0;			
	sum2_prog[i] = 0;		
	err_prog[i] = 0;
}

for(int i=0; i<M; i++)			
	s[i] = (1./100)*(101-sqrt(10201-10200*(r[i])));

///////////////////////////////Compute average for each of the N blocks and apply function data_blocking

for(int i=0; i<N; i++)	
{
	double sum = 0;
	for(int j=0; j<L; j++)
	{
		int k = i*L+j;
		sum += ((M_PI/2)*cos(M_PI*s[k]/2))/((100*1./51)*((101*1./100)-s[k]));
	}
	ave[i] = sum/(1.*L);
	ave2[i] = (ave[i])*(ave[i]);
}

data_blocking (ave, ave2, sum_prog, sum2_prog, err_prog, N);

///////////////////////////////Write to file

outfile.open("output_2.txt",ios::out);

if(!outfile)
{
	cerr << "Error in opening output file" << endl;
	return -1;
}

for(int i=0; i<N; i++)
	outfile << x[i] << "   "  << sum_prog[i]-1 << "   "  << err_prog[i] << endl;

outfile.close();







delete []ave;
delete []ave2;
delete []sum_prog;
delete []sum2_prog;
delete []err_prog;*/

rnd.SaveSeed();

return 0;
}

