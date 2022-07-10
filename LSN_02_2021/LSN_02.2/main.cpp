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




//////////////////////////////////////////////////////////////Initialize

fstream outfile;

double a = 1;		//Lattice constant
int M = 10000;		//Number of "throws" (copies of walker)
int Nblk = 100;		//blocks into which to sort the M walkers
int L = M/Nblk;		//block size
int steps = 100;	//Time steps

vector <double> x, y, z, r;
vector <double> ave, sum_prog, err_prog;
vector <double> r_avg, r_avg_err;
double rand_sign = 0;
double rand_dir = 0;
double sum = 0;
int sign = 0;

for(int i=0; i<M; i++)
{
	x.push_back(0);
	y.push_back(0);
	z.push_back(0);
        r.push_back(0);
}

for(int i=0; i<Nblk; i++)
{
	ave.push_back(0);
        sum_prog.push_back(0);
        err_prog.push_back(0);
}

for(int i=0; i<steps; i++)
{
        r_avg.push_back(0);
        r_avg_err.push_back(0);
}




//////////////////////////////////////////////////////////////PART 1 (cubic lattice)
    
for(int t=0; t<steps; t++)
{
	cout << "TIME STEP " << t << endl;

	for(int i=0; i<M; i++)
		r[i] = sqrt(pow(x[i],2) + pow(y[i],2) + pow(z[i],2));

	//data blocking on vector r
	int k=0;
	for(int blk=0; blk<Nblk; blk++)	
	{
		sum = 0;
		for(int j=0; j<L; j++)
		{
			k = blk*L+j;
			sum += r[k];
		}
       		ave[blk] = sum/L;
	}

	data_blocking(ave, sum_prog, err_prog);

	//keep track of last element in sum_prog, err_prog, which will be the quantity to plot for every time step
	r_avg[t] = sum_prog[Nblk-1];
	r_avg_err[t] = err_prog[Nblk-1];

	//move
	for(int i=0; i<M; i++)
	{
		//Decide whether particle moves forward (+1) or backward (-1)
		rand_sign = rnd.Rannyu();
		if(rand_sign<0.5) sign = 1;
		else sign = -1;

		//Decide direction of motion x/y/z
		rand_dir = rnd.Rannyu(0,3);
		if(rand_dir < 1)			x[i] = x[i] + sign*a;	//these (x,y,z) specify j-th move: always either (+-1,0,0) or (0,+-1,0) or (0,0,+-1)
		else if(rand_dir >= 1 && rand_dir < 2) 	y[i] = y[i] + sign*a;
		else					z[i] = z[i] + sign*a;
	}
}

outfile.open("discrete.txt",ios::out);

for(int t=0; t<steps; t++)
	outfile << t << " "  << r_avg[t] << " "  << r_avg_err[t] << endl;

outfile.close();




//////////////////////////////////////////////////////////////PART 2 (continuum)

double theta = 0;
double phi = 0;

for(int i=0; i<M; i++)
{
	x[i]=(0);
	y[i]=(0);
	z[i]=(0);
        r[i]=(0);
}

for(int i=0; i<Nblk; i++)
{
	ave[i]=(0);
        sum_prog[i]=(0);
        err_prog[i]=(0);
}

for(int i=0; i<steps; i++)
{
        r_avg[i]=(0);
        r_avg_err[i]=(0);
}


for(int t=0; t<steps; t++)
{
	cout << "TIME STEP " << t << endl;

	for(int i=0; i<M; i++)
	{
		r[i] = sqrt(pow(x[i],2) + pow(y[i],2) + pow(z[i],2));
		//cout << "part: " << i << " r=" <<  r[i] << "    ";
	}
	//cout << endl;


	//data blocking on vector r
	int k=0;
	for(int blk=0; blk<Nblk; blk++)	
	{
		sum = 0;
		for(int j=0; j<L; j++)
		{
			k = blk*L+j;
			sum += r[k];
		}
       		ave[blk] = sum/L;
	}

	data_blocking(ave, sum_prog, err_prog);
	//for(int tr=0; tr<Nblk; tr++) cout << sum_prog[tr] << " ";
	//cout << endl;

	//keep track of last element in sum_prog, err_prog, which will be the quantity to plot for every time step
	r_avg[t] = sum_prog[Nblk-1];
	r_avg_err[t] = err_prog[Nblk-1];
	//cout << "r avg " << r_avg[t] << endl;

	//move
	for(int i=0; i<M; i++)
	{
		//Uniform sampling of solid angle
		theta = 2*M_PI*rnd.Rannyu();
		phi = M_PI*rnd.Rannyu();
		
		x[i] = x[i] + cos(theta)*sin(phi);	
		y[i] = y[i] + sin(theta)*sin(phi);
		z[i] = z[i] + cos(phi);

		//cout << i << " sign: " << sign << "  dir: " << rand_dir << "  x: " << x[i] << "  y: " << y[i] << "  z: " << z[i] << endl;
	}
	//cout << "===============================================" << endl;
}

outfile.open("continuum.txt",ios::out);

for(int t=0; t<steps; t++)
	outfile << t << " "  << r_avg[t] << " "  << r_avg_err[t] << endl;

outfile.close();

rnd.SaveSeed();

return 0;
}

