#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include "random.h"
#include "functions.h"

using namespace std;

int main (int argc, char *argv[]){

ofstream outfile;

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

/////////////////////////////////////////////////////////////"Lay the table!"

const int n_bin = 100;		//box LxL where L=n_bin*d
const double d = 1.;		//lattice spacing
const double L = 0.9;		//needle lenght 
double r,s,theta,t,u,mod,a,b;	//coordinates
int N_run = 100;		//n° of MonteCarlo runs
int N_thr = 10000;		//n° throws in each experiment
int N_hit;			//counter

vector <double> pi; 		//vector that contains estimates of pi
vector <double> sum_prog;
vector <double> err_prog;

for(int k=0; k<N_run; k++)
{
	pi.push_back(0);    
	sum_prog.push_back(0); 
	err_prog.push_back(0); 
}   		

/////////////////////////////////////////////////////////////Perform N_run MonteCarlo runs

for(int k=0; k<N_run; k++)
{

	//Reset counter	
	N_hit=0;

	for(int i=0; i<N_thr; i++)
	{

		/*r = rnd.Rannyu(0,n_bin*d); 
		s = rnd.Rannyu(0,n_bin*d); 

		do
		{
			theta = rnd.Rannyu(0,2*M_PI);
			t = r + L*cos(theta);
			u = s + L*sin(theta);
		}
		while(t<0 || t>n_bin*d || u<0 || u>n_bin*d);*/

		//Launch coordinates of vertices: (r,s) casually in [0,n_bin*d], (t,u) with constraints 1) needle length = L 2) needle must be fully contained in box

	    	r=rnd.Rannyu(0,n_bin*d);
            	s=rnd.Rannyu(0,n_bin*d);
            
            	do
	    	{
            		a=(r-L)+2*L*rnd.Rannyu();
               		b=(s-L)+2*L*rnd.Rannyu();
                	mod=sqrt(pow(a-r,2)+pow(b-s,2));
           	}
		while(mod>L);	//while(mod>L || a<0 || a>n_bin*d || b<0 || b>n_bin*d)
				//Strange behaviour: works fine (estimated pi compatible with 3.14....) when surface effects are ignored, i.e. when while(mod>L) is used
				//Estimated pi no longer compatible with 3.14.... when surface effects are taken into account, i.e. while(mod>L || a<0 || a>n_bin*d || b<0 || b>n_bin*d)
            
		if(b>=s)
		{
			t=r+L*(a-r)/mod;
			u=s+L*sin(acos((a-r)/mod));
		}
		else
		{
			t=r+L*cos(2*M_PI-acos((a-r)/mod));
			u=s+L*sin(2*M_PI-acos((a-r)/mod));
		}
            
		//Run over "bins", locate vertex r in bin k, increment counter if t is located in one of the nearest-neighbour bins
		for(int k=0; k<n_bin; k++)
		{
			if( s>=(double)k*d & s<(double)((k+1)*d) ) 
			{
				if(	( u>=(double)((k+1)*d) & u<(double)((k+2)*d) ) || ( u>=(double)((k-1)*d) & u<(double)((k)*d) )		)
				N_hit++;
			}
		}
	}

	pi[k] = (2*L/d)*(N_thr*1./N_hit);
}

/////////////////////////////////////////////////////////////Data blocking and write to output

data_blocking(pi, sum_prog, err_prog);

outfile.open("output.txt",ios::out);

if(!outfile)
{
	cerr << "Error in opening output file" << endl;
	return -1;
}

    
for(int k=0; k<N_run; k++)
	outfile << k << " " << sum_prog[k] << " " << err_prog[k] << endl;

outfile.close();

rnd.SaveSeed();

return 0;
}

