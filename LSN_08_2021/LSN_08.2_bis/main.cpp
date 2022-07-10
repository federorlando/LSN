#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include "random.h"
#include "functions.h"

using namespace std;

int main (int argc, char *argv[]){

////////////////////////////////Prepare object rnd

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

////////////////////////////////Declare variables

double sigma, mi;
ofstream outfile_r, landscape, outfile_pdf;

int M = 1000000;			//Number of throws
int M_equil = 1000;			//Number of throws during equilibration phase
int N = 100;				//Number of blocks
int L = M/N;				//Number of throws per block

vector <double> x, x_equil;		

for(int i=0; i<M; i++) x.push_back(0);
for(int i=0; i<M_equil; i++) x_equil.push_back(0);

vector <double> ave, sum_prog, err_prog;
    
for(int i=0; i<N; i++)
{
	ave.push_back(0);
	sum_prog.push_back(0);
	err_prog.push_back(0);
}

double q_xy = 0;			//Function q(x|y)
double delta_uniform = 1;		//Width of Metropolis step (uniform T); start e.g. with delta = 1

ifstream ReadInput;
ReadInput.open("input.dat"); 
ReadInput >> sigma;
ReadInput >> mi;
ReadInput.close();

////////////////////////////////Estimate delta

int n_reject = 0;			//Counter for number of rejected moves
double delta_trial = 0.01;     		//Successive estimates for delta
double diff = 0.5;			//Deviation from 50% target; this will never be >0.5, so initial value may be safely set to 0.5
int N_moves = 1000;
//Histogram for the pdf
int nbins = 100;
double range = 6.;        		//pdf in range [-3.,3)
double bin_size = range/nbins;
vector <double> pdf;
for(int i=0; i<nbins; i++) pdf.push_back(0.);

for(int j=0; j<2000; j++)
{
	x[0] = 0;
        n_reject = 0;

        for(int i=1; i<N_moves; i++)	//Metropolis run (see further on for details)
	{
       		x[i] = x[i-1]-delta_trial+2*delta_trial*rnd.Rannyu();

		q_xy = psi2(x[i],sigma,mi)/psi2(x[i-1],sigma,mi);
		
		if(q_xy<1 && rnd.Rannyu()>q_xy)
		{
                   	x[i]=x[i-1];
                    	n_reject++;
           	}
        }
        
        if(abs(0.5-1.*n_reject/N_moves)<diff) 		//if rejection probability differs from 50% target less than a certain deviation...
	{
		delta_uniform = delta_trial;		//...store current value of delta_trial...
		diff = abs(0.5-1.*n_reject/N_moves);	//...and fine-tune deviation
        }
        delta_trial += 0.005;
}

cout << "Best estimate for delta:				" << delta_uniform << endl;
cout << "Rejection probability differs from 50% of amount:	" << diff << endl;

////////////////////////////////Equilibration

x_equil[0] = 0;
    
for(int i=1; i<M_equil; i++)
{
	//Randomize position according to uniform T
        x_equil[i] = x_equil[i-1]-delta_uniform+2*delta_uniform*rnd.Rannyu();

	//T is symmetric, so q(x|y) = p(x)/p(y) = |psi(x)|^2/|psi(y)|^2
        q_xy = psi2(x_equil[i],sigma,mi)/psi2(x_equil[i-1],sigma,mi);

	//if q>1 then alpha = min{1,q} = 1 and the move is accepted, so nothing more to do. If q<1 and casual number >q, the move is rejected
        if(q_xy<1 && rnd.Rannyu()>q_xy)
	{
		x_equil[i] = x_equil[i-1];
	}
}

////////////////////////////////Main run

x[0] = x_equil[M_equil-1];
n_reject = 0;
    
for(int i=1; i<M; i++)
{
	//Randomize position according to uniform T
        x[i] = x[i-1]-delta_uniform+2*delta_uniform*rnd.Rannyu();

	//T is symmetric, so q(x|y) = p(x)/p(y) = |psi(x)|^2/|psi(y)|^2
        q_xy = psi2(x[i],sigma,mi)/psi2(x[i-1],sigma,mi);

	//if q>1 then alpha = min{1,q} = 1 and the move is accepted, so nothing more to do. If q<1 and casual number >q, the move is rejected
        if(q_xy<1 && rnd.Rannyu()>q_xy)
	{
		x[i] = x[i-1];
                n_reject++;
	}

        //Filling the histogram
        for(int j=0;j<nbins;j++)
	{
		if(x[i] >= (-3+j*bin_size) && x[i] < (-3+(j+1)*bin_size))
		{
	                pdf[j] += 1.;
	                break;
		}
        }
}

cout << "Main run: fraction of rejected moves =			" << 1.*n_reject/M << endl;

////////////////////////////////Data blocking
    
int k=0;
for(int i=0; i<N; i++)
{
        for(int j=0; j<L; j++)
	{
		k = j + i*L;
		ave[i] += Hpsi_psi(x[k],sigma,mi);
        }

        ave[i]/=L;
}

data_blocking(ave,sum_prog,err_prog);

////////////////////////////////Write to file   

outfile_r.open("output_sigma_mi",ios::out);
for(int i=0; i<N; i++)
	outfile_r << i << "   "  << sum_prog[i] << "   "  << err_prog[i] << endl;

outfile_r.close();

landscape.open("E0_sigma_mi",ios::app);
landscape << sigma << "   "  << mi << "   "  << sum_prog[N-1] << endl;
landscape.close();

outfile_pdf.open("./pdf.txt");
    
//PDF normalized in jupyter NB
for(int i=0; i<nbins; i++)
        outfile_pdf <<-3+(i*bin_size) << " " << pdf[i] << endl;
outfile_pdf.close();

rnd.SaveSeed();

return 0;
}

