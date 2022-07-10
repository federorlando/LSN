#!/bin/bash

for ((alat=0.3;((alat<=0.4));alat=alat+0.1))
do 

#i=$(dc -e "3k $alat 0.001*p")

cat > main.cpp << EOF
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include "random.h"
#include "functions.h"

// INSTRUCTIONS
// The code is designed in order to sample with a fixed T (gaussian OR normal)
// To switch from one to the other:
// 1) in preventive Metropolis run, decomment "x[i] = x[i-1]-delta_trial+2*delta_trial*rnd.Rannyu()" and comment "x[i] = x[i-1]+rnd.Gauss(0,delta_gauss)" or v.v. (and same for y,z)
// 2) in preventive Metropolis run, decomment "delta_uniform = delta_trial" and comment "delta_gauss = delta_trial" or v.v.
// 3) in preventive Metropolis run, suitably decomment and comment cout lines 
// 4) in Metropolis run, repeat point 1
// 5) lastly, remember to change the name of the output files!!!!!

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

////////////////Declare variables

ofstream outfile_r;
ofstream outfile_points;

int M = 1000000;			//Number of throws
int N = 100;				//Number of blocks
int L = M/N;				//Number of throws per block

vector <double> x, y, z, r;		//Coordinates x,y,z and distance r (variable to be estimated)

for(int i=0; i<M; i++)
{
	x.push_back(0);
        y.push_back(0);
        z.push_back(0);
        r.push_back(0);
}

vector <double> ave, sum_prog, err_prog;
    
for(int i=0; i<N; i++)
{
	ave.push_back(0);
	sum_prog.push_back(0);
	err_prog.push_back(0);
}

double q_xy = 0;			//Function q(x|y)
double delta_uniform = 1;		//Width of Metropolis step (uniform T); start e.g. with delta = 1
double delta_gauss = 1;			//Width of Metropolis step (gaussian T); start e.g. with delta = 1





///////////////////////////////////////////////////////////////PART 1: GROUND STATE |0,0,0> WITH UNIFORM T(x|y) (for gaussian, find in the following suitable code lines to decomment)

////////////////Estimate delta with preventive Metropolis run

int n_reject = 0;			//Counter for number of rejected moves
double delta_trial = $alat;     	//Successive estimates for delta
double diff = 0.5;			//Deviation from 50% target; this will never be >0.5, so initial value may be safely set to 0.5

for(int j=0; j<100; j++)
{
	x[0] = 0;
	y[0] = 0;
	z[0] = 0; 
        n_reject = 0;

        for(int i=1; i<10000; i++)	//Metropolis run (see further on for details)
	{
       		x[i] = x[i-1]-delta_trial+2*delta_trial*rnd.Rannyu();
        	y[i] = y[i-1]-delta_trial+2*delta_trial*rnd.Rannyu();
        	z[i] = z[i-1]-delta_trial+2*delta_trial*rnd.Rannyu();
		/*x[i] = x[i-1]+rnd.Gauss(0,delta_gauss);
		y[i] = y[i-1]+rnd.Gauss(0,delta_gauss);
		z[i] = z[i-1]+rnd.Gauss(0,delta_gauss);*/

		q_xy = psi2_100(x[i],y[i],z[i])/psi2_100(x[i-1],y[i-1],z[i-1]);
		
		if(q_xy<1 && rnd.Rannyu()>q_xy)
		{
                   	x[i]=x[i-1];
                    	y[i]=y[i-1];
                    	z[i]=z[i-1];
                    	n_reject++;
           	}
        }
        
        if(abs(0.5-1.*n_reject/10000)<diff) 		//if rejection probability differs from 50% target less than a certain deviation...
	{
		delta_uniform = delta_trial;		//...store current value of delta_trial...
		//delta_gauss = delta_trial;
		diff = abs(0.5-1.*n_reject/10000);	//...and fine-tune deviation
        }
        delta_trial += 0.005;
}

cout << "PSI=|1,0,0>, T=uniform" << endl;
cout << "Best estimate for delta:				" << delta_uniform << endl;
cout << "Rejection probability differs from 50% of amount:	" << diff << endl << endl;

/*cout << "PSI=|1,0,0>, T=gaussian" << endl;
cout << "Best estimate for delta:				" << delta_gauss << endl;
cout << "Rejection probability differs from 50% of amount:	" << diff << endl << endl;*/

////////////////Build random walk Metropolis  

x[0] = 0;
y[0] = 0;
z[0] = 0;
r[0] = sqrt(pow(x[0],2)+pow(y[0],2)+pow(z[0],2));   
    
for(int i=1; i<M; i++)
{
	//Randomize position according to uniform (gaussian) T
        x[i] = x[i-1]-delta_uniform+2*delta_uniform*rnd.Rannyu();
        y[i] = y[i-1]-delta_uniform+2*delta_uniform*rnd.Rannyu();
        z[i] = z[i-1]-delta_uniform+2*delta_uniform*rnd.Rannyu();
	/*x[i] = x[i-1]+rnd.Gauss(0,delta_gauss);
	y[i] = y[i-1]+rnd.Gauss(0,delta_gauss);
	z[i] = z[i-1]+rnd.Gauss(0,delta_gauss);*/
	r[i] = sqrt(pow(x[i],2)+pow(y[i],2)+pow(z[i],2));

	//T is symmetric, so q(x|y) = p(x)/p(y) = |psi(x)|^2/|psi(y)|^2
        q_xy = psi2_100(x[i],y[i],z[i])/psi2_100(x[i-1],y[i-1],z[i-1]);

	//if q>1 then alpha = min{1,q} = 1 and the move is accepted, so nothing more to do. If q<1 and casual number >q, the move is rejected
        if(q_xy<1 && rnd.Rannyu()>q_xy)
	{
		x[i] = x[i-1];
		y[i] = y[i-1];
		z[i] = z[i-1];
		r[i] = r[i-1];
	}
}

////////////////Data blocking
    
int k=0;
for(int i=0; i<N; i++)
{
        for(int j=0; j<L; j++)
	{
		k = j + i*L;
		ave[i] += r[k];
        }

        ave[i]/=L;
}

data_blocking(ave,sum_prog,err_prog);

////////////////Write to file   

outfile_r.open("output_100_uniform.txt",ios::out);
for(int i=0; i<N; i++)
	outfile_r << i << "   "  << sum_prog[i]-1.5 << "   "  << err_prog[i] << endl;

outfile_points.open("output_100_uniform_points.txt",ios::out);
for(int i=0; i<M; i++)
        outfile_points << x[i] << "   " << y[i] << "   " << z[i] << endl;

outfile_r.close();
outfile_points.close();





///////////////////////////////////////////////////////////////PART 2: EXCITED STATE |2,1,0> WITH UNIFORM T(x|y) (for gaussian, find in the following suitable code lines to decomment)

for(int i=0; i<N; i++)
{
	ave[i]=0;
	sum_prog[i]=0;
	err_prog[i]=0;
}

////////////////Estimate delta with preventive Metropolis run

n_reject = 0;			//Counter for number of rejected moves
delta_trial = 3.1;     		//Successive estimates for delta
diff = 0.5;			//Deviation from 50% target; this will never be >0.5, so initial value may be safely set to 0.5

for(int j=0; j<100; j++)
{
	x[0] = 0;
	y[0] = 0;
	z[0] = 0; 
        n_reject = 0;

        for(int i=1; i<10000; i++)	//Metropolis run (see further on for details)
	{
       		x[i] = x[i-1]-delta_trial+2*delta_trial*rnd.Rannyu();
        	y[i] = y[i-1]-delta_trial+2*delta_trial*rnd.Rannyu();
        	z[i] = z[i-1]-delta_trial+2*delta_trial*rnd.Rannyu();
		/*x[i] = x[i-1]+rnd.Gauss(0,delta_gauss);
		y[i] = y[i-1]+rnd.Gauss(0,delta_gauss);
		z[i] = z[i-1]+rnd.Gauss(0,delta_gauss);*/

		q_xy = psi2_210(x[i],y[i],z[i])/psi2_210(x[i-1],y[i-1],z[i-1]);
		
		if(q_xy<1 && rnd.Rannyu()>q_xy)
		{
                   	x[i]=x[i-1];
                    	y[i]=y[i-1];
                    	z[i]=z[i-1];
                    	n_reject++;
           	}
        }
        
        if(abs(0.5-1.*n_reject/10000)<diff) 		//if rejection probability differs from 50% target less than a certain deviation...
	{
		delta_uniform = delta_trial;		//...store current value of delta_trial...
		//delta_gauss = delta_trial;
		diff = abs(0.5-1.*n_reject/10000);	//...and fine-tune deviation
        }
        delta_trial += 0.005;
}

cout << "PSI=|2,1,0>, T=uniform" << endl;
cout << "Best estimate for delta:				" << delta_uniform << endl;
cout << "Rejection probability differs from 50% of amount:	" << diff << endl << endl;

/*cout << "PSI=|2,1,0>, T=gaussian" << endl;
cout << "Best estimate for delta:				" << delta_gauss << endl;
cout << "Rejection probability differs from 50% of amount:	" << diff << endl << endl;*/

////////////////Build random walk Metropolis  

x[0] = 0;
y[0] = 0;
z[0] = 0;
r[0] = sqrt(pow(x[0],2)+pow(y[0],2)+pow(z[0],2));   
    
for(int i=1; i<M; i++)
{
	//Randomize position according to uniform T
        x[i] = x[i-1]-delta_uniform+2*delta_uniform*rnd.Rannyu();
        y[i] = y[i-1]-delta_uniform+2*delta_uniform*rnd.Rannyu();
        z[i] = z[i-1]-delta_uniform+2*delta_uniform*rnd.Rannyu();
	/*x[i] = x[i-1]+rnd.Gauss(0,delta_gauss);
	y[i] = y[i-1]+rnd.Gauss(0,delta_gauss);
	z[i] = z[i-1]+rnd.Gauss(0,delta_gauss);*/
	r[i] = sqrt(pow(x[i],2)+pow(y[i],2)+pow(z[i],2));

	//T is symmetric, so q(x|y) = p(x)/p(y) = |psi(x)|^2/|psi(y)|^2
        q_xy = psi2_210(x[i],y[i],z[i])/psi2_210(x[i-1],y[i-1],z[i-1]);

	//if q>1 then alpha = min{1,q} = 1 and the move is accepted, so nothing more to do. If q<1 and casual number >q, the move is rejected
        if(q_xy<1 && rnd.Rannyu()>q_xy)
	{
		x[i] = x[i-1];
		y[i] = y[i-1];
		z[i] = z[i-1];
		r[i] = r[i-1];
	}
}

////////////////Data blocking
    
k=0;
for(int i=0; i<N; i++)
{
        for(int j=0; j<L; j++)
	{
		k = j + i*L;
		ave[i] += r[k];
        }

        ave[i]/=L;
}

data_blocking(ave,sum_prog,err_prog);

////////////////Write to file   

outfile_r.open("output_210_uniform.txt",ios::out);
for(int i=0; i<N; i++)
	outfile_r << i << "   "  << sum_prog[i]-5 << "   "  << err_prog[i] << endl;

outfile_points.open("output_210_uniform_points.txt",ios::out);
for(int i=0; i<M; i++)
        outfile_points << x[i] << "   " << y[i] << "   " << z[i] << endl;

outfile_r.close();
outfile_points.close();

rnd.SaveSeed();

return 0;
}
EOF

make
./main.exe

#rgrep -e 'lattice parameter' -e ! cr.scf.out.GGA | \
      #rawk '/lattice/{i=$(NF-1)}/!/{print i, $(NF-1)}' >> cr.etot_vs_alat_GGA
	
echo alat $alat

done



