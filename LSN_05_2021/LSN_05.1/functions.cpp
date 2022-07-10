#include "functions.h"
#include <cmath>
#include <vector>

using namespace std;

double error(vector <double> v1, vector <double> v2, int n)
{
	if(n==0) return 0;
	else return sqrt((v1[n]-pow(v2[n],2))/n);
}



void data_blocking(vector <double> value, vector <double> &sum_prog, vector <double> &err_prog)
{
	vector <double> sum2_prog;
	int N=sum_prog.size();
    
	for(int i=0;i<N;i++)
	{
		sum_prog[i]=0;
		sum2_prog.push_back(0);
	}
    
	for(int i=0;i<N;i++)
	{
        	for(int j=0;j<i+1;j++)
		{
			sum_prog[i]=sum_prog[i]+value[j];
			sum2_prog[i]=sum2_prog[i]+pow(value[j],2);
		}

        	sum_prog[i]=sum_prog[i]/(i+1);
        	sum2_prog[i]=sum2_prog[i]/(i+1);
        	err_prog[i] = error(sum2_prog,sum_prog,i);
	}
}



double psi2_100(double x, double y, double z)			//in units of Bohr radius
{
	return exp(-2*sqrt(pow(x,2)+pow(y,2)+pow(z,2)))/M_PI;
}



double psi2_210(double x, double y, double z)			//in units of Bohr radius
{
	return pow(z,2)*exp(-sqrt(pow(x,2)+pow(y,2)+pow(z,2)))/(32*M_PI);
}


/*
void Metropolis_uniform(vector <double> x, vector <double> y, vector <double> z, vector <double> r, int i, double q_xy, double delta)
{
	//Randomize position according to uniform T
	x[i] = x[i-1] + rnd.Rannyu(x[i-1]-delta,x[i-1]+delta);
	y[i] = y[i-1] + rnd.Rannyu(x[i-1]-delta,x[i-1]+delta);
	z[i] = z[i-1] + rnd.Rannyu(x[i-1]-delta,x[i-1]+delta);
	r[i] = sqrt(pow(x[i],2)+pow(y[i],2)+pow(z[i],2));

	//T is symmetric, so q(x|y) = p(x)/p(y) = |psi(x)|^2/|psi(y)|^2
        q_xy = psi2_100(x[i],y[i],z[i])/psi2_100(x[i-1],y[i-1],z[i-1]);

	//if q>1 then alpha = min{1,q} = 1 and the move is accepted, so nothing more to do. If q<1 and casual number >q, the move is rejected
        if(q_xy<1)
	{
		if(rnd.Rannyu()>q_xy)
		{
			x[i] = x[i-1];
			y[i] = y[i-1];
			z[i] = z[i-1];
			r[i] = r[i-1];
		}
	}
}



void Metropolis_gauss(vector <double> x, vector <double> y, vector <double> z, vector <double> r, int i, double q_xy, double delta)
{
	//Randomize position according to uniform T
	x[i] = x[i-1] + rnd.Gauss(x[i-1],delta);
	y[i] = y[i-1] + rnd.Gauss(x[i-1],delta);
	z[i] = z[i-1] + rnd.Gauss(x[i-1],delta);
	r[i] = sqrt(pow(x[i],2)+pow(y[i],2)+pow(z[i],2));

	//T is symmetric, so q(x|y) = p(x)/p(y) = |psi(x)|^2/|psi(y)|^2
        q_xy = psi2_100(x[i],y[i],z[i])/psi2_100(x[i-1],y[i-1],z[i-1]);

	//if q>1 then alpha = min{1,q} = 1 and the move is accepted, so nothing more to do. If q<1 and casual number >q, the move is rejected
        if(q_xy<1)
	{
		if(rnd.Rannyu()>q_xy)
		{
			x[i] = x[i-1];
			y[i] = y[i-1];
			z[i] = z[i-1];
			r[i] = r[i-1];
		}
	}
}*/
 
