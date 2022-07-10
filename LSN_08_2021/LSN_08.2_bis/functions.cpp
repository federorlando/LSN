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



double psi2(double x, double sigma, double mi)			
{
	return exp(-(pow(x-mi,2))/(sigma*sigma)) + exp(-(pow(x+mi,2))/(sigma*sigma)) + 2*exp(-(x*x+mi*mi)/(sigma*sigma));
}



double Hpsi_psi(double x, double sigma, double mi)			
{
	double e_minus = exp(-(pow(x-mi,2))/(2*sigma*sigma));
	double e_plus = exp(-(pow(x+mi,2))/(2*sigma*sigma));
	double e_sum = e_plus + e_minus;
	double first = pow(x,4) - 2.5*pow(x,2) + (1-pow(((x-mi)/sigma),2))/(2*sigma*sigma);
	double second = pow(x,4) - 2.5*pow(x,2) + (1-pow(((x+mi)/sigma),2))/(2*sigma*sigma);

	return (e_minus*first + e_plus*second)/e_sum;
}



/*double Minimum(vector <double> value) 							
{
	double min = value[0];
	unsigned int N = value.size();
    
    	for (unsigned int i=0; i<N; i++) 
	{
        	if(value[i] < min) 
	  		min = value[i];	
    	}
    
	return min;    
}*/




