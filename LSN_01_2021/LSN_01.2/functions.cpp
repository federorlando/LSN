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



double average(vector <double> v, int n)
{
	double sum=0;

	for(int i=0; i<n; i++)
        	sum += v[i];

	return sum/(1.*n);
}





