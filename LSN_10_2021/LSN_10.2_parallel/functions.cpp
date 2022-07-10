#include "functions.h"
#include "TSP_position.h"
#include <cmath>
#include <vector>
#include <algorithm>

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



double TSP_cost(vector <int> o, TSP_position p[])						
{
	double cost = 0;
	int n = o.size();
	double t = 0;

	for(int i=1; i<n; i++)
	{
		t = sqrt(pow(p[o[i]-1].GetX()-p[o[i-1]-1].GetX(),2)+pow(p[o[i]-1].GetY()-p[o[i-1]-1].GetY(),2));	// += d(r_i, r_{i-1})
		//cout << << t << endl;
		cost += t;
	}
	
	cost += sqrt(pow(p[o[0]-1].GetX()-p[o[n-1]-1].GetX(),2)+pow(p[o[0]-1].GetY()-p[o[n-1]-1].GetY(),2)); 		//PBC: r_(n-1) = r_0
	return cost;
}



bool TSP_check(vector <int> o)
{
	int val = 0;
	for(int i=0; i<o.size(); i++)
	{
        	if(count(o.begin(),o.end(),i+1)>1)	//Returns the number of elements in the range [first,last) that compare equal to val.
		{
        		val=1;
        		break;
        	}
	}
	if(val==0) return true;     //Bonds respected
	else return false;          //Bonds not respected
}




void swap (vector <int> &vett, int a, int b)
{
	double box = vett[a];
	vett[a] = vett[b];
	vett[b] = box;
}




void shift (vector <int> &vett, int shift)			//E.g. [v0,v1,v2,v3,v4] shift=2 -> [v3,v4,v0,v1,v2]
{
	vector <int> box;
	
	if(shift<=vett.size()-1)
	{
		for(int i=0; i<vett.size(); i++)
			box.push_back(vett[i]);
	
		for(int i=1; i<vett.size()-shift; i++)
			vett[i+shift] = box[i];
	
		for(int i=1; i<(shift+1); i++)
			vett[i] = box[vett.size()-1-shift+i];
	}
}



void invert_order(vector <int> &vett, int a, int b)		//a included, b included. E.g. [v0,v1,v2,v3,v4] a=1, b=3 -> [v0,v3,v2,v1,v4]
{     
	vector <int> box;
	
	if((b-a+1)<=vett.size())
	{
		for(int i=0; i<(b-a+1); i++)
			box.push_back(vett[b-i]);
		
		for(int i=0; i<(b-a+1); i++)
        		vett[a+i] = box[i];
	}
}




void part_perm(vector <int> &vett, int a, int b, int c, int d)
{
	vector <int> box;

	if((a<=b) & (b<=c) & (c<=d))
	{
		for(int i=0; i<(b-a); i++)
        		box.push_back(vett[a+i]);
    
    		for(int i=0; i<(d-c); i++)
        		box.push_back(vett[c+i]);
    	
		random_shuffle(box.begin(),box.end());
    
		for(int i=0; i<(b-a); i++)
			vett[a+i] = box[i];

		for(int i=0; i<(d-c); i++)
        		vett[c+i] = box[b-a+i];
	}
}




void crossover(vector <int> &v1, vector <int> &v2, int cut)
{
	vector <int> box1, box2;

	for(int i=cut; i<v1.size(); i++)		//copy input vectors, starting from index = cut till end
	{
	        box1.push_back(v1[i]);
	        box2.push_back(v2[i]);
	}

	vector <int> pos;
	for(int i=0; i<box1.size(); i++)
	{
        	for(int j=0;j<v2.size();j++)
		{
        		if(box1[i]==v2[j])
			{
	        	        pos.push_back(j);
	        	        break;
        	    	}
        	}
    	}
    
	sort(pos.begin(),pos.end());

	for(int i=0;i<pos.size();i++)
	{
	        v1[cut+i]=v2[pos[i]];
	}
    
	for(auto el: pos) el=0;
    
	for(int i=0;i<box2.size();i++)
	{
	        for(int j=0;j<v1.size();j++)
		{
			if(box2[i]==v1[j])
			{
		                pos[i]=j;
		                break;
	            	}
	        }
	}
    
	sort(pos.begin(),pos.end());

	for(int i=0;i<pos.size();i++)
		v2[cut+i]=v1[pos[i]];
}



double expo(double y, double lambda)
{
	return -(1/lambda)*log(1-y);
}




int selection(vector <double> fitness, double rand)
{
	double min_fit = *min_element(fitness.begin(),fitness.end());	//minimum fitness
	double val = min_fit + expo(rand,10);				//to minimum fitness sum number exp. distributed
	int index = 0;
	double min_diff = abs(fitness[0]-val);

	for(int k=1; k<fitness.size(); k++)				//check which chromosome has fitness closest to (minimum fitness + above number)
	{
	        if(abs(fitness[k]-val) < min_diff)
		{
	        	min_diff=abs(fitness[k]-val);
	       		index=k;
	        }
    	}
    
	return index;							//select that one
}




