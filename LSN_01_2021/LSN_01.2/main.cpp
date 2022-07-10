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

vector <double> s;  //Standard dice
vector <double> e;  //Exponential dice
vector <double> l;  //Lorentzian dice
       
///////////////////////////////////////////////////N=1

for(int i=0; i<10000; i++)
{
        s.push_back(rnd.Rannyu());
        e.push_back(rnd.Exp(1));
        l.push_back(rnd.Lorentz(0,1));
}

outfile.open("output_1.txt",ios::out);

if(!outfile)
{
	cerr << "Error in opening output file" << endl;
	return -1;
}
  
for(int k=0; k<10000; k++)
        outfile << s[k] << "          " << e[k] << "          " << l[k] << endl;

outfile.close();

///////////////////////////////////////////////////N=2

//v1,v2,v3 contain random variables x_i to be averaged. In view of N=10, N=100 case (involving same procedure), these vectors are created with dimension 100 


vector <double> v1,v2,v3;

for(int i=0; i<100; i++)
{
        v1.push_back(0);
        v2.push_back(0);
        v3.push_back(0);
}
    
for(int i=0; i<10000; i++)
{
        for(int j=0; j<2; j++)
	{
		v1[j] = rnd.Rannyu();
		v2[j] = rnd.Exp(1);
		v3[j] = rnd.Lorentz(0,1);
        }

        s[i]=average(v1,2);
        e[i]=average(v2,2);
        l[i]=average(v3,2);
}
    
outfile.open("output_2.txt",ios::out);

if(!outfile)
{
	cerr << "Error in opening output file" << endl;
	return -1;
}
    
for(int k=0; k<10000; k++)
        outfile << s[k] << "          " << e[k] << "          " << l[k] << endl;

outfile.close();

///////////////////////////////////////////////////N=10

for(int i=0; i<10000; i++)
{
        for(int j=0; j<10; j++)
	{
		v1[j] = rnd.Rannyu();
		v2[j] = rnd.Exp(1);
		v3[j] = rnd.Lorentz(0,1);
        }

        s[i]=average(v1,10);
        e[i]=average(v2,10);
        l[i]=average(v3,10);
}
    
outfile.open("output_10.txt",ios::out);

if(!outfile)
{
	cerr << "Error in opening output file" << endl;
	return -1;
}
    
for(int k=0; k<10000; k++)
        outfile << s[k] << "          " << e[k] << "          " << l[k] << endl;

outfile.close();

///////////////////////////////////////////////////N=100

for(int i=0; i<10000; i++)
{
        for(int j=0; j<100; j++)
	{
		v1[j] = rnd.Rannyu();
		v2[j] = rnd.Exp(1);
		v3[j] = rnd.Lorentz(0,1);
        }

        s[i]=average(v1,100);
        e[i]=average(v2,100);
        l[i]=average(v3,100);
}
    
outfile.open("output_100.txt",ios::out);

if(!outfile)
{
	cerr << "Error in opening output file" << endl;
	return -1;
}
    
for(int k=0; k<10000; k++)
        outfile << s[k] << "          " << e[k] << "          " << l[k] << endl;

outfile.close();
  
rnd.SaveSeed();

return 0;
}

