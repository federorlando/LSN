/****************************************************************
*****************************************************************
é    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include <string>
#include <sstream>
#include <vector>
#include "Monte_Carlo_NVT.h"

using namespace std;

int main(){

Input();    //Inizialization

//////////////////////////////////////////////////Estimate delta

int n_reject = 0;			//Counter for number of rejected moves
double reject_rate;
double diff = 0.5;			//Deviation from 50% target; this will never be >0.5, so initial value may be safely set to 0.5
double delta_best = 1;
int N_moves = 1000;

for(int j=0; j<2000; j++)
{
	//cout << endl << "will now try with delta " << delta << endl;
	n_reject = 0;
	sigma = 0.625;
	mi = 0.8;

	for(int i=1; i<N_moves; i++)
	{
		accepted = 0;
		attempted = 0;
		Move();
		n_reject = n_reject + attempted - accepted;
		//cout << "Attempted " << attempted << endl;
		//cout << "Accepted " << accepted << endl;
		//cout << "Rejected " << n_reject << endl;
	}

	reject_rate = 1.*n_reject/N_moves;

	//cout << "Acceptance rate " << 1.*accepted/attempted << endl;
	//cout << "Rejection rate " << 1.*n_reject/attempted << endl;
	//cout << "Spread " << 0.5-1.*reject_rate << endl;
	//cout << "Diff " << diff << endl;

	if(abs(0.5-1.*reject_rate)<diff) 		//if rejection probability differs from 50% target less than a certain deviation...
	{
		delta_best = delta;			//...store current value of delta...
		diff = abs(0.5-1.*reject_rate);		//...and fine-tune deviation
        }
	
	delta += 0.005;

	//cout << "Delta best: " << delta_best << endl;
}

cout << endl << "Best estimate for delta:				" << delta_best << endl;
cout << "Rejection probability differs from 50% of amount:	" << diff << endl << endl;

delta = delta_best;

//////////////////////////////////////////////////Main run 

for(int iblk=1; iblk <= nblk; ++iblk)
{  
	Reset(iblk);    //Reset block averages

	for(int istep=1; istep <= nstep; ++istep)
	{
		Move();
		Measure();
		Accumulate();   //Update block averages
		inst_val_energy[(iblk-1)*nstep+istep] = walker[iv];
		inst_val_sigma[(iblk-1)*nstep+istep] = walker[is];
		inst_val_mi[(iblk-1)*nstep+istep] = walker[im];
	}

	Averages(iblk);     //Print results for current block
}

ConfFinal();    //Write final configuration

ofstream WriteInstEne;
WriteInstEne.open("inst_energy_iter"+to_string(niter));
for (int i=0; i<nstep*nblk; ++i) WriteInstEne << inst_val_energy[i] << endl;
WriteInstEne.close();

ofstream outfile_points;
outfile_points.open("points_iter"+to_string(niter));
for(int i=0; i<nblk*nstep; i++) outfile_points << inst_val_sigma[i] << "   " << inst_val_mi[i] << endl;
outfile_points.close();


////////////////////////////////////////////////////////////////////
//		   ESTIMATION OF BEST SIGMA, MI			  //
////////////////////////////////////////////////////////////////////

/*
int M = 1000;
int Nbin = 100;
double L = 0.1;
double binsize = L/Nbin;
double lowlim_sigma = 0.55;
double lowlim_mi = 0.75;
//double uplim_sigma = 0.65;
//double uplim_mi = 0.85;

vector <vector<int>> Matrix(Nbin, vector<int>(Nbin));
vector <double> sigma;
vector <double> mi;
double value_mi, value_sigma;
double Max_value, Max_sigma, Max_mi, index_i, index_j; 

for(int i=0; i<M; i++)
{
        sigma.push_back(0);
	mi.push_back(0);
}

ifstream ReadPoints;
ReadPoints.open("./points_iter10000");
for(int t=0; t<M; t++) ReadPoints >> sigma[t] >> mi[t];
ReadPoints.close();


for(int i=0; i<M; i++)
{
	for(int j=0; j<M; j++)
	{
        	value_sigma = sigma[i];
        	value_mi = mi[j];
        	for(int k=0; k<Nbin; k++)
		{
            		for(int l=0; l<Nbin; l++)
			{
                		if ((value_sigma >= (lowlim_sigma+k*binsize)) & (value_sigma < (lowlim_sigma+(k+1)*binsize)) & (value_mi >= (lowlim_mi+l*binsize)) & (value_mi < (lowlim_mi+(l+1)*binsize))) Matrix[k][l]=Matrix[k][l]+1;
			}
		}
	}
}

double runner = Matrix[0][0];
    
for(int i=0; i<Nbin; i++) 
{
	for(int j=0; j<Nbin; j++)
	{
        	if(Matrix[i][j] > runner)
		{ 
			runner = Matrix[i][j];	
			index_i = i;
			index_j = j;
		}
	}
}
    
Max_value = runner;
Max_sigma = lowlim_sigma+index_i*binsize;
Max_mi = lowlim_mi+index_j*binsize;

cout << Max_value << "   " << Max_sigma << "   " << Max_mi << endl;
*/

return 0;
}










void Input(void){
    ifstream ReadInput, ReadConf;
    
    //Read seed for random numbers
    int p1, p2;
    ifstream Primes("Primes");
    Primes >> p1 >> p2;
    Primes.close();

    ifstream input("seed.in");
    input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd.SetRandom(seed,p1,p2);
    input.close();
  
    //Read input informations
    ReadInput.open("input.dat");

    ReadInput >> temp;
    ReadInput >> niter;
    temp = temp/(1.*niter);
    beta = 1.0/temp;
    cout << "Temperature = " << temp << endl;
    cout << "Iteration index = " << niter << endl;

    ReadInput >> vol;

    box = pow(vol,1.0/2.0);
    cout << "Volume of the simulation box = " << vol << endl;
    cout << "Edge of the simulation box = " << box << endl;

    ReadInput >> delta;
    ReadInput >> nblk;
    ReadInput >> nstep;

    cout << "Metropolis delta = " << delta << endl;
    cout << "Number of blocks = " << nblk << endl;
    cout << "Number of steps in one block = " << nstep << endl << endl;

    ReadInput >> restart_mode;
    if(restart_mode==0) cout << "Restart mode: from scratch " << endl;
    else cout << "Restart mode: from previous run " << endl;

    //Prepare arrays for measurements
    iv = 0; //Potential energy
    is = 1;	
    im = 2;
 
    n_props = 3; //Number of observables

    //Read initial configuration
    if(restart_mode==0)				
    {
        ReadInput >> sigma >> mi;
        sigma = Pbc(sigma);/////////////////////////////////////////////////////////////////
	cout << "Sigma = " << sigma << endl;
        mi = Pbc(mi);/////////////////////////////////////////////////////////////////
	cout << "Mi = " << mi << endl;
    }

    else				
    {
        cout << "Read initial configuration from file config.final " << endl << endl;
        ReadConf.open("config.final");

        ReadConf >> sigma >> mi;
        sigma = Pbc(sigma);/////////////////////////////////////////////////////////////////
        mi = Pbc(mi);/////////////////////////////////////////////////////////////////
        
	ReadConf.close();
    }

    ReadInput.close();
  
    //Evaluate energy of the initial configuration
    Measure();

    //Print initial values for the energy
    cout << "Initial 'energy'  = " << walker[iv] << endl;
    cout << "----------------------------" << endl << endl;  
}










void Move(void)
{
	double p, energy_old, energy_new;
	double sigma_old, mi_old, sigma_new, mi_new;

	//Old
	sigma_old = sigma;
	mi_old = mi;
	energy_old = Eval_energy(sigma_old, mi_old);
	//cout << "Sigma old " << sigma_old << endl;
	//cout << "Mi old " << sigma_old << endl;
	//cout << "Energy old " << energy_old << endl;

	//New
	sigma_new = Pbc( sigma_old + delta*(rnd.Rannyu() - 0.5) );
	mi_new = Pbc( mi_old + delta*(rnd.Rannyu() - 0.5) );        
	energy_new = Eval_energy(sigma_new, mi_new);
	//cout << "Sigma new " << sigma_new << endl;
	//cout << "Mi new " << sigma_new << endl;
	//cout << "Energy new " << energy_new << endl;
        
	//Metropolis test
	p = exp(beta*(energy_old-energy_new));

	if(p >= rnd.Rannyu())
	{
        	//Update
		sigma = sigma_new;
		mi = mi_new;
            	accepted = accepted + 1.0;
		//cout << "ACCEPTED" << endl << endl;
	}

	//else cout << "REJECTED" << endl << endl;

        attempted = attempted + 1.0;
}










double Eval_energy(double sigma, double mi)
{
	double coef_8 = 3*(1+exp(-(mi*mi)/(sigma*sigma)));
	double coef_6 = 12*mi*mi - (5./3.)*coef_8;
	double coef_4 = 2*mi*mi*(2*mi*mi-5);
	double coef_2 = coef_8/3;
	double coef_0 = -2*mi*mi*exp(-(mi*mi)/(sigma*sigma));
	double den = 4*pow(sigma,4)*coef_2;
	double value = 0;

	value = (coef_8*pow(sigma,8) + coef_6*pow(sigma,6) + coef_4*pow(sigma,4) + coef_2*pow(sigma,2) + coef_0)/den;
 
	return value;
}










void Measure()
{
	walker[iv] = Eval_energy(sigma, mi);
	walker[is] = sigma;
	walker[im] = mi;
}










void Reset(int iblk){ //Reset block averages
    if(iblk == 1){
        for(int i=0; i<n_props; ++i){
            glob_av[i] = 0;
            glob_av2[i] = 0;
        }
    }
    
    for(int i=0; i<n_props; ++i){
        blk_av[i] = 0;
    }
    blk_norm = 0;
    attempted = 0;
    accepted = 0;
}










void Accumulate(void){      //Update block averages
    for(int i=0; i<n_props; ++i){
        blk_av[i] = blk_av[i] + walker[i];
    }
    blk_norm = blk_norm + 1.0;
}










void Averages(int iblk){    //Print results for current block

    ofstream Ene;
    //const int wd = 20;
    
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl;
      
    //Potential energy per particle 
    Ene.open("output.epot_iter"+to_string(niter),ios::app);
    stima_u = blk_av[iv]/blk_norm;  	
    glob_av[iv] += stima_u;
    glob_av2[iv] += stima_u*stima_u;
    err_u = Error(glob_av[iv],glob_av2[iv],iblk);
    Ene << iblk <<  " " << stima_u << " " << glob_av[iv]/(double)iblk << " " << err_u << endl;
    Ene.close();
    
    cout << "----------------------------" << endl << endl;    
}










void ConfFinal(void){
    ofstream WriteConf, WriteSeed;
    
    cout << "Print final configuration to file config.final " << endl << endl;
    WriteConf.open("config.final");
    WriteConf << sigma << "   " <<  mi << "   " << endl;
    WriteConf.close();
    
    rnd.SaveSeed();
}

















double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    //return r - box * rint(r/box);
    return r - box * rint((r/box)-0.5);	//for interval centred on box/2, avoids negative values for sigma and mi
}










double Error(double sum, double sum2, int iblk){
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}





/*
void Maximum(vector <double> value) 							
{
	double max = value[0];
	unsigned int N = value.size();
    
    	for (unsigned int i=0; i<N; i++) 
	{
        	if(value[i] > max)
		{ 
	  		max = value[i];	
			index_i = i;
			//index_j = j;
		}
    	}
    
	Max_value = max;
	Max_sigma = index_i;
	//Max_mi = index_j;    
}*/

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
