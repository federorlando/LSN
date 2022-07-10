/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
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
#include <string>
#include <sstream>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main(){ 

Input(); //Inizialization

stringstream stream;
stream.precision(1);
stream << fixed;
stream << temp;
string str = stream.str();
string command = "mkdir -p ";
namedir = "T_" + str;    
command = command + namedir;
system(command.c_str());

for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
{
	Reset(iblk);   //Reset block averages

	for(int istep=1; istep <= nstep; ++istep)
	{
		Move(metro);
		Measure();
		Accumulate(); //Update block averages
		inst_val_energy[(iblk-1)*nstep+istep] = walker[iu]/(double)nspin;
		inst_val_mag[(iblk-1)*nstep+istep] = walker[im]/(double)nspin;
	}

    	Averages(iblk);   //Print results for current block
}

ConfFinal(); //Write final configuration

ofstream WriteInstEne;
WriteInstEne.open("inst_energy_T"+str+"altseed");
for (int i=0; i<nstep*nblk; ++i) WriteInstEne << i << "        " << inst_val_energy[i] << endl;
WriteInstEne.close();
ofstream WriteInstMag;
WriteInstMag.open("inst_mag_T"+str+"altseed");
for (int i=0; i<nstep*nblk; ++i) WriteInstMag << i << "        " << inst_val_mag[i] << endl;
WriteInstMag.close();

Measure();
cout << "Final energy (per DOF) = " << walker[iu]/(double)nspin << endl;
cout << "Final magnetic moment (per DOF) = " << walker[im]/(double)nspin << endl << endl;
cout << "***********************************" << endl << endl;

return 0;
}










void Input(void)
{
  ifstream ReadInput;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs
  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;

  ReadInput >> nblk;

  ReadInput >> nstep;

  ReadInput >> restart_mode;
  if(restart_mode==0) cout << "Restart mode: from scratch " << endl;
  else cout << "Restart mode: from previous run " << endl;

  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl;
  cout << "Overall steps = " << nblk*nstep << endl << endl;
  ReadInput.close();


//Prepare arrays for measurements
  iu = 0; //Energy
  iusquared = 1; //Squared energy for heat capacity
  im = 2; //Magnetization
  imsquared = 3; //Squared magnetization for susceptibility
 
  n_props = 4; //Number of observables

//initial configuration T=infinity (disorder)

if(restart_mode==0)				
{
	for (int i=0; i<nspin; ++i)
	{
	if(rnd.Rannyu() >= 0.5) s[i] = 1;
	else s[i] = -1;
	//cout << s[i] << endl;
	}
}

else 		
{
	ifstream ImportConf;
	ImportConf.open("config.final");
	if (ImportConf.fail()) 
	{
		cerr << "****************************************************" << endl;
		cerr << "ERROR: RESTART FROM PREVIOUS RUN, BUT FILE config.final ABSENT" << endl;
		cerr << "****************************************************" << endl;
		return;
	}
	cout << "Initial config" << endl;
	for (int i=0; i<nspin; ++i) 
	{
		ImportConf >> s[i];
		cout << s[i] << "   ";
	}
	ImportConf.close();
	cout << endl;
}

//initial configuration T=0 (order)
  //for (int i=0; i<nspin; ++i) s[i] = 1;	//all up
  //for (int i=0; i<nspin; ++i) s[i] = -1;	//all down

//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial energy (per DOF) = " << walker[iu]/(double)nspin << endl;
  cout << "Initial magnetic moment (per DOF) = " << walker[im]/(double)nspin << endl << endl;
}










void Move(int metro)
{
	int o;
	double energy_old, energy_new, sm, alpha;
	double energy_up, energy_down;

	for(int i=0; i<nspin; ++i)
	{
		//Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
		o = (int)(rnd.Rannyu()*nspin);

		if(metro==1) 	//Metropolis
  		{
			attempted++;

			sm = s[o];
			energy_old = Boltzmann(sm, o);		
			s[o] = -s[o]; 				
			sm = s[o];
			energy_new = Boltzmann(sm, o);
			
			alpha = exp(-(energy_new-energy_old)*beta);
			if(alpha<1 && rnd.Rannyu()>alpha) s[o] = -s[o];
			else accepted++; 
		}

		else 		//Gibbs 
		{
			energy_up = Boltzmann(1,o);
			energy_down = Boltzmann(-1,o);
			alpha = 1/(1+exp(-beta*(energy_down-energy_up)));
			if(rnd.Rannyu()<alpha) s[o]=1;
			else s[o]=-1;
 		}
  	}
}












double Boltzmann(int sm, int ip)
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}











void Measure()
{
	double u = 0.0, m = 0.0, u_squared = 0.0, m_squared = 0.0;

	//cycle over spins
	for (int i=0; i<nspin; ++i)
	{
		u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
		m += s[i];
	}

	u_squared = pow(u,2);		//will be useful to compute C_V
	m_squared = pow(m,2);		//will be useful to compute susceptibility
    
	walker[iu] = u;
	walker[iusquared] = u_squared;
	walker[im] = m;
	walker[imsquared] = m_squared;
}










void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}










void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}










void Averages(int iblk) //Print results for current block
{
    
	ofstream Ene, Heat, Mag, Chi;
	const int wd=20;
    
	cout << "Block number " << iblk << endl;
	if(metro==1) cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    	if(h==0.0) Ene.open(namedir+"/output.ene."+to_string(metro),ios::app);		
	else Ene.open(namedir+"/output.ene.h_ext."+to_string(metro),ios::app);
	//Ene.open("output.ene.0",ios::app);
	stima_u = blk_av[iu]/blk_norm/(double)nspin; 								//Internal energy per DOF
	glob_av[iu]  += stima_u;
	glob_av2[iu] += stima_u*stima_u;
	err_u = Error(glob_av[iu],glob_av2[iu],iblk);
	Ene << setw(wd) << iblk << setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
	Ene.close();

	if(h==0.0)Heat.open(namedir+"/output.cv."+to_string(metro),ios::app);
	else Heat.open(namedir+"/output.cv.h_ext."+to_string(metro),ios::app);
	//Heat.open("output.cv.0",ios::app);
	stima_c = pow(beta,2)*(blk_av[iusquared]/blk_norm - pow(blk_av[iu]/blk_norm,2))/(double) nspin;  	//Heat capacity per DOF (specific heat) = K_B*(beta)^2*[<H^2>-<H>^2] 
    	glob_av[iusquared]  += stima_c;
    	glob_av2[iusquared] += stima_c*stima_c;
   	err_c = Error(glob_av[iusquared],glob_av2[iusquared],iblk);
	Heat << setw(wd) << iblk << setw(wd) << stima_c << setw(wd) << glob_av[iusquared]/(double)iblk << setw(wd) << err_c << endl;
	Heat.close();

	if(h==0.0) Mag.open(namedir+"/output.mag."+to_string(metro),ios::app);
	else Mag.open(namedir+"/output.mag.h_ext."+to_string(metro),ios::app);
	//Mag.open("output.mag.0",ios::app);
	stima_m = blk_av[im]/blk_norm/(double)nspin; 								//Magnetization per DOF (average magnetic moment)
	glob_av[im]  += stima_m;
	glob_av2[im] += stima_m*stima_m;
	err_m = Error(glob_av[im],glob_av2[im],iblk);
	Mag << setw(wd) << iblk << setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
	Mag.close();

	if(h==0.0) Chi.open(namedir+"/output.chi."+to_string(metro),ios::app);
	else Chi.open(namedir+"/output.chi.h_ext."+to_string(metro),ios::app);
	//Chi.open("output.chi.0",ios::app);
	stima_x = beta*blk_av[imsquared]/blk_norm/(double) nspin;    						//Magnetic susceptibility = beta*[<M^2>-<M>^2] 
	glob_av[imsquared]  += stima_x;
	glob_av2[imsquared] += stima_x*stima_x;
	err_x = Error(glob_av[imsquared],glob_av2[imsquared],iblk);
	Chi << setw(wd) << iblk << setw(wd) << stima_x << setw(wd) << glob_av[imsquared]/(double)iblk << setw(wd) << err_x << endl;
	Chi.close();

	cout << "----------------------------" << endl;
}










void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}










int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}










double Error(double sum, double sum2, int iblk)
{
    if(iblk==1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
