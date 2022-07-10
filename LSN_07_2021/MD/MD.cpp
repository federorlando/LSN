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
#include "MD.h"

using namespace std;

int main(){

int nconf = 1;
Input();    //Inizialization

for(int iblk=1; iblk <= nblk; ++iblk)
{  
	Reset(iblk);    //Reset block averages

	for(int istep=1; istep <= nstep; ++istep)
	{
		Move();
		Measure();
		Accumulate();   //Update block averages
		inst_val_energy[(iblk-1)*nstep+istep] = (walker[iv]/npart) + vtail;

		if(nconf==(nstep*nblk-1))
		{
			ConfFinal();
			//system("mv config.final config.old");   
			int systemRet = system("mv config.final config.old");	//Restart from penultimate configuration
			if(systemRet == -1) cerr << "The system method failed" << endl; 
        	}

		nconf += 1;
	}

	Averages(iblk);     //Print results for current block
}

ConfFinal();    //Write final configuration

ofstream WriteInstEne;
WriteInstEne.open("inst_energy.dat");
for (int i=0; i<nstep*nblk; ++i) WriteInstEne << inst_val_energy[i] << endl;
WriteInstEne.close();

return 0;
}








/************************************************************************************************/

void Input(void){
    ifstream ReadInput, ReadConf, ReadConfOld;
    
    cout << "Classic Lennard-Jones fluid        " << endl;
    cout << "Molecular dynamics            " << endl << endl;
    cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
    cout << "The program uses Lennard-Jones units " << endl;

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
    cout << "Initial temperature = " << temp << endl;

    ReadInput >> npart;
    cout << "Number of particles = " << npart << endl;

    ReadInput >> rho;
    cout << "Density of particles = " << rho << endl;
    vol = (double)npart/rho;
    box = pow(vol,1.0/3.0);
    cout << "Volume of the simulation box = " << vol << endl;
    cout << "Edge of the simulation box = " << box << endl;

    ReadInput >> rcut;
    cout << "Cutoff of the interatomic potential = " << rcut << endl << endl;
    
    //Tail corrections for potential energy and pressure
    vtail = (8.0*pi*rho)/(9.0*pow(rcut,9)) - (8.0*pi*rho)/(3.0*pow(rcut,3));
    ptail = (32.0*pi*rho)/(9.0*pow(rcut,9)) - (16.0*pi*rho)/(3.0*pow(rcut,3));
    cout << "Tail correction for the potential energy = " << vtail << endl;
    cout << "Tail correction for the virial = " << ptail << endl;

    ReadInput >> delta;
    ReadInput >> nblk;
    ReadInput >> nstep;

    cout << "The program integrates Newton equations with the Verlet method " << endl;
    cout << "Time step = " << delta << endl;
    cout << "Number of blocks = " << nblk << endl;
    cout << "Number of steps in one block = " << nstep << endl;
    cout << "Number of steps overall = " << nstep*nblk << endl << endl;

    ReadInput >> restart_mode;
    if(restart_mode==0) cout << "Restart mode: from initial config " << endl;
    else cout << "Restart mode: from previous run " << endl;
    ReadInput.close();

    //Prepare arrays for measurements
    iv = 0; //Potential energy
    iw = 1; //Virial
 
    n_props = 2; //Number of observables

    //Prepare measurement of g(r)
    igofr = 2;
    nbins = 100;
    n_props = n_props + nbins;		//first two components of walker are energy and virial, the rest is the histogram
    bin_size = (box/2.0)/(double)nbins;

    double sumv2 = 0.0, fs;



/////////////////////////////////////////////////////////////////////////////////////////////////////

if(restart_mode==false) //at the end one has I) set of coordinates x,y,z[i] read from config.0 II) set of random velocities vx,vy,vz[i] III) set of extrapolated old positions xold,yold,zold[i]
{
	//Read initial configuration
	cout << "Read initial configuration from file config.0 " << endl << endl;
	ReadConf.open("config.0");
	for (int i=0; i<npart; ++i)
	{
		ReadConf >> x[i] >> y[i] >> z[i];
		x[i] = x[i] * box;
		y[i] = y[i] * box;
		z[i] = z[i] * box;
	}
	ReadConf.close();

	//Prepare initial velocities
	cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
	double sumv[3] = {0.0, 0.0, 0.0};
	for (int i=0; i<npart; ++i)
	{
		vx[i] = rand()/double(RAND_MAX) - 0.5;
		vy[i] = rand()/double(RAND_MAX) - 0.5;
		vz[i] = rand()/double(RAND_MAX) - 0.5;

		sumv[0] += vx[i];
		sumv[1] += vy[i];
		sumv[2] += vz[i];
	}

	for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;

	for (int i=0; i<npart; ++i)
	{
		vx[i] = vx[i] - sumv[0];
		vy[i] = vy[i] - sumv[1];
		vz[i] = vz[i] - sumv[2];

		sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
	}

	sumv2 /= (double)npart;

	fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
	for (int i=0; i<npart; ++i)
	{
		vx[i] *= fs;
		vy[i] *= fs;
		vz[i] *= fs;

		xold[i] = Pbc(x[i] - vx[i] * delta);
		yold[i] = Pbc(y[i] - vy[i] * delta);
		zold[i] = Pbc(z[i] - vz[i] * delta);
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

if(restart_mode==true)
{           
	double a[npart], b[npart], c[npart];
        
	//Rename last two config files, to be read in following run as 'old.0' and 'old.final'
	system("rm -f old.0");
	system("rm -f old.final");
	system("mv config.old old.final");	//config.old was penultimate config from previous run
	system("mv config.final old.0");	//config.final was last config from previous run	
        
	//Read configuration 'old.0' (formerly 'config.final', which corresponded to nstep)
	cout << "Read actual configuration from file old.0 " << endl;
	ReadConf.open("old.0");
	for (int i=0; i<npart; ++i)
	{
		ReadConf >> x[i] >> y[i] >> z[i];
		x[i] = x[i] * box;
		y[i] = y[i] * box;
		z[i] = z[i] * box;
        }
	ReadConf.close();
        
	//Read configuration 'old.final' (formerly 'config.old', which corresponded to nstep-1)
	cout << "Read old configuration from file old.final " << endl << endl;
	ReadConfOld.open("old.final");
	for (int i=0; i<npart; ++i)
	{
		ReadConfOld >> xold[i] >> yold[i] >> zold[i];
		xold[i] = xold[i] * box;
		yold[i] = yold[i] * box;
		zold[i] = zold[i] * box;
            
		a[i]=xold[i];
		b[i]=yold[i];
		c[i]=zold[i];
	}
	ReadConfOld.close();
        
	Move();     //One step of the Verlet algorithm
        
	for (int i=0; i<npart; ++i)		//compute velocities NOT at t (already done by Move()!) but at t+(dt/2)
	{
		vx_halfdt[i]=Pbc((x[i]-xold[i]))/(1.*delta);
		vy_halfdt[i]=Pbc((y[i]-yold[i]))/(1.*delta);
		vz_halfdt[i]=Pbc((z[i]-zold[i]))/(1.*delta);
            
		sumv2 += vx_halfdt[i]*vx_halfdt[i] + vy_halfdt[i]*vy_halfdt[i] + vz_halfdt[i]*vz_halfdt[i];
        }
        
	sumv2 /= (double)(npart);
        
	fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor
    
	for (int i=0; i<npart; ++i)
	{
		vx[i] *= fs;		//scale VELOCITIES AT t
		vy[i] *= fs;
		vz[i] *= fs;

		xold[i] = Pbc(x[i] - vx[i] * delta);	//new estimate of r(t), to be contained in xold for routine Move to start correctly
		yold[i] = Pbc(y[i] - vy[i] * delta);	//x,y,z[i], which contains r(t+dt), must not be touched any longer!
		zold[i] = Pbc(z[i] - vz[i] * delta);
        }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
  
    //Evaluate potential energy and virial of the initial configuration
    Measure();

    //Print initial values for the potential energy and virial
    cout << "Initial potential energy (with tail corrections) = " << walker[iv]/(double)npart + vtail << endl;
    cout << "Virial (with tail corrections) = " << walker[iw]/(double)npart + ptail << endl;
    cout << "Pressure (with tail corrections) = " << rho * temp + (walker[iw] + (double)npart * ptail) / vol << endl << endl;
}






/************************************************************************************************/

void Move(void)
{
	double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

	for(int i=0; i<npart; ++i)
	{
		fx[i] = Force(i,0);
		fy[i] = Force(i,1);
		fz[i] = Force(i,2);
	}

	for(int i=0; i<npart; ++i)
	{
		xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
		ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
		znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

		vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
		vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
		vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

		xold[i] = x[i];
		yold[i] = y[i];
		zold[i] = z[i];

		x[i] = xnew;
		y[i] = ynew;
		z[i] = znew;
	}

	return;
}





/************************************************************************************************/

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}





/************************************************************************************************/

void Measure()
{
	int bin=0;
	double v = 0.0, w = 0.0;
	double vij, wij;
	double dx, dy, dz, dr;
    
	//reset the hystogram of g(r)
	for (int k=igofr; k<igofr+nbins; ++k) walker[k]=0.0;
    
	//cycle over pairs of particles
	for (int i=0; i<npart-1; ++i)
	{
        	for (int j=i+1; j<npart; ++j)
		{
			// distance i-j in pbc
			dx = Pbc(x[i] - x[j]);
			dy = Pbc(y[i] - y[j]);
			dz = Pbc(z[i] - z[j]);
            
			dr = dx*dx + dy*dy + dz*dz;
			dr = sqrt(dr);
            
			//update of the histogram of g(r)
			for(int k=0; k<nbins; k++)
			{
        	        	if(dr >= k*bin_size && dr < (k+1)*bin_size)
				{
					bin = k;                //Select bin to act upon
					walker[bin+2] += 2.;    //The bin content is increased by 2
					break;                  //To quickly exit the for cycle on bins
        	        	}
			}
            
			if(dr < rcut)
			{
				vij = 1.0/pow(dr,12) - 1.0/pow(dr,6);
        	        	wij = 1.0/pow(dr,12) - 0.5/pow(dr,6);
        	        	// contribution to energy and virial
        	        	v += vij;
        	        	w += wij;
			}
        	}
	}
    
	walker[iv] = 4.0 * v;
	walker[iw] = 48.0 * w / 3.0;
    
	//Histogram normalized inside the block average
}





/************************************************************************************************/

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
}





/************************************************************************************************/

void Accumulate(void){      //Update block averages
    for(int i=0; i<n_props; ++i){
        blk_av[i] = blk_av[i] + walker[i];
    }
    blk_norm = blk_norm + 1.0;
}





/************************************************************************************************/

void Averages(int iblk){    //Print results for current block

    double V, gdir;
    ofstream Gofr, Gave, Ene, Pres;
    const int wd = 20;
    
    cout << "Block number " << iblk << endl;
      
    //Potential energy per particle 
    Ene.open("output.epot.0",ios::app);
    stima_u = blk_av[iv]/blk_norm/(double)npart + vtail;  //Potential energy
    glob_av[iv] += stima_u;
    glob_av2[iv] += stima_u*stima_u;
    err_u = Error(glob_av[iv],glob_av2[iv],iblk);
    Ene << iblk <<  " " << stima_u << " " << glob_av[iv]/(double)iblk << " " << err_u << endl;
    Ene.close();

    //Pressure
    Pres.open("output.pres.0",ios::app);
    stima_pres = rho * temp + (blk_av[iw]/blk_norm + ptail*(double)npart)/ vol;  //Pressure
    glob_av[iw] += stima_pres;
    glob_av2[iw] += stima_pres*stima_pres;
    err_pres = Error(glob_av[iw],glob_av2[iw],iblk);
    Pres << iblk <<  " " << stima_pres << " " << glob_av[iw]/(double)iblk << " " << err_pres << endl;
    Pres.close();
    
    //In output.gofr.0: first column is block number and the other 100 columns are values of g(r=i*bin_size), i=1...100
    Gofr.open("output.gofr.0",ios::app);  
    Gofr << iblk << " ";
    for(int i=0; i<nbins; i++)
    {
        V = rho*npart*(4/3)*pi*(pow((i+1)*bin_size,3)-pow(i*bin_size,3));	//Normalization factor
        gdir = (blk_av[igofr+i]/blk_norm)/V;                            	//g(r=i*bin_size) for the i-block
        Gofr << " "  << gdir << " ";                                	//g(r) -> Gofr
        glob_av[igofr+i] += gdir;                                       	//To be used for the total average
        glob_av2[igofr+i] += pow(gdir,2);
    }
    Gofr << endl;
    Gofr.close();

    //g(r) -> Gave
    Gave.open("output.gave.0",ios::app);   
    if (Gave.is_open()) cerr << "Everything OK" << endl;
    else cerr << "PROBLEM: Unable to open Gave" << endl; 
    if(iblk==nblk)
    {
	for(int i=igofr;i<igofr+nbins;i++)
	{
		err_gdir = Error(glob_av[i],glob_av2[i],iblk);
		//We print the value of r and the corresponding average value of g(r) on the total number of blocks used
		Gave << (i-igofr)*bin_size << " "  << glob_av[i]/(double)iblk << " "  << err_gdir << endl;
        }
    }

    Gave.close();
    
    cout << "----------------------------" << endl << endl;    
}





/************************************************************************************************/

void ConfFinal(void){
    ofstream WriteConf, WriteSeed;
    
    cout << "Print final configuration to file config.final " << endl << endl;
    WriteConf.open("config.final");
    for (int i=0; i<npart; ++i){
        WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
    }
    WriteConf.close();
    
    rnd.SaveSeed();
}





/************************************************************************************************/

void ConfXYZ(int nconf){ //Write configuration in .xyz format
    ofstream WriteXYZ;
    
    WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
    WriteXYZ << npart << endl;
    WriteXYZ << "This is only a comment!" << endl;
    for (int i=0; i<npart; ++i){
        WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
    }
    WriteXYZ.close();
}





/************************************************************************************************/

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}





/************************************************************************************************/

double Error(double sum, double sum2, int iblk){
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
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
