/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <vector>
#include <string>
#include "random.h"
#include "MolDyn_NVE.h"

using namespace std;

int main(){
 
Input();             			//Inizialization
int nconf = 1;
int n_blocks = 500;

for(int istep=1; istep <= nstep; ++istep)
{
	Move();           		//Move particles with Verlet algorithm
	if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
        Measure();     			//Properties measurement 

        if(nconf==(nstep-1))
	{
		ConfFinal();
		//system("mv config.final config.old");   
		int systemRet = system("mv config.final config.old");	//Restart from penultimate configuration
		if(systemRet == -1) cerr << "The system method failed" << endl; 
        }

        nconf += 1;
}

ConfFinal();        			//Write final configuration to restart
Blocking_average(n_blocks);

return 0;
}



/********************************************************************************************************************************/

void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf,ReadConfOld;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open("input.dat"); //Read input

  ReadInput >> temp;
  cout << "Initial temperature = " << temp << endl;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;
  ReadInput >> restart_mode;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl;

  if(restart_mode==false)  cout << "Restart from: initial configuration" << endl << endl;	//true = 1, false = 0
  else cout << endl << "Restart from: previous run" << endl;

  ReadInput.close();

  //Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables

  double sumv2 = 0.0, fs;


//////////////////////////////////////////////////////////////////////////////////////////////////////////////

if(restart_mode==false) //at the end one has I) set of coordinates x,y,z[i] read from config.0 II) set of random velocities vx,vy,vz[i] III) set of extrapolated old positions xold,yold,zold[i]
{
	////Rename last two config files, to be read in following run as 'old.0' and 'old.final'
	//system("rm -f old.0");
	//system("rm -f old.final");
	//system("mv config.old old.final");
	//system("mv config.final old.0");

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


//////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
            
		//a[i]=xold[i];
		//b[i]=yold[i];
		//c[i]=zold[i];
            
		//xold[i] = Pbc(x[i] - 2.0 * vx[i] * delta);
		//yold[i] = Pbc(y[i] - 2.0 * vy[i] * delta);
		//zold[i] = Pbc(z[i] - 2.0 * vz[i] * delta);
            
		//I valori di xold, yold e zold trovati dall'algoritmo di Verlet corrispondono esattamente ai valori di x, y e z caricati dal file old.0
		//x[i]=a[i];
		//y[i]=b[i];
		//z[i]=c[i];
        }
}





   return;
}



/********************************************************************************************************************************/

void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

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



/********************************************************************************************************************************/

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



/********************************************************************************************************************************/

void Measure(){ //Properties measurement
  //int bin;
  double v, t, w, vij, wij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp, Pres;

  Epot.open("output_epot.dat",ios::app);
  Ekin.open("output_ekin.dat",ios::app);
  Temp.open("output_temp.dat",ios::app);
  Etot.open("output_etot.dat",ios::app);
  Pres.open("output_pres.dat",ios::app);

  v = 0.0; //reset observables
  t = 0.0;
  w = 0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i)
  {
    for (int j=i+1; j<npart; ++j)
    {

     dx = Pbc( xold[i] - xold[j] ); // here I use old configurations [old = r(t)]
     dy = Pbc( yold[i] - yold[j] ); // to be compatible with EKin which uses v(t)
     dz = Pbc( zold[i] - zold[j] ); // => EPot should be computed with r(t)

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

     if(dr < rcut)
     {
	//Potential energy
       	vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
       	v += vij;

 	//Pressure
     	wij=48.0/pow(dr,12) - 24.0/pow(dr,6);
     	w += wij;
     }
    }          
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    stima_pot = v/(double)npart; //Potential energy per particle
    stima_kin = t/(double)npart; //Kinetic energy per particle
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_etot = (t+v)/(double)npart; //Total energy per particle
    stima_pres = rho*stima_temp + w/(double)(3*vol);  //Pressure

    Epot << stima_pot  << endl;
    Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;
    Pres << stima_pres << endl;

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();
    Pres.close();

    return;
}



/********************************************************************************************************************************/

void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
  return;
}

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

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}



/********************************************************************************************************************************/

void Blocking_average(int nblocks)
{
	vector <double> ve,vu,vk,vt,vp,sum_prog,err_prog;
    
	for(int i=0;i<nblocks;i++)
	{
		ve.push_back(0);
		vk.push_back(0);
		vu.push_back(0);
		vt.push_back(0);
		vp.push_back(0);
        
		sum_prog.push_back(0);
		err_prog.push_back(0);
    	}
    
	ifstream input_e, input_k, input_u, input_t, input_p;    
	input_e.open("output_etot.dat");
	input_k.open("output_ekin.dat");
	input_u.open("output_epot.dat");
	input_t.open("output_temp.dat");
	input_p.open("output_pres.dat");
    
	//Di seguito si memorizzano le varie misure di etot, ekin, epot e temp.
	double val = 0;
	int L = nstep/nblocks;
	for(int i=0; i<nblocks; i++)
	{
        	for(int j=0; j<L; j++)
		{
			input_e >> val;
			ve[i] += val;
			input_k >> val;
			vk[i] += val;
			input_u >> val;
			vu[i] += val;
			input_t >> val;
			vt[i] += val;
			input_p >> val;
			vp[i] += val;
		}
	ve[i] /= L;
	vk[i] /= L;
	vu[i] /= L;
	vt[i] /= L;
	vp[i] /= L;
    	}
    
	input_e.close();
	input_k.close();
	input_u.close();
	input_t.close();
	input_p.close();    
    
	ofstream output_e;
	ofstream output_k;
	ofstream output_u;
	ofstream output_t;
	ofstream output_p;
    
	output_e.open("ave_etot.dat");
	output_k.open("ave_ekin.dat");
	output_u.open("ave_epot.dat");
	output_t.open("ave_temp.dat");
	output_p.open("ave_pres.dat");
    
	//Total energy
	data_blocking(ve,sum_prog,err_prog);
    	for(int i=0; i<nblocks; i++)
		output_e << i << " " << sum_prog[i] << " " << err_prog[i] << endl;
    
	//Kinetic energy
	data_blocking(vk,sum_prog,err_prog);
    	for(int i=0; i<nblocks; i++)
		output_k << i << " " << sum_prog[i] << " " << err_prog[i] << endl;
    
	//Potential energy
	data_blocking(vu,sum_prog,err_prog);
    	for(int i=0; i<nblocks; i++)
		output_u << i << " " << sum_prog[i] << " " << err_prog[i] << endl;
    
	//Temperature
	data_blocking(vt,sum_prog,err_prog);
    	for(int i=0; i<nblocks; i++)
		output_t << i << " " << sum_prog[i] << " " << err_prog[i] << endl;

	//Pressure
	data_blocking(vp,sum_prog,err_prog);
    	for(int i=0; i<nblocks; i++)
		output_p << i << " " << sum_prog[i] << " " << err_prog[i] << endl;
       
	output_e.close();
	output_k.close();
	output_u.close();
	output_t.close();
	output_p.close();
    
	return;
}



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

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
