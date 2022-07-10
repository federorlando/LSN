/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __fluid_
#define __fluid_

//Random numbers
#include "random.h"
int seed[4];
Random rnd;

//parameters, observables
const int m_props=1000;
int n_props, iv, iw, igofr;
double vtail,ptail,bin_size,nbins,sd;
double walker[m_props];
const int inst_val_dim=10000000;
double inst_val_energy[inst_val_dim];

// averages
double blk_av[m_props],blk_norm;
double glob_av[m_props],glob_av2[m_props];
double stima_u,stima_pres,err_u,err_pres,err_gdir;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part];
double vx[m_part],vy[m_part],vz[m_part];
double vx_halfdt[m_part],vy_halfdt[m_part],vz_halfdt[m_part];
double xold[m_part],yold[m_part],zold[m_part];

// thermodynamical state
int npart;
double beta,temp,vol,rho,box,rcut;

// simulation
int nstep, nblk, restart_mode;
double delta;

//pigreco
const double pi=3.1415927;

//functions
void Input(void);
void Move(void);
double Force(int, int);
void Measure(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void ConfFinal(void);
void ConfXYZ(int);
double Pbc(double);
double Error(double,double,int);

#endif

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
