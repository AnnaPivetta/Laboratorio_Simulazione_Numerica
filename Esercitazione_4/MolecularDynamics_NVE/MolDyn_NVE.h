/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include "random.h"
Random rnd;
int restart;
double fs;
int nblocks;
void Accumula();
void Blocco(int b);
double ave_pot, ave_kin, ave_temp, ave_etot,ave_pressione;
double pressione;
int L;
int elemento;
double sigma;
double eKb;

//parameters, observables
const int m_props=4;
int n_props;

const int nbins=100;
double bin_size, blk_norm;
double blk_av[nbins], glob_av[nbins], glob_av2[nbins], walker[nbins];

int iv,ik,it,ie;
double stima_pot, stima_kin, stima_etot, stima_temp, stima_pressione, stima_gdir;
double media_pot, media_kin, media_etot, media_temp, media_pressione;
double media_pot2, media_kin2, media_etot2, media_temp2, media_pressione2;

// averages
double acc,att;
double temp_b;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part],vp[m_part],vs[m_part],va[m_part];  //velocità precedente, successiva, attuale

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

// simulation
int nstep, iprint, seed;
double delta;

//functions
void Input(void);
void Move(void);
void ConfFinal(void);
void ConfStart(char* filename);
void ConfXYZ(int);
void Measure(void);
double Force(int, int);
double Pbc(double);
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
