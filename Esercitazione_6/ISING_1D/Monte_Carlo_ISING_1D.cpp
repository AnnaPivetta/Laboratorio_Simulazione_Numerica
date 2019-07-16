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
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main()
{ conta=0;
  Input(); //Inizialization
  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move(metro);
      Measure();
      Accumulate(); //Update block averages
    }
    Averages(iblk);   //Print results for current block
  }
  ConfFinal(); //Write final configuration

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
 
  
//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> restart;

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

  ReadInput >> nblk;

  ReadInput >> nstep;   //numero di stap PER BLOCCO	
  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();

	cout<<"RESTART= "<<restart<<endl;
  if (restart==0){
	int p1, p2;
	ifstream Primes("Primes");
	Primes >> p1 >> p2 ;
	Primes.close();

	ifstream input("seed.in");
	input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
	rnd.SetRandom(seed,p1,p2);
	input.close();
    }

   if (restart==1){
	int p1, p2;
	ifstream Primes("Primes");
	Primes >> p1 >> p2 ;
	Primes.close();

	ifstream input("seed.in");
	input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
	rnd.SetRandom(seed,p1,p2);
	input.close();
    }
//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables

//initial configuration

  if(restart==0){
	  system ("rm output.heat.0 output.ene.0 output.mag.0 output.chi.0");
	  for (int i=0; i<nspin; ++i)
	  {
	    if(rnd.Rannyu() >= 0.5) s[i] = 1;
	    else s[i] = -1;
	  }
  }
  if (restart==1){
	cout<<"RESTART VALE 1!"<<endl;
	ifstream in;
	in.open("config.final");  
	if(in.fail()){
    		cout << endl << "Errore apertura file config.final, si riparte come da restart=0" << endl;
		for (int i=0; i<nspin; ++i){	
		    if(rnd.Rannyu() >= 0.5) {s[i] = 1;}
		    else {s[i] = -1;}
	  	}
 	}    
	for (int i=0; i<nspin; ++i){
  	  in >> s[i];
  	}

  }

//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}


void Move(int metro)
{
  int o;
  double energy_old, energy_new;
  double delta_energy;

  for(int i=0; i<nspin; ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);

    if(metro==1) //Metropolis
    {	attempted++;
 	energy_old=Boltzmann(s[o],o);
	energy_new=Boltzmann(-s[o],o);	
	delta_energy=energy_new - energy_old;
	double alpha=min(1.0,exp(-beta*delta_energy));
	double p=rnd.Rannyu(0,1);
	if (p<=alpha) {s[o]=-s[o];
			accepted++;
			}
	
    }
    else //Gibbs sampling
    {
	double delta_E=Boltzmann(-1,o)-Boltzmann(1,o);
	double p_up= 1.0/(1+exp(-(delta_E*beta)));
	double p=rnd.Rannyu(0,1);
	if (p<=p_up) {s[o]=+1;}
	if (p>p_up) {s[o]=-1;}
    }
  }
}

double Boltzmann(int sm, int ip) //ip contraddistingue lo spin che si è mosso
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}

void Measure()
{
 
  double u = 0.0; //energia e magnetizzazione
  double u2 = 0.0;
  double spin= 0.0;
  

//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
     u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]); //energia istantanea della configurazione
     u2+= (-J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)])) * (-J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]));
     spin+=s[i];
// INCLUDE YOUR CODE HERE
     
  }
  walker[iu] = u;
  walker[ic] = u*u;           //nb salvo solo il quadrato dell'energia!
//  if (h==0) {walker[im]=0.;}
  walker[im]=spin;
  walker[ix]=walker[im]*walker[im];  //nb ho salvato solo il quadrato della magnetizzazione;
// INCLUDE YOUR CODE HERE
}


void Reset(int iblk) //Reset block averages
{
   if(iblk==nblk) {
	cout<<"ultimo blocco= "<<glob_av[iu]<<endl;}
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

   for(int i=0; i<n_props; ++i)   //ciclo sulle prprietà che devo calcolare
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{
    
   ofstream Ene, Heat, Mag, Chi;
 
    
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << (double)(accepted/(double)(attempted)) << endl << endl;
//cout <<"accettati= "<<accepted<<endl;
//cout<<"totali= "<<attempted<<endl;
    
    Ene.open("output.ene.0",ios::app);
    stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy stima del valo medio nel blocco
    glob_av[iu]  += stima_u;
    glob_av2[iu] += stima_u*stima_u;
    err_u=Error(glob_av[iu],glob_av2[iu],iblk);
    Ene<<iblk<<"	"<<glob_av[iu]/(double)iblk<<"	"<<err_u<<endl;
    //Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
    Ene.close();

    Heat.open("output.heat.0",ios::app);
    stima_c=beta*beta*( (blk_av[ic]/blk_norm) - nspin*nspin*stima_u*stima_u )/((double)(nspin));
    glob_av[ic]  += stima_c;
    glob_av2[ic] += stima_c*stima_c;
    err_c=Error(glob_av[ic],glob_av2[ic],iblk);
    //Heat << setw(wd) << iblk <<  setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
    Heat<<iblk<<"	"<<glob_av[ic]/(double)iblk<<"	"<<err_c<<endl;
    Heat.close();

    Mag.open("output.mag.0",ios::app);
    stima_m=blk_av[im]/blk_norm/(double)nspin;
    glob_av[im]  += stima_m;
    glob_av2[im] += stima_m*stima_m;
    err_m=Error(glob_av[im],glob_av2[im],iblk);
    //Mag << setw(wd) << iblk <<  setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
    Mag<<iblk<<"	"<<glob_av[im]/(double)iblk<<"	"<<err_m<<endl;
    Mag.close();

    Chi.open("output.chi.0",ios::app);
    stima_x=beta*(blk_av[ix]/blk_norm)/((double)(nspin));
    glob_av[ix]  += stima_x;
    glob_av2[ix] += stima_x*stima_x;
    err_x=Error(glob_av[ix],glob_av2[ix],iblk);
    //Chi << setw(wd) << iblk <<  setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
    Chi<<iblk<<"	"<<glob_av[ix]/(double)iblk<<"	"<<err_x<<endl;
    Chi.close();

// INCLUDE YOUR CODE HERE

    cout << "----------------------------" << endl << endl;
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
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk); //non c'è /sqrt(N-1) perchè N sarà abbastanza grande
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
