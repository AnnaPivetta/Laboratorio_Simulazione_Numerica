#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <string>  
#include "MolDyn_NVE.h"
#include "random.h"

/*
bisogna prima fare girare il codice mettendo restart=0 e nsteps=500,
poi restart=1 e nstep=500 e farlo girare finchè non si è soddisfatti della temperatura

*/
using namespace std;

int main(int argc, char *argv[]){ 
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

  Input();             //Inizialization
  int nconf = 1;
  L=nstep/nblocks;
  ave_pot=0;
  ave_kin=0;
  ave_temp=0;
  ave_etot=0;
  ave_pressione=0;
  media_pot=0;
  media_kin=0;
  media_etot=0;
  media_temp=0;
  media_pressione=0;
  media_pot2=0;
  media_kin2=0;
  media_etot2=0;
  media_temp2=0;
  media_pressione2=0;
  int block=1;
  system("rm ave_gdir.out");
  for(int i=0; i<nbins; ++i){
  	blk_av[i] = 0;
	walker[i] = 0;
  }
  blk_norm = 0;
  for(int istep=1; istep <= nstep; ++istep){
     Move();           //Move particles with Verlet algorithm
     if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
     if(istep%10 == 0){
        Measure();     //Properties measurement
//        ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
	Accumula();
        nconf += 1;
     }
     if (istep%L == 0) {
	Blocco(block);
        block++;
        cout<<"blocco: "<<block<<endl;
     }
  }
 
  ConfFinal();         //Write final configuration to restart
  cout<<"temperatura attuale: "<<stima_temp<<endl;

  return 0;
}


void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;
  //double ep, ek, pr, et, vir;
  int restart;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open("input.dat"); //Read input
  ReadInput >> restart;
  cout<<"restart="<<restart<<endl;

  ReadInput >> temp;

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
  ReadInput >> nblocks;
  ReadInput >> iprint;
  ReadInput >> elemento;

  //per il calcolo della g(r)
 
  bin_size = (double)(box/2.0)/((double)(nbins));
  for(int i=0; i<nbins; i++){
	glob_av[i]=0;
	glob_av2[i]=0;
  }

  if (elemento==1) { //Argon
	sigma=0.34*pow(10,-9);
	eKb=120;
  }
  if (elemento==2) { //Krypton
	sigma=0.364*pow(10,-9);
	eKb=164;
  }
  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables
  if(restart==0){
 	 system("rm ave_ekin.out ave_epot.out ave_temp.out ave_etot.out ave_pressione.out");	
	//Read initial configuration
	  cout << "Read initial configuration from file config.0 " << endl << endl;
	  ReadConf.open("config.0");
	  for (int i=0; i<npart; ++i){
	    ReadConf >> x[i] >> y[i] >> z[i];
	    x[i] = x[i] * box;
	    y[i] = y[i] * box;
	    z[i] = z[i] * box;
	  }
	  ReadConf.close();

	//Prepare initial velocities
	   cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
	   double sumv[3] = {0.0, 0.0, 0.0};
	   for (int i=0; i<npart; ++i){
	     vx[i] = rnd.Rannyu(0,1); 
	     vy[i] = rnd.Rannyu(0,1);
	     vz[i] = rnd.Rannyu(0,1);

	     sumv[0] += vx[i];
	     sumv[1] += vy[i];
	     sumv[2] += vz[i];
	   }
	   for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
	   double sumv2 = 0.0;
	   for (int i=0; i<npart; ++i){
	     vx[i] = vx[i] - sumv[0];
	     vy[i] = vy[i] - sumv[1];
	     vz[i] = vz[i] - sumv[2];

	     sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
	   }
	   sumv2 /= (double)npart;
	   fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
	   for (int i=0; i<npart; ++i){
	     vx[i] *= fs;
	     vy[i] *= fs;
	     vz[i] *= fs;

	     xold[i] = x[i] - vx[i] * delta;  //calcolo delle ipotetiche x vecchie
	     yold[i] = y[i] - vy[i] * delta;
	     zold[i] = z[i] - vz[i] * delta;
	   }
 }

 if(restart==1) {
	  ifstream in;
	  in.open("config.final");
	  for (int i=0; i<npart; ++i){
	    in >> x[i] >> y[i] >> z[i];
	    x[i] = x[i] * box;
	    y[i] = y[i] * box;
	    z[i] = z[i] * box;
	  }
	  in.close();
	  in.open("config.pre");
	  for (int i=0; i<npart; ++i){
	    in >> xold[i] >> yold[i] >> zold[i];
	    xold[i] = xold[i] * box;
	    yold[i] = yold[i] * box;
	    zold[i] = zold[i] * box;
	  }
	  in.close();
		
          Move();
	double vquadro=0;
	 for(int i=0; i<npart; ++i){
		vquadro+=(Pbc(x[i]-xold[i])*Pbc(x[i]-xold[i])/(delta*delta)+Pbc(y[i]-yold[i])*Pbc(y[i]-yold[i])/(delta*delta)+Pbc(z[i]-zold[i])*Pbc(z[i]-zold[i])/(delta*delta))/(double)npart;
	}	                          
	 double T=vquadro/3.;
	 cout<< "Temperatura con restart" << T << endl;
	 fs=sqrt(temp/T);
	 for(int i=0; i<npart; ++i){
	 vx[i]=vx[i]*fs;
	 vy[i]=vy[i]*fs;
	 vz[i]=vz[i]*fs;
	 }
	 for(int i=0; i<npart; ++i){
	   xold[i]=Pbc(x[i]-vx[i]*delta);
	   yold[i]=Pbc(y[i]-vy[i]*delta);
	   zold[i]=Pbc(z[i]-vz[i]*delta);
	 }
 } 
   return;
}

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

void Measure(){ //Properties measurement
  int bin;
  double v, t, vij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp, Press;

  Epot.open("output_epot.dat",ios::app);
  Ekin.open("output_ekin.dat",ios::app);
  Temp.open("output_temp.dat",ios::app);
  Etot.open("output_etot.dat",ios::app);
  Press.open("output_press.dat",ios::app);

  v = 0.0; //reset observables
  t = 0.0;
  for(int i=0; i<nbins; i++){
	walker[i]=0;
  }
//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( x[i] - x[j] );
     dy = Pbc( y[i] - y[j] );
     dz = Pbc( z[i] - z[j] );

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);
//per il calcolo della g(r)
     if (dr<box/2.0){
     	bin=(int)(dr/bin_size);
      	walker[bin]=walker[bin]+2;
     }

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);

//Potential energy
       v += vij;
     }
    }          
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    stima_pot = v/(double)npart; //Potential energy
    stima_kin = t/(double)npart; //Kinetic energy
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_etot = (t+v)/(double)npart; //Total enery

    v=0;
    for (int i=0; i<npart-1; ++i){
    	for (int j=i+1; j<npart; ++j){

     		dx = Pbc( x[i] - x[j] );
     		dy = Pbc( y[i] - y[j] );
     		dz = Pbc( z[i] - z[j] );

     		dr = dx*dx + dy*dy + dz*dz;
     		dr = sqrt(dr);

     if(dr < rcut){
       vij = 48.0/pow(dr,12) - 24.0/pow(dr,6);
       v += vij;
     }
    }          
  }

    stima_pressione=rho*stima_temp+(1/(3.*(double)vol))*(v);   //Pressure

    Epot << stima_pot  << endl;
    Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;
    Press << stima_pressione <<endl;

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();
    Press.close();
    return;
}
void Accumula(){
	ave_pot+=stima_pot;
	ave_kin+=stima_kin;
	ave_temp+=stima_temp;
	ave_etot+=stima_etot;
	ave_pressione+=stima_pressione;
	for(int i=0; i<nbins; ++i){
     		blk_av[i] = blk_av[i] + walker[i];
   	}
   	blk_norm = blk_norm + 1.0;
}


void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
  if(restart==0) {
	WriteConf.open("config.pre");
	for (int i=0; i<npart; ++i){
   	 WriteConf << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
        }
  WriteConf.close();
  }
  return;
}

void ConfStart(char* filename){ //Write final configuration
  ofstream WriteConf;

  cout << "Print start configuration to file config.start " << endl << endl;
  WriteConf.open(filename);

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
void Blocco(int b){
	double stddev;
	ofstream out;
	out.open("ave_epot.out", ios::app);
 	media_pot+=eKb*1.380649*pow(10,-23)*ave_pot/(double)(L/10.);	
	media_pot2+=(eKb*1.380649*pow(10,-23)*ave_pot/(double)(L/10.))*(eKb*1.380649*pow(10,-23)*ave_pot/(double)(L/10.));
	if (b==1) {out<<b<<"	"<<media_pot<<"	"<<0<<endl;}
	else{stddev=sqrt(((media_pot2/double(b))-media_pot*media_pot/((double)b*b))/(double)(b-1));
	out<<b<<"	"<<media_pot/(double)b<<"	"<<stddev<<endl;
	}
	out.close();
	
	out.open("ave_ekin.out", ios::app);
	media_kin+=eKb*1.380649*pow(10,-23)*ave_kin/(double)(L/10);	
	media_kin2+=(eKb*1.380649*pow(10,-23)*ave_kin/(double)((L/10)))*(eKb*1.380649*pow(10,-23)*ave_kin/(double)((L/10)));
	if (b==1) {out<<b<<"	"<<media_kin<<"	"<<0<<endl;}
	else {stddev=sqrt((media_kin2/((double)b)-media_kin*media_kin/((double)b*b))/(double)(b-1));
	out<<b<<"	"<<media_kin/(double)b<<"	"<<stddev<<endl;
	}
	out.close();

	out.open("ave_temp.out", ios::app);
	media_temp+=eKb*ave_temp/(double)((L/10));	
	media_temp2+=(eKb*ave_temp/(double)((L/10)))*(eKb*ave_temp/(double)((L/10)));
	if (b==1) {out<<b<<"	"<<media_temp<<"	"<<0<<endl;}
	else {stddev=sqrt((media_temp2/(double)(b)-media_temp*media_temp/(double)(b*b))/(double)(b-1));
	out<<b<<"	"<<media_temp/(double)b<<"	"<<stddev<<endl;
	}
	out.close();

	out.open("ave_etot.out", ios::app);
	media_etot+=eKb*1.380649*pow(10,-23)*ave_etot/(double)(L/10);	
	media_etot2+=(eKb*1.380649*pow(10,-23)*ave_etot/(double)(L/10))*(eKb*1.380649*pow(10,-23)*ave_etot/(double)(L/10));
	if (b==1) {out<<b<<"	"<<media_etot<<"	"<<0<<endl;}
	else {stddev=sqrt((media_etot2/((double)b)-media_etot*media_etot/((double)b*b))/(double)(b-1));
	out<<b<<"	"<<media_etot/(double)b<<"	"<<stddev<<endl;
	}
	out.close();
	
	out.open("ave_pressione.out", ios::app);
	media_pressione+=(eKb/(sigma*sigma*sigma))*1.380649*pow(10,-23)*ave_pressione/(double)((L/10));	
	media_pressione2+=((eKb/(sigma*sigma*sigma))*1.380649*pow(10,-23)*ave_pressione/(double)((L/10)))*((eKb/(sigma*sigma*sigma))*1.380649*pow(10,-23)*ave_pressione/(double)((L/10)));
	if (b==1) {out<<b<<"	"<<media_pressione<<"	"<<0<<endl;}
	else{stddev=sqrt((media_pressione2/(double)b-media_pressione*media_pressione/((double)b*b))/double(b-1));
	out<<b<<"	"<<media_pressione/(double)b<<"	"<<stddev<<endl;
	}
	out.close();

	for (int i=0; i<nbins; i++){
	stima_gdir=(double)(3.0/(rho*npart*4.0*M_PI*((bin_size*i)*(bin_size*i)*(bin_size*i)-(bin_size*(i-1))*(bin_size*(i-1))*(bin_size*(i-1)))))*blk_av[i]/blk_norm;
	glob_av[i]+=stima_gdir;
	glob_av2[i]+=stima_gdir*stima_gdir;
    	}

        if (b==nblocks) {
		out.open("ave_gdir.out");
		for(int i=0; i<nbins; i++){
		stddev=sqrt((glob_av2[i]/(double)b-glob_av[i]*glob_av[i]/((double)b*b))/double(b-1));
	    	out<<(bin_size*(i)+bin_size/2.0)*sigma<<"	"<<glob_av[i]/(double)b<<"	"<<stddev<<endl;
		}
		out.close();
	}
	
	ave_pot=0;
	ave_kin=0;
	ave_temp=0;
	ave_etot=0;
	ave_pressione=0;
	for(int i=0; i<nbins; ++i){
     		blk_av[i] = 0;
   	}
   	blk_norm = 0;
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
