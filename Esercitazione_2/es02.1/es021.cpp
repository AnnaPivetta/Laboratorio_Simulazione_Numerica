#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

using namespace std;
 
int main (int argc, char *argv[]){

   Random rnd;
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

  /* for(int i=0; i<20; i++){
      cout << rnd.Rannyu() << endl;
   }*/
//----------------------------PARTE UNO, uniform sampling------------------------------------------------------------
	
	int M = 100000;
	int N= 100;
	int L= M/N; //numero di punti per ogni blocco
	double integral_ave=0;
	double integral2_ave=0;
	double integral;
	double sum=0;
	double x;
	double f;
	
	ofstream out;
	out.open ("integraleunif.dat");
	

	for (int k=0; k<N; k++) {
		sum=0;
		for (int i=0; i<L; i++) {
			x =rnd.Rannyu(0,1);
			f = M_PI*0.5*cos(x*M_PI*0.5);
			sum+=f;
		}
	
		integral=(sum/(double)L);  //media del blocco
		integral_ave+=integral;
		integral2_ave+=integral*integral;
		if (k==0) {out<<(L*(k+1))<<","<<integral_ave/(double)(k+1)<<","<<0<<endl;}
		else {out<<(L*(k+1))<<","<<integral_ave/(double)(k+1)<<","<<sqrt( (integral2_ave/(double)(k+1)-(integral_ave*integral_ave)/(double)((k+1)*(k+1)) )/((double)k) )<<endl;}
	}	
	out.close();

//---------------------------PARTE DUE, importance sampling------------------------------------------

	
	out.open ("integraleimp.dat");
	integral_ave=0;
	integral2_ave=0;
	double x0;
	
	for (int k=0; k<N; k++) {
		sum=0;
		for (int i=0; i<L; i++) {
			x0 =rnd.Rannyu(0,1);
			x=1-sqrt(1-x0);
			f = (M_PI*0.5*cos(x*M_PI*0.5))/(-2*x+2);
			sum+=f;
		}
	
		integral=(sum/(double)L);  //media del blocco
		integral_ave+=integral;
		integral2_ave+=integral*integral;
		if (k==0) {out<<(L*(k+1))<<","<<integral_ave/(double)(k+1)<<","<<0<<endl;}
		else {out<<(L*(k+1))<<","<<integral_ave/(double)(k+1)<<","<<sqrt( (integral2_ave/(double)(k+1)-(integral_ave*integral_ave)/(double)((k+1)*(k+1)) )/((double)k) )<<endl;}
	}	
	out.close();


   rnd.SaveSeed();
   return 0;
}
