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

   /*for(int i=0; i<20; i++){
      cout << rnd.Rannyu() << endl;
   }*/
//------------------------------------------------------------------------------------------------------------
	int M=10000;
	int N=100;
	int L=M/N;
	double S0=100; //asset price
	double T=1; //delivery time
	double K=100; //strike price
	double r=0.1; //risk-free interest rate
	double sigma=0.25; //volatilità

	double S=0;
	double maxC1=0;
	double C1=0;
	double C1_ave=0;
	double C1_2=0;
	double maxP1=0;
	double P1=0;
	double P1_ave=0;
	double P1_2=0;
	double maxC2=0;
	double C2=0;
	double C2_ave=0;
	double C2_2=0;
	double maxP2=0;
	double P2=0;
	double P2_ave=0;
	double P2_2=0;
	
	int nsteps = 100;
	double step=T/(double)nsteps;
	double Snew=S0;
	double z=0;	
	
	ofstream out_C1, out_C2, out_P1, out_P2;
	out_C1.open ("C1.dat");
	out_C2.open ("C2.dat");
	out_P1.open ("P1.dat");
	out_P2.open ("P2.dat");
	
	for (int i=0; i<N; i++) {
		C1=0;
		P1=0;
		C2=0;
		P2=0;
		for (int k=0; k<L; k++){
			S=S0*exp( (r-0.5*sigma*sigma)*T+sigma*rnd.Gauss(0,T) );
			if(0>(S-K)) {maxC1=0;}
			else {maxC1=S-K;}
			C1+=(exp(r*T*-1))*maxC1/(double)L;
			if(0>(K-S)) {maxP1=0;}
			else {maxP1=K-S;}
			P1+=(exp(r*T*-1))*maxP1/(double)L;
			
			Snew=S0;
			z=0;	
			for (int j=0; j<nsteps-1; j++){
				z=rnd.Gauss(0,1);
				S=Snew;
				Snew=S*exp((r-0.5*sigma*sigma)*step+sigma*z*sqrt(step) );
			}	
			if(0>(Snew-K)) {maxC2=0;}
			else {maxC2=Snew-K;}
			C2+=(exp(r*T*-1))*maxC2/(double)L;
			if(0>(K-Snew)) {maxP2=0;}
			else {maxP2=K-Snew;}
			P2+=(exp(r*T*-1))*maxP2/(double)L;
		}
		C1_ave+=C1;
		C2_ave+=C2;
		P1_ave+=P1;
		P2_ave+=P2;
		C1_2+=C1*C1;
		C2_2+=C2*C2;
		P1_2+=P1*P1;
		P2_2+=P2*P2;
	
		if (i==0) {
			out_C1<<(L*(i+1))<<","<<C1_ave/(double)(i+1)<<","<<0<<endl;
			out_C2<<(L*(i+1))<<","<<C2_ave/(double)(i+1)<<","<<0<<endl;	
			out_P2<<(L*(i+1))<<","<<P2_ave/(double)(i+1)<<","<<0<<endl;
			out_P1<<(L*(i+1))<<","<<P1_ave/(double)(i+1)<<","<<0<<endl;
		}
		else {
			out_C1<<(L*(i+1))<<","<<C1_ave/(double)(i+1)<<","<<sqrt( (C1_2/(double)(i+1)-(C1_ave*C1_ave)/(double)((i+1)*(i+1)) )/((double)i) )<<endl;
			out_C2<<(L*(i+1))<<","<<C2_ave/(double)(i+1)<<","<<sqrt( (C2_2/(double)(i+1)-(C2_ave*C2_ave)/(double)((i+1)*(i+1)) )/((double)i) )<<endl;
			out_P1<<(L*(i+1))<<","<<P1_ave/(double)(i+1)<<","<<sqrt( (P1_2/(double)(i+1)-(P1_ave*P1_ave)/(double)((i+1)*(i+1)) )/((double)i) )<<endl;
			out_P2<<(L*(i+1))<<","<<P2_ave/(double)(i+1)<<","<<sqrt( (P2_2/(double)(i+1)-(P2_ave*P2_ave)/(double)((i+1)*(i+1)) )/((double)i) )<<endl;

		}				
	}

	out_C1.close();
	out_C2.close();
	out_P1.close();
	out_P2.close();

   rnd.SaveSeed();
   return 0;
}
