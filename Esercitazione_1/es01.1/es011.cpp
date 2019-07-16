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
//----------------------------PARTE ZERO------------------------------------------------------------
	ofstream out;
	out.open("grafico.dat");
	
	int M=100000;	//numero totale di lanci
	int N=100;	//numero di blocchi
	int L=int(M/N);	//numero di lanci per blocco
	double randr=0.0;
	double r_ave=0.0;
	double r2_ave=0.0;	
	double r_b=0;
	for (int i=0; i<N; i++) {

		r_b=0;
		for (int k=0; k<L; k++) {
			randr=rnd.Rannyu(0,1);
			r_b+=randr/(double)L;
		}
		 
		r_ave+=r_b;
		r2_ave+=r_b*r_b;

		if (i==0) {out<<(L*(i+1))<<","<<r_ave/(double)(i+1)<<","<<0<<endl;}
		else {out<<(L*(i+1))<<","<<r_ave/(double)(i+1)<<","<<sqrt( (r2_ave/(double)(i+1)-(r_ave*r_ave)/(double)((i+1)*(i+1)) )/((double)i) )<<endl;}
	}	
	
	out.close();
//-------------------------PARTE UNO---------------------------------------------------------------

	//tengo gli stessi M, N, L di prima

	out.open ("grafico1.dat");
	
	double randsigma=0.0;
	double sigma_ave=0.0;
	double sigma2_ave=0.0;	
	double sigma_b=0;
	for (int i=0; i<N; i++) {

		sigma_b=0;
		for (int k=0; k<L; k++) {
			randsigma=pow(rnd.Rannyu(0,1)-0.5,2);
			sigma_b+=randsigma/(double)L;
		}
		 
		sigma_ave+=sigma_b;
		sigma2_ave+=sigma_b*sigma_b;

		if (i==0) {out<<(L*(i+1))<<","<<sigma_ave/(double)(i+1)<<","<<0<<endl;}
		else {out<<(L*(i+1))<<","<<sigma_ave/(double)(i+1)<<","<<sqrt( (sigma2_ave/(double)(i+1)-(sigma_ave*sigma_ave)/(double)((i+1)*(i+1)) )/((double)i) )<<endl;}
	}	
	
	out.close();

//-----------------------------PARTE DUE----------------------------------------------------------------
	out.open("grafico2.dat");
	
	M=100;
	double number=0.0;
	double O=0;
	double E=10000/(double)M;
	double chi=0;
	
	for (int k=0; k<100; k++) { 
	chi=0;	
		for (int i=1; i<101; i++){ //scorre sui bin
			O=0;
			for (int j=0; j<10000; j++){
				number=rnd.Rannyu(0,1);
				if((double)(i-1)/(double)M<number) {
					if (number<(double)i/(double)M) {O++;}
				}			
			}
		chi+=pow((O-E),2)/(double)E;	
		}
	out<<k+1<<","<<chi<<endl;
	}

	out.close();

   rnd.SaveSeed();
   return 0;
}
