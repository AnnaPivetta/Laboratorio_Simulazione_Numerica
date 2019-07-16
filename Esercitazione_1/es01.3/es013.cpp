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
//----------------------------ESPERIMENTO DI BUFFON------------------------------------------------------------
	
	//genero un numero distribuito uniformemente tra 0 e 2pi
	double x=0;
	double y=0;
	double theta=0;
	double pos=0;
	double l=1;
	double d=1.6;
	double pos1=0;
	int hit=0;
	int M=10000;
	int N=100;
	int L=M/N;
	double pi_b=0.0;
	double pi_ave=0.0;
	double pi2_ave=0.0;

	ofstream out;
	out.open ("pi.dat");

	for (int j=0; j<N; j++){

		pi_b=0;
		hit=0;	
		for (int i=0; i<L; i++){
			do{ 
				x=rnd.Rannyu(-1, 1);
				y=rnd.Rannyu(0,1);

			} while (pow(x,2)+pow(y,2)>1);
	
			theta=2*acos(x/sqrt(pow(x,2)+pow(y,2)));  //theta angolo tra zero e 2*pi
			pos=rnd.Rannyu(0,100000);
			pos1=pos+l*sin(theta);

			if ( pos1<( (int)(pos/d)*d) || pos1>(((int)(pos/d)+1)*d)) {hit++;}
		}

	pi_b=2*l*(double)L/((double)hit*d);   //qui sto usando hit del blocco
	pi_ave+=pi_b;
	pi2_ave+=pi_b*pi_b;
	if (j==0) {out<<(L*(j+1))<<","<<pi_ave/(double)(j+1)<<","<<0<<endl;}
		else {out<<(L*(j+1))<<","<<pi_ave/(double)(j+1)<<","<<sqrt( ((pi2_ave/(double)(j+1)-(pi_ave*pi_ave)/((double)(j+1)*(j+1))) )/(double)j )<<endl;}
	}
	
   rnd.SaveSeed();
   return 0;
}
