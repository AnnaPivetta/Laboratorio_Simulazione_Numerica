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
//--------------------------------------------------------------------------------------------------------


	double u=0.0;
	double accumula=0.0;
	int nbins=50;
	double estremi[nbins];
	double occorrenze[nbins];
	double larghezza;
	double inf;
	double sup;
	ofstream out;

	//---------------------------------------uniforme-------------------------------------
	for (int i=0; i<102; i++){
		if(i==1 || i==2 || i==10 || i==100){
			if (i==1) {inf=-0.1; sup=1.1; out.open("unifN1.dat");}
			if (i==2) {inf=-0.1; sup=1.1; out.open("unifN2.dat");}
			if (i==10){inf=0.1; sup=0.9; out.open("unifN10.dat");}
			if(i==100){inf=0.35; sup=0.65; out.open("unifN100.dat");}
			larghezza=(sup-inf)/(double)nbins;
			for(int k=0; k<nbins; k++) {
				estremi[k]=inf+k*larghezza;
				occorrenze[k]=0;
			}

			for (int j=0; j<10000; j++){    
				accumula=0;	
				for(int p=0; p<i; p++){     //generero Sn	
					u=rnd.Rannyu(0,1);
					accumula+=u;				
				}
				for (int t=0; t<nbins-1; t++) {  //assegno Sn a un bin

					if (estremi[t]<=(accumula/(double)i) && (estremi[t+1]>accumula/(double)i)) {occorrenze[t]++;}
					
				}
				if ( (accumula/(double)i)>=estremi[nbins-1] && (accumula/(double)i)<sup) {occorrenze[nbins-1]++;}
			}
			for(int p=0; p<nbins; p++) {  //stampo su file di output
				out<<estremi[p]<<"	"<<occorrenze[p]<<endl;
			}		
			out.close();				
		}
	}
	
	//---------------------------------------------esponenziale------------------------------------------------
	for (int i=0; i<102; i++){
		if(i==1 || i==2 || i==10 || i==100){
			if (i==1) {inf=-0.5; sup=5.7; out.open("expN1.dat");}
			if (i==2) {inf=-0.5; sup=5.7; out.open("expN2.dat");}
			if (i==10){inf=-0.5; sup=3.0; out.open("expN10.dat");}
			if(i==100){inf= 0.4; sup=1.6; out.open("expN100.dat");}
			larghezza=(sup-inf)/(double)nbins;
			for(int k=0; k<nbins; k++) {
				estremi[k]=inf+k*larghezza; //esremi sinistri dei bin
				occorrenze[k]=0;
			}

			for (int j=0; j<10000; j++){    
				accumula=0;	
				for(int p=0; p<i; p++){     //generero Sn	
					u=rnd.Exp(1);
					accumula+=u;				
				}
				for (int t=0; t<nbins-1; t++) {  //assegno Sn a un bin

					if (estremi[t]<=(accumula/(double)i) && (estremi[t+1]>accumula/(double)i)) {occorrenze[t]++;}
					
				}
				if (accumula/(double)i>=estremi[nbins-1] && (accumula/(double)i)<sup) {occorrenze[nbins-1]++;}
			}
			for(int p=0; p<nbins; p++) {  //stampo su file di output
				out<<estremi[p]<<"	"<<occorrenze[p]<<endl;
			}		
			out.close();				
		}
	}

	//---------------------------------------------------Cauchy-Lorentz--------------------------------------------
	
	nbins=100;
	double estremi_1[nbins];
	double occorrenze_1[nbins];
	for (int i=0; i<102; i++){
		if(i!=1 && i!=2 && i!=10 && i!=100) continue;
		
			if (i==1) {inf=-22.0; sup=22.0; out.open("LorentzN1.dat");}
			if (i==2) {inf=-22.0; sup=22.0; out.open("LorentzN2.dat");}
			if (i==10){inf=-22.0; sup=22.0; out.open("LorentzN10.dat");}
			if(i==100){inf=-22.0; sup=22.0; out.open("LorentzN100.dat");}
			larghezza=(sup-inf)/(double)nbins;
			for(int k=0; k<nbins; k++) {
				estremi_1[k]=inf+k*larghezza;
				occorrenze_1[k]=0;
			}
	
			for (int j=0; j<10000; j++){    
				accumula=0;	
				for(int p=0; p<i; p++){     //generero Sn	
					u=rnd.Lorentz(1,0);
					accumula+=u;				
				}
				for (int t=0; t<nbins-1; t++) {  //assegno Sn a un bin
					if (estremi_1[t]<=(accumula/(double)i) && (estremi_1[t+1]>accumula/(double)i)) {occorrenze_1[t]++;}
				}
				if ( (accumula/(double)i)>=estremi_1[nbins-1] && (accumula/(double)i)<sup) {occorrenze_1[nbins-1]++;}
				
			}
			for(int p=0; p<nbins; p++) {  //stampo su file di output
				out<<estremi_1[p]<<"	"<<occorrenze_1[p]<<endl;
			}		
			out.close();				
		
	}

   rnd.SaveSeed();
   return 0;
}

