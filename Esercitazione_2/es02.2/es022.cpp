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
//----------------------------RANDOM WALK 3D DISCRETO------------------------------------------------------------
	
	int simulazioni=10000;
	int N=100; //numero di step per ogni simulazione
	double dir=0;
	double sgn=0;
	double x=0;
	double y=0;
	double z=0;
	int blocchi=100;
	int L=simulazioni/blocchi;
	double sum=0;

	double** r2blocchi= new double*[blocchi];    //dichiaro una matrice 100x100  
	for (int i=0; i<blocchi; i++) {
		r2blocchi[i]=new double [N];
	}

	for (int i=0; i<blocchi; i++){ 	     //inizializzo la matrice a zero
		for (int j=0; j<N; j++){
		r2blocchi[i][j]=0;
		}
	}

	for (int k=0; k<blocchi; k++){

		for (int j=0; j<L; j++){

			x=0;
			y=0;
			z=0;

			for (int i=0; i<N; i++){
				dir=rnd.Rannyu(0,3);
				sgn=rnd.Rannyu(0,2);
	
				if(dir<1){
					if(sgn<1) {x+=-1;}
					if(sgn>1) {x+=1;}
				}
				if(dir<2 && dir>1){
					if(sgn<1) {y+=-1;}
					if(sgn>1) {y+=1;}
				}
				if(dir<3 && dir>2){
					if(sgn<1) {z+=-1;}
					if(sgn>1) {z+=1;}
				}
				r2blocchi[k][i]+=(pow(x,2)+pow(y,2)+pow(z,2))/(double)L;
			}
		}
	}

	ofstream out;
	out.open ("randomwalk.dat");
	sum=0;
	double media=0;
	double sum2=0;
	double var=0;
	double stddev=0;

	for (int j=0; j<N; j++){
		sum=0;
		sum2=0;
		for (int i=0; i<blocchi; i++){	
			sum+=r2blocchi[i][j];
			sum2+=pow(r2blocchi[i][j],2);
		}
		media=sum/(double)blocchi;
		var=sum2/(double)blocchi-pow(media,2);
		stddev=sqrt(var/(int)(blocchi-1));
		out<<j+1<<","<<sqrt(media)<<","<<(1/sqrt(media))*stddev<<endl;
	}

//-------------------------------------RANDOM WALK 3D CONTINUO----------------------------------------------
	double theta=0;
	double fi=0;
	double a=0;

	for (int i=0; i<blocchi; i++){ 	     //inizializzo la matrice a zero
		for (int j=0; j<N; j++){
		r2blocchi[i][j]=0;
		}
	}

	for (int k=0; k<blocchi; k++){

		for (int j=0; j<L; j++){

			x=0;
			y=0;
			z=0;

			for (int i=0; i<N; i++){
			a=rnd.Rannyu(0,1);
			theta=acos(1-2*a); //genero theta tra 0 e pi
			fi=rnd.Rannyu(0,2*M_PI); //genero fi tra 0 e 2pi
			x+=sin(theta)*cos(fi);
			y+=sin(theta)*sin(fi);
			z+=cos(theta);
			r2blocchi[k][i]+=(pow(x,2)+pow(y,2)+pow(z,2))/(double)L;
			}
		}
	}
	out.close();

	//--------------------------------------istogramma per verificare andamento diffusivo-------------
	out.open("posfinal.dat");
	/*for(int i=0; i<blocchi; i++){
		out<<sqrt(r2blocchi[i][99])<<endl;
	}*/
	for (int i=0; i<10000; i++){   //aggiungo altri 1000 punti all'istogramma, per avere più statistica
			x=0;
			y=0;
			z=0;
			for (int j=0; j<10000; j++){
			a=rnd.Rannyu(0,1);
			theta=acos(1-2*a); //genero theta tra 0 e pi
			fi=rnd.Rannyu(0,2*M_PI); //genero fi tra 0 e 2pi
			x+=sin(theta)*cos(fi);
			y+=sin(theta)*sin(fi);
			z+=cos(theta);
			}
			out<<x<<endl;
	}
	out.close();
	//------------------------------------------------------------------------------------------------

	out.open ("randomwalk.dat");
	sum=0;
	media=0;
	sum2=0;
	var=0;
	stddev=0;

	for (int j=0; j<N; j++){
		sum=0;
		sum2=0;
		for (int i=0; i<blocchi; i++){	
			sum+=r2blocchi[i][j];
			sum2+=pow(r2blocchi[i][j],2);
		}
		media=sum/(double)blocchi;
		var=sum2/(double)blocchi-pow(media,2);
		stddev=sqrt(var/(int)(blocchi-1));
		out<<j+1<<","<<sqrt(media)<<","<<(1/sqrt(media))*stddev<<endl;
	}
	
   rnd.SaveSeed();
   return 0;
}
