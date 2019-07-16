#include "classi.h" 
#include <cstdlib>
#include <iostream>
#include <cmath>

using namespace std;

int Individuo::Pbc(int indice){
	int appo=indice-_dim;
	if (indice<_dim) {return indice;}	
	else {
		do {appo=appo-_dim;}
		while(appo>=_dim);
		return appo;
	}
}
void Individuo::Pair_Mutation(){
	int indice_1=(int)(rnd.Rannyu(0,_dim));
	int indice_2=(int)(rnd.Rannyu(0,_dim));
	int appo=_geni[indice_1];
	_geni[indice_1]=_geni[indice_2];
	Set_Elemento(appo, indice_2);
}
void Individuo::Simple_Shift(){
	int n=(int)(rnd.Rannyu(0,_dim)); //l'elemento n-esimo diventa il primo del percorso
	rotate(_geni.begin(),_geni.begin()+n,_geni.end());
}
void Individuo::Inversion(){
	int inizio=(int)(rnd.Rannyu(0,_dim));
	int quanti=(int)(rnd.Rannyu(0,_dim));
	reverse(_geni.begin()+inizio, _geni.begin()+Pbc(inizio+quanti));
}
void Individuo::Permutation(){
	int estremo_1=(int)(rnd.Rannyu(0,_dim/2));
	int lunghezza_segmento=(int)(rnd.Rannyu(0,Pbc(_dim/2 -estremo_1)));	        
	int estremo_2=(int)(rnd.Rannyu(_dim/2, Pbc(_dim-lunghezza_segmento)));
	swap_ranges(_geni.begin()+(estremo_1),_geni.begin()+estremo_1+lunghezza_segmento,_geni.begin()+estremo_2);	
}
int Individuo::Get_Dimensione(){
	return _dim;
}
Individuo::Individuo(){
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
  _dim=0;
  _lunghezza=0;
}
Individuo::~Individuo (){}
Individuo::Individuo(int dimensione){
	_dim=dimensione;
}
void Individuo::Inizializza(int dimensione){
	_dim=dimensione;
	_geni.resize(_dim);
	for(int i=0; i<_dim; i++){
		_geni[i]=i;
	}
	random_shuffle (_geni.begin(), _geni.end() );
}

void Individuo::Calcola_lunghezza(double* x, double* y){
	_lunghezza=0;
	double xA, xB, yA, yB;
	for (int j=0; j<_dim-1; j++){
		xA=x[_geni[j+1]];
		xB=x[_geni[j]];
		yA=y[_geni[j+1]];
		yB=y[_geni[j]];
		_lunghezza+=((xA-xB)*(xA-xB)+(yA-yB)*(yA-yB));
	}
	_lunghezza+=((x[_geni[_dim-1]]-x[_geni[0]])*(x[_geni[_dim-1]]-x[_geni[0]])+(y[_geni[_dim-1]]-y[_geni[0]])*(y[_geni[_dim-1]]-y[_geni[0]]));
	
}
void Individuo::Set_Elemento(double elemento, int posizione){
	_geni[posizione]=elemento;
}
int Individuo::Get_Elemento(int posizione){
	return _geni[posizione];
}

double Individuo::Get_lunghezza(double* x, double*y){
	_lunghezza=0;
	double xA, xB, yA, yB;
	for (int j=0; j<_dim-1; j++){
		xA=x[_geni[j+1]];
		xB=x[_geni[j]];
		yA=y[_geni[j+1]];
		yB=y[_geni[j]];
		_lunghezza+=((xA-xB)*(xA-xB)+(yA-yB)*(yA-yB));
	}
	_lunghezza+=((x[_geni[_dim-1]]-x[_geni[0]])*(x[_geni[_dim-1]]-x[_geni[0]])+(y[_geni[_dim-1]]-y[_geni[0]])*(y[_geni[_dim-1]]-y[_geni[0]]));
	return _lunghezza;
}

double Individuo::Check (){
	int ripetizioni=0;
	for (int i=0; i<_dim; i++){
		for(int j=0; j<_dim; j++){
			if(_geni[j]==_geni[i] && i!=j) {ripetizioni++;}
		}
	}
	cout<<"CHECK SINGOLO= "<<ripetizioni<<endl;
	return ripetizioni;
}
void Individuo::Set_dimensione(int dimensione){
	_dim=dimensione;
	_geni.resize(_dim);
}
vector<int>::iterator Individuo::Get_start(){
	return _geni.begin();
}
vector<int>::iterator Individuo::Get_end(){
	return _geni.end();
}






