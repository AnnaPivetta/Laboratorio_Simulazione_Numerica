#include "classi.h" 
#include <cstdlib>
#include <iostream>
#include <cmath>

using namespace std;

Individuo::Individuo(){
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

double Individuo::Get_lunghezza(){
	return _lunghezza;
}

double Individuo::Check (){
	int ripetizioni=0;
	for (int i=0; i<_dim; i++){
		for(int j=0; j<_dim; j++){
			if(_geni[j]==_geni[i] && i!=j) {ripetizioni++;}
		}
	}
	//cout<<"CHECK SINGOLO= "<<ripetizioni<<endl;
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
Popolazione::~Popolazione (){}
Popolazione::Popolazione() {
_dim=0;
}
Popolazione::Popolazione(int dimensione) { //setta la popolazione di dimensione dim*dim e gli individuui di dimensione dim
	_dim=dimensione;
	_popolazione.resize(_dim*_dim);
	for(int i=0; i<_dim*_dim; i++){
		_popolazione[i].Set_dimensione(_dim);
		_popolazione[i].Inizializza(_dim);
}
}
/*bool Individuo::comp(Individuo a, Individuo b) {
	return (a.Get_lunghezza()<b.Get_lunghezza());
}
	
void Popolazione::Ordina() { //ordina _popolazione dall'individuo con _lunghezza più breve a quello con _lunghezza pià lunga
	sort(_popolazione.begin(), _popolazione.end(), _popolazionecomp);
}*/

void Allele::Set_citta(int citta){
	_citta=citta;
}
void Allele::Set_indice(int indice){
	_indice=indice;
}
int Allele::Get_citta(){
	return _citta;
}
int Allele::Get_indice(){
	return _indice;
}	
	







