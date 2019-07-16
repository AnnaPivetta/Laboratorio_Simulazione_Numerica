#ifndef _CLASSI_H_
#define _CLASSI_H_
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <iterator>
#include "random.h"
using namespace std;

class Individuo {

	public:
	Individuo();
	Individuo(int);
	~Individuo ();
	void Set_dimensione(int);
	void Calcola_lunghezza(double*, double*);
	double Get_lunghezza(double*, double*);
	void Inizializza(int);
	double Check(); //ritorna zero se l'individuo rispetta i constrain e un numero diverso da zero altrimenti
	void Set_Elemento(double, int);
	int Get_Elemento(int);
	int Get_Dimensione();
	vector<int>::iterator Get_start(); //ho bisogno che il tipo di ritorno sia un iteratore per usare la finzione "rotate" 
	vector<int>::iterator Get_end();
	//bool comp (Individuo, Individuo);
	void Pair_Mutation();
	void Simple_Shift();
	void Inversion();
	void Permutation();
	int Pbc(int);
	private:
	Random rnd; //serve per le mutazioni
	vector<int> _geni;
	double _lunghezza;
	int _dim;
};

/*
class Allele {
	
	public:
	void Set_citta(int);
	void Set_indice(int);
	int Get_citta();
	int Get_indice();

	private:
	int _citta;
	int _indice;
};*/
#endif
