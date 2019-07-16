#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <time.h>  //per usare il generatore di numeri casuali del c++
#include "random.h"
#include "classi.h"
//#include "classi.cpp"

using namespace std;
Random rnd;
double* selection_array;
double number;
int selected_indice;
Individuo Elitario;

int forma;
double p_inversion;
double p_pair_mutation;
double p_simple_shift;
double p_shift;
double p_permutation;
double p_crossover;
double p_elitario;
int n_generazioni;

int indice_taglio;
int indice_padre;
int indice_madre;
double normalizzazione;
int n_citta; //numero di città, per ora lo scrivo io poi si potrebbe leggere da file di input
int controllo; //0 se i constrain sono rispettati, diverso da zero altrimenti 

//void Inizializzazione();
//double lunghezza(); //calcola la lunghezza di un percorso
//int check (vector<vector<int> >);


