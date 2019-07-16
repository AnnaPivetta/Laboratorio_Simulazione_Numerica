#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include "random.h"

using namespace std;
Random rnd;
double mu;
double sigma;
int nbins;
//----per fare l'istogramma----------
double*occorrenze;
double*estremi;
double sup, inf, larghezza;
//----------------------------------
int conta;
int L;
double dens_prob; //osservabile densità di probabilità
double O; //osservabile energia
double psi_quadro (double, double, double);
double psi (double, double, double);
double psi_der2 (double, double, double);
void Input ();
void Passo ();
void Accumula();
void Blocks(int blocco);
double x,delta, throws, blocks;
double xtest;
double energia, energia2;
double energia_blocco;


void Equilibrazione ();
void Istogramma(double);

