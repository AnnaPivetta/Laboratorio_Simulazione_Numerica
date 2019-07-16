#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include "random.h"

using namespace std;
Random rnd;
int conta;
int eq;
int distribuzione;
double psi1S (double);
void Input ();
void Passo ();
void Accumula();
void Blocks(int blocco);
string nome_stato;
double x,y,z,delta, throws, blocks;
double xtest, ytest, ztest;
double r, stddev;
double X, Y, Z;
double media, media2;
double psi2P (double, double, double);
void Equilibrazione ();
//void Grafico3D (int);
