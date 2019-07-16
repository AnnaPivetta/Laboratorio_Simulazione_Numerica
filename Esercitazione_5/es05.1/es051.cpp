#include "es051.h"

int main (int argc, char *argv[]){

 
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

//----------------------------  STATI 1S e 2P IDROGENO  ------------------------------------------------------------
	
	Input();
	if (nome_stato=="1S") system("rm 1Sposizioni.dat grafico3D1S.dat");						
	if (nome_stato=="2P") system("rm 2Pposizioni.dat grafico3D2P.dat");
	conta=0;
	eq=0;
	Equilibrazione();
	cout<<"ho modificato delta "<<eq<<" volte"<<endl;
	X=0;
	Y=0;
	Z=0;	
	r=0;
	media=0;
	media2=0;
	int L=throws/blocks;
	for (int i=0; i<blocks; i++) {
		for (int j=0; j<L; j++) {
		Passo();
		Accumula();
		}
	Blocks(i+1);
	}
	cout<<"frazione di mosse accettate= "<<(double)conta/(double)throws<<endl;

   rnd.SaveSeed();
   return 0;
}

double psi1S (double r) {
	return (1.0/M_PI)*exp(-2.0*r);
}
double psi2P (double x, double y, double z) {
	double radius=sqrt(x*x+y*y+z*z);
	return ((1.0/(32.0*M_PI))*radius*radius*exp(-radius)*(z*z/(radius*radius)));
}
void Input (){
	ifstream in;
	in.open ("input.dat");
	in>>delta;
	in>>x;
	in>>y;
	in>>z;
	in>>throws;
	in>>blocks;
	in>>nome_stato;
	in>>distribuzione;
}

void Passo () {
	if (distribuzione==1) {
		xtest=rnd.Rannyu(x-delta, x+delta);
		ytest=rnd.Rannyu(y-delta, y+delta);
		ztest=rnd.Rannyu(z-delta, z+delta);
	}
	if (distribuzione==2) {
		xtest=rnd.Gauss(x, delta);
		ytest=rnd.Gauss(y, delta);
		ztest=rnd.Gauss(z, delta);
	}
	double alpha;
	if (nome_stato=="1S") {
		double q=sqrt(x*x+y*y+z*z);
		double s=sqrt(xtest*xtest+ytest*ytest+ztest*ztest);
		alpha=min(1.,psi1S(s)/psi1S(q));
	}
	if (nome_stato=="2P") {alpha=min(1.,psi2P(xtest, ytest, ztest)/psi2P(x, y, z));}	
	double p=rnd.Rannyu(0,1);
	if (p<=alpha) {
		x=xtest;
		y=ytest;
		z=ztest;
		conta++;	
	}
	//grafici 3D
	if (nome_stato=="1S"){
	ofstream out1;
	out1.open("grafico3D1S.dat", ios::app);
	out1<<x<<"	"<<y<<"	"<<z<<endl;
	out1.close();
	}
	if (nome_stato=="2P"){
	ofstream out1;
	out1.open("grafico3D2P.dat", ios::app);
	out1<<x<<"	"<<y<<"	"<<z<<endl;
	out1.close();
	}
	
}
void Accumula (){
	r+=sqrt(x*x+y*y+z*z);
}
void Blocks (int blocco){
	int L=throws/blocks;
	ofstream out;
	if (nome_stato=="1S") {out.open ("1Sposizioni.dat", ios::app);}
	if (nome_stato=="2P") {out.open ("2Pposizioni.dat", ios::app);}
	media+=r/((double)(L));
	media2+=r/((double)(L))*(r/((double)(L)));
	if (blocco==1) {stddev=0;}
	else {stddev=sqrt((media2/((double)(blocco))-(media*media/((double)(blocco*blocco))))/((double)(blocco-1)));}
	out<<blocco<<"	"<<media/((double)(blocco))<<"	"<<stddev<<endl;
	out.close();
	r=0;
}	
void Equilibrazione(){
	double lanci=1000;
	for (int j=0; j<lanci; j++) {
			Passo();
		}
	double rapporto=(double)conta/(double)lanci;
	if (rapporto>=0.48 && rapporto<=0.52) {return;} 
	do{	
		conta=0;
		delta+=0.1;
		for (int j=0; j<lanci; j++) {
			Passo();
		}
		rapporto=(double)conta/(double)lanci;
		eq++;
	}
	while (rapporto<=0.48 || rapporto>=0.52); 
}	
