#include "es081.h"

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

//----------------------------  Es 8: VARIATIONAL MONTE CARLO	  ------------------------------------------------------------

	system ("rm blocchi.dat isto.dat");
	Input();
	conta=0;
	Equilibrazione();	
	energia=0;
	energia2=0;
	L=throws/blocks;
	conta=0;
	for (int i=0; i<blocks; i++) {
		energia_blocco=0;
		for (int j=0; j<L; j++) {
			Passo();
			Istogramma(x);
			Accumula();
		}
		Blocks(i+1);
	}
	cout<<"frazione di mosse accettate= "<<(double)conta/(double)throws<<endl;
	cout<<"delta="<<delta<<endl;
	ofstream out_isto;
	out_isto.open("isto.dat");
	for(int p=0; p<nbins; p++) {  //stampo su file di output per istogramma
		out_isto<<estremi[p]<<"	"<<occorrenze[p]/((double)(larghezza*throws))<<endl;
	}
   out_isto.close();
   delete [] occorrenze;
   delete [] estremi;

   rnd.SaveSeed();
   return 0;
}

double psi_quadro (double r,double mu,double sigma) { //ritorna il modulo quadro della psi trial
	return (exp(-(r-mu)*(r-mu)/(2*sigma*sigma))+exp(-(r+mu)*(r+mu)/(2*sigma*sigma)))*(exp(-(r-mu)*(r-mu)/(2*sigma*sigma))+exp(-(r+mu)*(r+mu)/(2*sigma*sigma)));
}

double psi(double r,double mu,double sigma) { //ritorna psi trial
	return (exp(-(r-mu)*(r-mu)/(2*sigma*sigma))+exp(-(r+mu)*(r+mu)/(2*sigma*sigma)));
}

double psi_der2 (double r,double mu,double sigma){ //ritorna la derivata seconda di psi trial
	return (exp(-(r-mu)*(r-mu)/(2*sigma*sigma)))*(((mu-r)/(sigma*sigma))*((mu-r)/(sigma*sigma))-1/(sigma*sigma)) +(exp(-(r+mu)*(r+mu)/(2*sigma*sigma)))*(((-mu-r)/(sigma*sigma))*((-mu-r)/(sigma*sigma))-1/(sigma*sigma));
}
void Input (){
	ifstream in;
	in.open ("input.dat");
	in>>mu;
	in>>sigma;
	in>>delta;
	in>>x;   //punto di partenza 
	in>>throws;
	in>>blocks;
	in>>nbins;
	occorrenze=new double [nbins];
	sup=3.0;
	inf=-3.0;
	estremi= new double [nbins];
	larghezza=(sup-inf)/(double)nbins;
	for(int k=0; k<nbins; k++) { //inizializzo a zero il vettore occorrenze, per fare l'istogramma
		estremi[k]=inf+k*larghezza;
		occorrenze[k]=0;
	}
	
}

void Passo () { //Metropolis
	xtest=rnd.Rannyu(x-delta, x+delta);
	double alpha;
	alpha=min(1.,psi_quadro(xtest, mu, sigma)/psi_quadro(x, mu, sigma)); 
	double p=rnd.Rannyu(0,1);
	if (p<=alpha) {
		x=xtest;
		conta++;	
	}
	O=-0.5*psi_der2(x, mu, sigma)+((x*x*x*x)-(2.5)*x*x)*psi(x,mu,sigma);
	O=O/psi(x, mu, sigma);
}
void Accumula (){
	energia_blocco+=O;
}
void Blocks (int blocco){

	energia+=energia_blocco/(double)L;
	energia2+=(energia_blocco/(double)L)*(energia_blocco/(double)L);
	ofstream out;
	out.open("blocchi.dat", ios::app);
	double media=energia/((double)blocco);
	double stddev=0;
	if (blocco==1) {stddev=0;}
	else {stddev=sqrt((energia2/((double)blocco)-media*media)/((double)blocco-1));}
	out<<blocco<<"	"<<media<<"	"<<stddev<<endl;
}	
void Equilibrazione(){
	double lanci=5000;
	for (int j=0; j<lanci; j++) {
			Passo();
		}
	double rapporto=(double)conta/(double)lanci;
	if (rapporto>=0.47 && rapporto<=0.53) {return;} 
	do{	
		conta=0;
		delta+=0.02;
		for (int j=0; j<lanci; j++) {
			Passo();
		}
		rapporto=(double)conta/(double)lanci;
	}
	while (rapporto<=0.47 || rapporto>=0.53); 
}
void Istogramma(double x){	
	for (int t=0; t<nbins-1; t++) {  //assegno x a un bin
		if (estremi[t]<=(x) && (estremi[t+1]>x)) {occorrenze[t]++;}
	}
	if ( x>=estremi[nbins-1] && x<sup) {occorrenze[nbins-1]++;}	
}

