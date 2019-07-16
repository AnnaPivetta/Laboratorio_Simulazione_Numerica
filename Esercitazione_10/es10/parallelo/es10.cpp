#include "es10.h"
#include "mpi.h"


double Boltzman(Individuo, Individuo);
void Input();
void passo_1(Individuo&);
void passo_2(Individuo&);
void passo_3(Individuo&);
void passo_4(Individuo&);
using namespace std;
 
int main (int argc, char *argv[]){

   MPI::Init(argc,argv);
   int size = MPI::COMM_WORLD.Get_size();
   int rank = MPI::COMM_WORLD.Get_rank();

   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
	for(int i=0; i<rank+34; i++){
      		Primes >> p1 >> p2 ;
	}
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
	
	Input();
	//system("rm lunghezza.dat");
	x_citta=new double [n_citta];
	y_citta=new double [n_citta];

	if (rank==0){    //calcolo le coordinate solo per il primo nodo e le comunico agli altri, così ottimizzo lo stesso percorso e non size percorsi diversi
		if (forma==0){
			for (int i=0; i<n_citta; i++){
				double angolo=rnd.Rannyu(0,2*M_PI);
				x_citta[i]=cos(angolo);
				y_citta[i]=sin(angolo);
			}
		}
		if (forma==1){
			for (int i=0; i<n_citta; i++){
				x_citta[i]=rnd.Rannyu(0,1);
				y_citta[i]=rnd.Rannyu(0,1);
			}
		}
	}

	MPI_Bcast ( x_citta, n_citta, MPI_DOUBLE, 0, MPI::COMM_WORLD);
	MPI_Bcast ( y_citta, n_citta, MPI_DOUBLE, 0, MPI::COMM_WORLD);

//	cout<<"rank= "<<rank<<", x_1= "<<x_citta[1]<<"y_1= "<<y_citta[1]<<", x_7= "<<x_citta[7]<<"y_7= "<<y_citta[7]<<endl;
	
	Individuo percorso;
	double lunghezza;	
	percorso.Inizializza(n_citta);
	vett=new double [n_blocchi];
	for (int i=0; i<n_blocchi; i++){
		T=T*0.9;
		step_blocco=step_blocco_iniziale*(i+1);
		accettazione=0;
			for(int j=0; j<step_blocco; j++){
			passo_1(percorso);
			passo_2(percorso);
			passo_3(percorso);
			passo_4(percorso);
		}
		//cout<<"t= "<<T<<endl;
		
		vett[i]=percorso.Get_lunghezza(x_citta,y_citta);		
	}
	lunghezza=percorso.Get_lunghezza(x_citta,y_citta);
	
	Dati prop;
	prop.l=lunghezza;
	prop.r=rank;
 	Dati best;
	best.l=0.;
	best.r=0;
	
	MPI_Reduce(&prop, &best, 1, MPI_DOUBLE_INT, MPI_MINLOC, 0, MPI::COMM_WORLD );
	MPI_Bcast ( &best, 2, MPI_DOUBLE_INT, 0, MPI::COMM_WORLD);

	//cout<<"rank= "<<rank<<endl;
	if (rank==0) cout<<rank<<", lunghezza= "<<best.l<<",best rank= "<<best.r<<endl;
	
	if (rank==best.r) {
		
		ofstream out;
		out.open("lunghezza.dat");
		for (int i=0; i<n_blocchi; i++){
			out<<i+1<<"	"<<vett[i]<<endl;
		}
		out.close();
		out.open("percorso_migliore.dat");
		for (int i=0; i<n_citta; i++){
			out<<x_citta[percorso.Get_Elemento(i)]<<"	"<<y_citta[percorso.Get_Elemento(i)]<<endl;
		}
		out.close();
	}

	
	delete [] x_citta;
	delete [] y_citta;
	delete [] vett;
//--------------------------------------------------------------------------------------------------------	
   if (rank==0){rnd.SaveSeed();}
   MPI::Finalize();
   return 0;
}

void Input(){
	ifstream in;	
	in.open ("input.dat");
	in>>n_citta;
	in>>forma;
	in>>T;
	in>>n_blocchi;
	in>>step_blocco_iniziale;
}



double Boltzman(Individuo p, Individuo c){
	return exp(-(1.0/(double)(T))*(c.Get_lunghezza(x_citta, y_citta)-p.Get_lunghezza(x_citta, y_citta)));
}

void passo_1(Individuo& p){
	Individuo copia;
	copia.Inizializza(p.Get_Dimensione());
	copia=p;
	copia.Pair_Mutation();
	double alpha;
	alpha=min(1.,(double)(Boltzman(p, copia)));
	double a=rnd.Rannyu(0,1);
	if (a<=alpha) {
	p=copia;
	accettazione++;
	}
}
void passo_2(Individuo& p){
	Individuo copia;
	copia.Inizializza(p.Get_Dimensione());
	copia=p;
	copia.Simple_Shift();
	double alpha;
	alpha=min(1.,(double)(Boltzman(p, copia)));
	double a=rnd.Rannyu(0,1);
	if (a<=alpha) {
	p=copia;
	accettazione++;
	}
}
void passo_3(Individuo& p){
	Individuo copia;
	copia.Inizializza(p.Get_Dimensione());
	copia=p;
	copia.Inversion();
	double alpha;
	alpha=min(1.,(double)(Boltzman(p, copia)));
	double a=rnd.Rannyu(0,1);
	if (a<=alpha) {
	p=copia;
	accettazione++;
	}
}
void passo_4(Individuo& p){
	Individuo copia;
	copia.Inizializza(p.Get_Dimensione());
	copia=p;
	copia.Permutation();
	double alpha;
	alpha=min(1.,(double)(Boltzman(p, copia)));
	double a=rnd.Rannyu(0,1);
	if (a<=alpha) {
	p=copia;
	accettazione++;
	}
}
