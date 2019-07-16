#include "es091.h"

bool comp (Individuo, Individuo);
void Ordina (vector<Individuo>&);
void Inizializza_p(vector<Individuo>&);
void Calcola_lunghezza(vector<Individuo>&, double*, double*);
double Check (vector<Individuo>);
int Selection(vector<Individuo>);
void Crossover(vector<Individuo>&,vector<Individuo>);
void Pair_Mutation(vector<Individuo>&);
void Simple_shift(vector<Individuo>&);
void Shift(vector<Individuo>&);
void Inversion (vector<Individuo>&);
void Permutation (vector<Individuo>&);
int Pbc(int);
bool comp_Alleli (Allele, Allele);
void Ordina_Alleli (vector<Allele>&);
void Stampa_best(vector<Individuo>, int);	
void Stampa_media(vector<Individuo>, int);
void Pre_Selezione(vector<Individuo>);
void Set_Elitario(vector<Individuo>);
void Inserisci_Elitario(vector<Individuo>&);
void Print_percorso(vector<Individuo>, double*, double*);
void Input();
using namespace std;
 
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

//----------------------------PROBLEMA DEL COMMESSO VIAGGIATORE-----------------------------------------
	
	Input();
	system ("rm best_mean.dat best_path.dat percorso_migliore.dat");
	//----------------città distribuite su un cerchio----------------------------
	double* x_citta=new double [n_citta];
	double* y_citta=new double [n_citta];
	if(forma==0){		//citta distribuite sul cerchio se forma=0
		for (int i=0; i<n_citta; i++){
		double angolo=rnd.Rannyu(0,2*M_PI);
		x_citta[i]=cos(angolo);
		y_citta[i]=sin(angolo);
		}
	}
	if (forma==1){		//citta distribuite nel quadrato se forma=1
		for (int i=0; i<n_citta; i++){
		x_citta[i]=rnd.Rannyu(0,1);
		y_citta[i]=rnd.Rannyu(0,1);
		}
	}
	if(forma!=0 && forma!=1) {cout<<"ATTENZIONE! il secondo campo del file di input può prendere solo zero o uno!"<<endl;}	
	double inserimento_elitario;
	selection_array=new double [n_citta*n_citta];
	vector<Individuo> popolazione(n_citta*n_citta);
	Inizializza_p(popolazione);
	vector<Individuo> popolazione_appoggio(n_citta*n_citta);
	popolazione_appoggio=popolazione; 
	Elitario.Inizializza(n_citta);
	Elitario.Calcola_lunghezza(x_citta, y_citta);
	double check_popolazione=Check(popolazione);
	Calcola_lunghezza(popolazione, x_citta, y_citta);
	Ordina(popolazione);
	double mutazione;
	
	for (int i=0; i<n_generazioni; i++){
		
		check_popolazione=Check(popolazione);
		cout<<i+1<<")"<<"check popolazione= "<<check_popolazione<<endl;			
		popolazione_appoggio=popolazione;
		Pre_Selezione(popolazione_appoggio);
		Crossover(popolazione, popolazione_appoggio); //crossover (avviene nel 70% dei casi)
		
		for (int j=0; j<n_citta*n_citta; j++){    //mutazioni
			mutazione=rnd.Rannyu(0,1);
			if(mutazione<p_inversion){Inversion(popolazione);}
			if(mutazione>0.5 and mutazione<0.5+p_pair_mutation){Pair_Mutation(popolazione);}
			if(mutazione>0.1 and mutazione<0.1+p_simple_shift){Simple_shift(popolazione);}
			if(mutazione>0.3 and mutazione<0.3+p_shift){Shift(popolazione);}
			if(mutazione>0.7 and mutazione<0.7+p_permutation){Permutation(popolazione);}
		}
		Calcola_lunghezza(popolazione, x_citta, y_citta);
		Ordina(popolazione);	
		Set_Elitario(popolazione);
		inserimento_elitario=rnd.Rannyu(0,1);
		if(inserimento_elitario<=p_elitario) {
			Inserisci_Elitario(popolazione);	
			Ordina(popolazione);
		}	
		Stampa_best(popolazione,i);
		Stampa_media(popolazione,i);
	}
	
	Print_percorso(popolazione, x_citta, y_citta);
		
	delete [] x_citta;
	delete [] y_citta;
	delete [] selection_array;

//--------------------------------------------------------------------------------------------------------	
   rnd.SaveSeed();
   return 0;
}
void Input(){
	ifstream in;
	in.open("input.dat");
	in>>n_citta;
	in>>forma;
	in>>n_generazioni;
	in>>p_crossover;
	in>>p_inversion;
	in>>p_pair_mutation;
	in>>p_simple_shift;
	in>>p_shift;
	in>>p_permutation;
	in>>p_elitario;
}

bool comp(Individuo a, Individuo b) {
	return (a.Get_lunghezza()<b.Get_lunghezza());
}
bool comp_Alleli(Allele a, Allele b) {
	return (a.Get_indice()<b.Get_indice());
}
void Ordina (vector<Individuo>& popolazione){
	sort(popolazione.begin(), popolazione.end(), comp);
}
void Inizializza_p(vector<Individuo>& popolazione){
	for (int i=0; i<(int)(popolazione.size()); i++){
		//Individuo percorso (n_citta);
		//percorso.Inizializza(n_citta);
		popolazione[i].Inizializza(n_citta);
	}
}
void Calcola_lunghezza(vector<Individuo>& popolazione, double* x, double*y){
	for (int i=0; i<(int)(popolazione.size()); i++){
		popolazione[i].Calcola_lunghezza(x,y);
	}
}
double Check(vector<Individuo> popolazione){
	double check=0;
	for(int i=0; i<(int)(popolazione.size()); i++){
		check+=popolazione[i].Check();
	}
	if(check!=0) {
		cout<<"gli individui non rispettano i constrain!"<<endl;
		return 0;
	}
	else {cout<<"gli individui rispettano i constrain"<<endl;
	return check;
	}
}
int Pbc(int indice){
	int appo=indice-n_citta;
	if (indice<n_citta) {return indice;}	
	else {
		do {appo=appo-n_citta;}
		while(appo>=n_citta);
		return appo;
	}
}
void Pre_Selezione(vector <Individuo> popolazione) {
	normalizzazione=0;
	//cout<<"inizializzo normalizzazione"<<endl;
	for(int i=0; i<(int)popolazione.size(); i++){
		normalizzazione+=1/(double)(popolazione[i].Get_lunghezza());
	}
	//cout<<"calcolo normalizzazione"<<endl;
	selection_array[0]=(1.0/(double)(popolazione[0].Get_lunghezza()))/(double)(normalizzazione);
	for(int i=1;i<(int)popolazione.size();i++){
		selection_array[i]=selection_array[i-1]+(1.0/(double)(popolazione[i].Get_lunghezza()))/(double)(normalizzazione);
	}
}	
int Selection(vector<Individuo> popolazione){
	
	number=rnd.Rannyu(0,1);
	selected_indice=0;
	if (number<=selection_array[0]) {selected_indice=0;}
	for(int i=1; i<(int)popolazione.size(); i++){
		if (number>selection_array[i-1] && number<selection_array[i]) {selected_indice=i;}
	}
	return selected_indice;	
}
void Crossover(vector<Individuo>& popolazione, vector<Individuo> popolazione_appoggio){ //vorrei farla in modo che ricambi TUTTA la generazion	
	for (int p=0; p<(int)popolazione.size()-1; p+=2) {
		
		indice_padre=Selection(popolazione_appoggio);
		indice_madre=Selection(popolazione_appoggio);
		double estrazione=rnd.Rannyu(0,1);
		if (estrazione>p_crossover) {	
			popolazione[p]=popolazione_appoggio[indice_padre];
			popolazione[p+1]=popolazione_appoggio[indice_madre];
		}

		if (estrazione<=p_crossover){

			indice_taglio=(int)(rnd.Rannyu(0,n_citta)); //per scegliere dove tagliare
	
			for(int i=0; i<=indice_taglio; i++){
				popolazione[p].Set_Elemento(popolazione_appoggio[indice_padre].Get_Elemento(i),i);
				popolazione[p+1].Set_Elemento(popolazione_appoggio[indice_madre].Get_Elemento(i),i);	
			}
			vector<Allele> geni_padre;  //salvo i geni mancanti del padre con associato l'indice in cui compaiono nella madre
			geni_padre.resize(n_citta-(indice_taglio+1));
			vector<Allele> geni_madre; //salvo i geni mancanti della madre con associato l'indice in cui compaiono nel padre
			geni_madre.resize(n_citta-(indice_taglio+1));
			int k=0;
			for (int i=indice_taglio+1; i<n_citta; i++){
				geni_padre[k].Set_citta(popolazione_appoggio[indice_padre].Get_Elemento(i));
				geni_madre[k].Set_citta(popolazione_appoggio[indice_madre].Get_Elemento(i));
				for (int j=0; j<n_citta; j++){
					if(popolazione_appoggio[indice_padre].Get_Elemento(i)==popolazione_appoggio[indice_madre].Get_Elemento(j)) {
						geni_padre[k].Set_indice(j);		
					}
					if(popolazione_appoggio[indice_madre].Get_Elemento(i)==popolazione_appoggio[indice_padre].Get_Elemento(j)){
						geni_madre[k].Set_indice(j);
					}
				}		
				k++;
			}
			Ordina_Alleli(geni_padre);
			Ordina_Alleli(geni_madre);
			k=0;
			for(int i=indice_taglio+1; i<n_citta; i++){
				popolazione[p].Set_Elemento(geni_padre[k].Get_citta(),i);
				popolazione[p+1].Set_Elemento(geni_madre[k].Get_citta(),i);
				k++;
			}
		}
	}
	return;
}
	

void Pair_Mutation(vector<Individuo>& popolazione){ //estrae gli individui a cui fare la mutazione
	int individuo=(int)(rnd.Rannyu(0,popolazione.size()));
	int indice_1=(int)(rnd.Rannyu(0,n_citta));
	int indice_2=(int)(rnd.Rannyu(0,n_citta));
	int appo=popolazione[individuo].Get_Elemento(indice_1);
	popolazione[individuo].Set_Elemento(popolazione[individuo].Get_Elemento(indice_2),indice_1);
	popolazione[individuo].Set_Elemento(appo, indice_2);	
}
void Simple_shift(vector<Individuo>& popolazione){
	int individuo=(int)(rnd.Rannyu(0,popolazione.size()));	
	int n=(int)(rnd.Rannyu(0,n_citta)); //l'elemento n-esimo diventa il primo del percorso
	rotate(popolazione[individuo].Get_start(),popolazione[individuo].Get_start()+n,popolazione[individuo].Get_end());
}
void Shift(vector<Individuo>& popolazione){
	int individuo=(int)(rnd.Rannyu(0,popolazione.size()));
	int posizione_1=(int)(rnd.Rannyu(0,n_citta/2));  //primo elemento del segmento da spostare
	int lunghezza_segmento=(int)(rnd.Rannyu(0,n_citta/3));
	int posizione_2=(int)(rnd.Rannyu(Pbc(posizione_1+lunghezza_segmento+1),n_citta));  //posizione dove si troverà l'ultimo elemento del segmento spostato
	rotate(popolazione[individuo].Get_start()+posizione_1, popolazione[individuo].Get_start()+Pbc(posizione_2-lunghezza_segmento), popolazione[individuo].Get_start()+Pbc(posizione_2));
}
void Permutation (vector<Individuo>& popolazione){
	
	int individuo=(int)(rnd.Rannyu(0,popolazione.size()));
	int estremo_1=(int)(rnd.Rannyu(0,n_citta/2));
	int lunghezza_segmento=(int)(rnd.Rannyu(0,Pbc(n_citta/2 -estremo_1)));	        
	int estremo_2=(int)(rnd.Rannyu(n_citta/2, Pbc(n_citta-lunghezza_segmento)));
	swap_ranges(popolazione[individuo].Get_start()+(estremo_1), popolazione[individuo].Get_start()+estremo_1+lunghezza_segmento, popolazione[individuo].Get_start()+estremo_2);
	
}
void Inversion (vector<Individuo>& popolazione) {
	int individuo=(int)(rnd.Rannyu(0,popolazione.size()));
	int inizio=(int)(rnd.Rannyu(0,n_citta));
	int quanti=(int)(rnd.Rannyu(0,n_citta)); //lunghezza del segmento da invertire
	reverse(popolazione[individuo].Get_start()+inizio, popolazione[individuo].Get_start()+Pbc(inizio+quanti));	
}
void Ordina_Alleli (vector<Allele>& vettore){
	sort(vettore.begin(), vettore.end(), comp_Alleli);
}

void Stampa_best(vector<Individuo> popolazione, int indice){
	ofstream out;
	out.open("best_path.dat", ios::app);
	out<<indice<<"	"<<popolazione[0].Get_lunghezza()<<endl;
}

void Stampa_media(vector<Individuo> popolazione, int indice){
	ofstream out;
	out.open("best_mean.dat", ios::app);	
	double media=0;
	for (int i=0; i<(int)popolazione.size()/((int)2); i++){
		media+=popolazione[i].Get_lunghezza();
	}
	out<<indice<<"	"<<media*2/(double)popolazione.size()<<endl;
}
void Set_Elitario(vector<Individuo> popolazione){  //uso questo metodo quando la popolazione è già stata ordinata
	if(popolazione[0].Get_lunghezza()<Elitario.Get_lunghezza()){
		Elitario=popolazione[0];
	}
}
void Inserisci_Elitario(vector<Individuo>& popolazione){
	popolazione[n_citta-1]=Elitario;
}	
void Print_percorso(vector<Individuo> popolazione, double* x, double* y){
	ofstream out;
	out.open("percorso_migliore.dat");
	for (int i=0; i<n_citta; i++){
		out<<x[popolazione[0].Get_Elemento(i)]<<"	"<<y[popolazione[0].Get_Elemento(i)]<<endl;
	}
	out.close();
}			
