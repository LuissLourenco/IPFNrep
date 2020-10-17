#include "../src/DataAnalysis.cpp"
#include <dirent.h>
#include <string>
#include <vector>
#include "TLine.h"

using namespace std;

double dt=2e-3;

typedef struct simulation_data{

	double px0 = -10;
	double kdamp = 0;
	int T = 5e3;
	int N = T/dt;
	int pri = 10;
	int wave_type = 3;
	double tfwhm = 50;
	double stable = T*10;
	double Eo = 10;
	double delta = 0;
	double w0 = 5;
	double lambda = 1;
	int l = 1;
	int p = 0;

	double rmin = 1;
	double rmax = 2.1;
	double dr = 1;

	double phimin = 0;
	double phimax = 180;
	double dphi = 5;

	string directory = "../outputs/Data";

	void print(){
		cout << endl << endl << endl;
		cout << "Simulation <" << directory << ">" << endl;
		cout << "r from " << rmin << " to " << rmax << " with increment " << dr << endl;
		cout << "phi from " << phimin << " to " << phimax << " with increment " << dphi << endl;
		cout << "px0        = "<< px0       << endl;
		cout << "kdamp      = "<< kdamp     << endl;
		cout << "T          = "<< T         << endl;
		cout << "N          = "<< N         << endl;
		cout << "pri        = "<< pri       << endl;
		cout << "wave_type  = "<< wave_type << endl;
		cout << "tfwhm      = "<< tfwhm     << endl;
		cout << "stable     = "<< stable    << endl;
		cout << "Eo         = "<< Eo        << endl;
		cout << "delta      = "<< delta     << endl;
		cout << "w0         = "<< w0        << endl;
		cout << "lambda     = "<< lambda    << endl;
		cout << "l          = "<< l         << endl;
		cout << "p          = "<< p         << endl;
		cout << endl << endl << endl;
	}

}simulation_data;

void run_simulations(simulation_data data, int n_terminals){

	system("g++ -o2 single_run.cpp -lm -o single_run `root-config --cflags --glibs`");

	int process = 0;
	
	double r;
	double phi;
	double px0 = data.px0;
	double kdamp = data.kdamp;
	int T = data.T;
	int N = data.N;
	int pri = data.pri;
	int wave_type = data.wave_type;
	double tfwhm = data.tfwhm;
	double stable = data.stable;
	double Eo = data.Eo;
	double delta = data.delta;
	double w0 = data.w0;
	double lambda = data.lambda;
	int l = data.l;
	int p = data.p;

	char cmd[512];

	double rmin = data.rmin;
	double rmax = data.rmax;
	double dr = data.dr;

	double phimin = data.phimin;
	double phimax = data.phimax;
	double dphi = data.dphi;
	
	phi = phimin;
	while(phi < phimax){
		r = rmin;
		while(r < rmax){
			if(process % n_terminals == 0){
				sprintf(cmd, "./single_run %i %.14e %.14e %.14e %.14e %i %i %i %i %.14e %.14e %.14e %.14e %.14e %.14e %i %i ; rm InputToBatch%i.txt", 
					process, r, phi, px0, kdamp, T, N, pri, wave_type, tfwhm, stable , Eo, delta, w0, lambda, l, p, process);
				system(cmd);
				cout << "PROCESS = " << process << "\tr = " << r << "\tphi = " << phi << endl;
			}else{
				sprintf(cmd, "gnome-terminal --tab -- bash -ic './single_run %i %.14e %.14e %.14e %.14e %i %i %i %i %.14e %.14e %.14e %.14e %.14e %.14e %i %i ; rm InputToBatch%i.txt'", 
				//sprintf(cmd, "gnome-terminal -- bash -ic 'xdotool getactivewindow windowminimize; ./single_run %i %.14e %.14e %.14e %.14e %i %i %i %i %.14e %.14e %.14e %.14e %.14e %.14e %i %i ; rm InputToBatch%i.txt'", 
					process, r, phi, px0, kdamp, T, N, pri, wave_type, tfwhm, stable , Eo, delta, w0, lambda, l, p, process);
				system(cmd);
			}
			process++;
			r+=dr;
		}
		phi+=dphi;
	}

	cout << "Sleeping!" << endl;
	system("sleep 30");
	//system("read line");
	
	sprintf(cmd, "rm %s/*", data.directory.c_str());
	system(cmd);
	sprintf(cmd, "mkdir %s", data.directory.c_str());
	system(cmd);
	sprintf(cmd, "mv Out* %s", data.directory.c_str());
	system(cmd);
	
}


int main(){

	/*
	CHANGE THE PARAMETERS OF THE SIMULATION HERE
	struct simulation_data HAS THE PARAMETERS
	run_simulations(data, n_terminals)
	n_terminals IS THE NUMBER OS TABS THAT WILL BE RAN AT THE SAME TIME
	*/

	

	/*
	
	simulation_data sim1;

	//a0=5,10,15,20,25,30    
	//px0=20,25,30
	//pl=01
	for(int i=1; i<=6; i++){
		for(int j=1; j<=3; j++){

			sim1.Eo = 5.*i;
			sim1.px0 = -15. - j*5.;
			char aux[128];
			sprintf(aux, "../outputs/PicardLindelof_a0_%02.0lf_px0_%02.0lf/", sim1.Eo, sim1.px0);
			sim1.directory = aux;
			sim1.print();
			run_simulations(sim1,6);
		}
	}

	//a0=10   px0=-10   pl=02
	sim1.Eo = 10;
	sim1.px0 = -10;
	sim1.l = 2;
	sim1.phimax = 181;
	char aux[128];
	sprintf(aux, "../outputs/PicardLindelof02_a0_%02.0lf_px0_%02.0lf/", sim1.Eo, sim1.px0);
	sim1.directory = aux;
	sim1.print();
	run_simulations(sim1,6);



	//a0=5(baixo)   px0=-5   pl=01
	sim1.Eo = 5;
	sim1.px0 = -5;
	sim1.l = 1;
	sim1.phimax = 181;
	char aux1[128];
	sprintf(aux1, "../outputs/PicardLindelof01_a0_%02.0lf_px0_%02.0lf/", sim1.Eo, sim1.px0);
	sim1.directory = aux1;
	sim1.print();
	run_simulations(sim1,6);

	*/
	

	simulation_data sim1;
	sim1.Eo = 5;
	sim1.Eo = -10;
	sim1.directory = "../outputs/Testephi0_a0_5_p0_10/";
	run_simulations(sim1,6);
}