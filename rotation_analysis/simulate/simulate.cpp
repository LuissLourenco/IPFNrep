#include "../src/DataAnalysis.cpp"
#include <dirent.h>
#include <string>
#include <vector>
#include "TLine.h"

using namespace std;

double dt=2e-3;
bool NHECOS = true;

typedef struct simulation_data{

	double px0 = -20;
	double kdamp = 0;
	int T = 20e3;
	int N = T/dt;
	int pri = 10;
	int wave_type = 3;
	double tfwhm = 50;
	double stable = T;
	double Eo = 10;
	double delta = 0;
	double w0 = 5;
	double lambda = 1;
	int l = 1;
	int p = 0;

	double rmin = 1;
	double rmax = 4;
	double dr = 0.1;

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
	if(NHECOS) system("sleep 30");
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

	simulation_data sim;
	
	/*for(int i=1; i<=2; i++){
		for(int j=3; j<=3; j++){
			sim.Eo = (double)(5*i);
			sim.px0 = -(double)(15+j*5);
			char aux[128];
			sprintf(aux, "../outputs/Data01_a0_%02.lf_p0_%02.lf/",sim.Eo, -sim.px0);
			sim.directory = aux;
			run_simulations(sim,5);
		}
	}*/

	NHECOS = false;
	
	sim.rmin = 3;
	sim.rmax = 3.01;
	sim.phimin = 0;
	sim.phimax = 1;

	sim.Eo = 30;
	sim.px0 = -20;

	sim.directory = "../outputs/meme/";
	run_simulations(sim,1);





	
	
	return 0;
}