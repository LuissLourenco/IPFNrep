#include "../../src/DataAnalysis.cpp"
#include <dirent.h>
#include <string>
#include <vector>
#include "TLine.h"

using namespace std;

typedef struct simulation_data{

	double px0 = -10;
	double kdamp = 0;
	int T = 5000;
	int N = 1500000;
	int pri = 10;
	int wave_type = 3;
	double tfwhm = 50;
	double stable = 1000000000;
	double Eo = 50;
	double delta = 0;
	double w0 = 5;
	double lambda = 1;
	int l = 0;
	int p = 0;

	double rmin = 0;
	double rmax = 4.1;
	double dr = 0.1;

	double phimin = 0;
	double phimax = 181;
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
				sprintf(cmd, "gnome-terminal -- bash -ic 'xdotool getactivewindow windowminimize; ./single_run %i %.14e %.14e %.14e %.14e %i %i %i %i %.14e %.14e %.14e %.14e %.14e %.14e %i %i ; rm InputToBatch%i.txt'", 
					process, r, phi, px0, kdamp, T, N, pri, wave_type, tfwhm, stable , Eo, delta, w0, lambda, l, p, process);
				system(cmd);
			}
			process++;
			r+=dr;
		}
		phi+=dphi;
	}

	system("sleep 10");
	
	sprintf(cmd, "rm -rvf %s", data.directory.c_str());
	system(cmd);
	sprintf(cmd, "mkdir %s", data.directory.c_str());
	system(cmd);
	sprintf(cmd, "mv Out* %s", data.directory.c_str());
	system(cmd);
	
}


int main(){
/*
	// A0 = 70
	simulation_data sim1;
	sim1.Eo = 70;
	sim1.phimax = 91;
	sim1.directory = "../outputs/Data_a70/";
	sim1.print();
	run_simulations(sim1, 3);

	// px = -20
	simulation_data sim2;
	sim2.px = -20;
	sim2.phimax = 91;
	sim2.directory = "../outputs/Data_px-20/";
	sim2.print();
	run_simulations(sim2, 3);
*/
	// p = 0, l = 2
	simulation_data sim3;
	sim3.p = 0;
	sim3.l = 2;
	sim3.directory = "../outputs/Data_p0l2/";
	sim3.print();
	run_simulations(sim3, 3);



}