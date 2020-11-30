#include "../src/DataAnalysis.cpp"
#include <dirent.h>
#include <string>
#include <vector>
#include "TLine.h"

using namespace std;

typedef struct simulation_data{

	double r = 2;
	double phi = 0;
	double px0 = -20;
	double kdamp = 0;
	int T = 500;
	int N = 100000;
	int pri = 10;
	int wave_type = 3;
	double tfwhm = 50;
	double stable = 100000;
	double Eo = 10;
	double delta = 0;
	double w0 = 5;
	double lambda = 1;
	int l = 1;
	int p = 0;

	void print(){
		char cmd[1024];
		sprintf(cmd, "Sim: r=%1.2lf; \u03A6=%1.2lf; p=%1.2lf; E=%1.2lf; T=%2d; N=%2d; pri=%2d; kdamp=%1.2lf; l=%2d; p=%2d; rise=%1.2lf; stable=%1.2lf; w0=%1.2lf; lambda=%1.2lf; wt=%1d", r, phi, px0, Eo, T, N, pri, kdamp, l, p, tfwhm, stable, w0, lambda, wave_type);
		cout << cmd << endl;
	}
	void print_short(){
		char cmd[1024];
		sprintf(cmd, "| r %1.2lf| \u03A6 %1.2lf| p %1.2lf| E %1.2lf| T %2d| N %2d| \u03BC %1.2lf| lp %1d%1d| \u21E7 %4.0lf| \u21E8 %4.0lf", r, phi, px0, Eo, T, N, kdamp, l, p, tfwhm, stable);
		cout << cmd << endl;
	}


}simulation_data;

void run_simulations(vector<simulation_data> data, string directory, int n_terminals){

	system("g++ -o2 single_run.cpp -lm -o single_run `root-config --cflags --glibs`");
	
	double r, phi, px0, kdamp, tfwhm, stable, Eo, delta, w0, lambda;
	int T, N, pri, wave_type, l, p;
	
	char cmd[512];

	FILE *data_logger;

	data_logger = fopen("Logger","w");
	
	fprintf(data_logger, "process\tr\tphi\tpx0\tkdamp\tT\tN\tpri\twave_type\ttfwhm\tstable\tEo\tdelta\tw0\tlambda\tl\tp");

	for(int process = 0; process < data.size(); process++){

		r = data[process].r;
		phi = data[process].phi;
		px0 = data[process].px0;
		kdamp = data[process].kdamp;
		T = data[process].T;
		N = data[process].N;
		pri = data[process].pri;
		wave_type = data[process].wave_type;
		tfwhm = data[process].tfwhm;
		stable = data[process].stable;
		Eo = data[process].Eo;
		delta = data[process].delta;
		w0 = data[process].w0;
		lambda = data[process].lambda;
		l = data[process].l;
		p = data[process].p;

		fprintf(data_logger, "\n%05d\t", process);
		fprintf(data_logger, "%.5e\t", r);
		fprintf(data_logger, "%.5e\t", phi);
		fprintf(data_logger, "%.5e\t", px0);
		fprintf(data_logger, "%.5e\t", kdamp);
		fprintf(data_logger, "%d\t", T);
		fprintf(data_logger, "%d\t", N);
		fprintf(data_logger, "%d\t", pri);
		fprintf(data_logger, "%d\t", wave_type);
		fprintf(data_logger, "%.5e\t", tfwhm);
		fprintf(data_logger, "%.5e\t", stable);
		fprintf(data_logger, "%.5e\t", Eo);
		fprintf(data_logger, "%.5e\t", delta);
		fprintf(data_logger, "%.5e\t", w0);
		fprintf(data_logger, "%.5e\t", lambda);
		fprintf(data_logger, "%d\t", l);
		fprintf(data_logger, "%d", p);

		sprintf(cmd, "%4d of %4ld \u262d  ", process+1, data.size());
		cout << cmd;
		data[process].print_short();

		if(process % n_terminals == 0){
			sprintf(cmd, "./single_run %i %.14e %.14e %.14e %.14e %i %i %i %i %.14e %.14e %.14e %.14e %.14e %.14e %i %i ; rm InputToBatch%i.txt", 
				process, r, phi, px0, kdamp, T, N, pri, wave_type, tfwhm, stable , Eo, delta, w0, lambda, l, p, process);
			system(cmd);
		}else{
			sprintf(cmd, "gnome-terminal --tab -- bash -ic './single_run %i %.14e %.14e %.14e %.14e %i %i %i %i %.14e %.14e %.14e %.14e %.14e %.14e %i %i ; rm InputToBatch%i.txt'", 
			//sprintf(cmd, "gnome-terminal -- bash -ic 'xdotool getactivewindow windowminimize; ./single_run %i %.14e %.14e %.14e %.14e %i %i %i %i %.14e %.14e %.14e %.14e %.14e %.14e %i %i ; rm InputToBatch%i.txt'", 
				process, r, phi, px0, kdamp, T, N, pri, wave_type, tfwhm, stable , Eo, delta, w0, lambda, l, p, process);
			system(cmd);
		}
	}

	fclose(data_logger);

	cout << "Sleeping!" << endl;
	system("sleep 10");
	//system("read line");
	
	sprintf(cmd, "rm %s/*", directory.c_str());
	system(cmd);
	sprintf(cmd, "mkdir %s", directory.c_str());
	system(cmd);
	sprintf(cmd, "mv Out* %s", directory.c_str());
	system(cmd);
	sprintf(cmd, "mv Logger %s", directory.c_str());
	system(cmd);
	
}

int main(){

	vector<simulation_data> data;
	simulation_data single;

	single.r = 2;
	single.phi = 0;
	single.px0 = -20;
	single.kdamp = 0;
	single.T = 1000;
	single.N = 1000000;
	single.pri = 20;
	single.wave_type = 3;
	single.tfwhm = 50;
	single.stable = 100000;
	single.Eo = 10;
	single.delta = 0;
	single.w0 = 5;
	single.lambda = 1;
	single.l = 1;
	single.p = 0;

	double px0_min = -10;
	double px0_max = -0.2;
	double dpx0 = 0.2;

	double a0_min = 0.2;
	double a0_max = 20;
	double da0 = 0.2;

	double a0, px0;

	px0 = px0_min;
	while(px0 <= px0_max){
		a0 = -1.5*px0;
		while(a0 <= a0_max){
			single.Eo = a0;
			single.px0 = px0;
			data.push_back(single);
			a0 += da0;
		}
		px0 += dpx0;
	}

	run_simulations(data, "../outputs/Data_Finner", 6);

	return 0;
}