#include "DataAnalysis.cpp"
#include <dirent.h>
#include <string>
#include <vector>
#include "TLine.h"

using namespace std;


void run_simulations(double a0_in, double T_in, int N_in){

	system("g++ -o2 macro_laguerre.cpp -lm -o macro_laguerre `root-config --cflags --glibs`");

	int process = 0;
	
	double r = 0.1;
	double phi = 15;
	double px0 = -10;
	double kdamp = 0;
	int T = T_in;
	int N = N_in;
	int pri = 10;
	int wave_type = 3;
	double tfwhm = 50;
	double stable = 1000000000;
	double Eo = a0_in;
	double delta = 0;
	double w0 = 5;
	double lambda = 1;
	int l = 0;
	int p = 0;

	char cmd[512];

	double rmin = 0;
	double rmax = 4;
	double dr = 0.1;

	double phimin = 0;
	double phimax = 91;
	double dphi = 5;

	phi = phimin;
	while(phi < phimax){
		r = rmin;
		while(r < rmax){
			if(process % 3 == 0){
				sprintf(cmd, "./macro_laguerre %i %.14e %.14e %.14e %.14e %i %i %i %i %.14e %.14e %.14e %.14e %.14e %.14e %i %i ; rm InputToBatch%i.txt", 
					process, r, phi, px0, kdamp, T, N, pri, wave_type, tfwhm, stable , Eo, delta, w0, lambda, l, p, process);
				system(cmd);
				cout << "PROCESS = " << process << "\tr = " << r << "\tphi = " << phi << endl;
			}else{
				sprintf(cmd, "gnome-terminal -- bash -ic 'xdotool getactivewindow windowminimize; ./macro_laguerre %i %.14e %.14e %.14e %.14e %i %i %i %i %.14e %.14e %.14e %.14e %.14e %.14e %i %i ; rm InputToBatch%i.txt'", 
					process, r, phi, px0, kdamp, T, N, pri, wave_type, tfwhm, stable , Eo, delta, w0, lambda, l, p, process);
				system(cmd);
			}
			process++;
			r+=dr;
		}
		phi+=dphi;
	}

	system("sleep 10");

	sprintf(cmd, "rm -rf Data_p0l0_%i/", Eo);
	system(cmd);
	sprintf(cmd, "mkdir Data_p0l0_%i/", Eo);
	system(cmd);
	sprintf(cmd, "mv Out* Data_p0l0_%i/", Eo);
	system(cmd);
	
}

void mood(string* files){

	int n_files = 1;

	int* n_cols = new int[n_files];
	int* n_points = new int[n_files];
	double** values;

	DataSet* T = new DataSet[n_files];
	DataSet* X = new DataSet[n_files];
	DataSet* Y = new DataSet[n_files];
	DataSet* Z = new DataSet[n_files];
	DataSet* PX = new DataSet[n_files];
	DataSet* PY = new DataSet[n_files];
	DataSet* PZ = new DataSet[n_files];
	DataSet* GAMMA = new DataSet[n_files];
	DataSet* THETA = new DataSet[n_files];
	DataSet* PTHETA = new DataSet[n_files];

	for(int i = 0; i < n_files; i++){
		values = ReadFile(files[i], &n_cols[i], &n_points[i], false, false);
		T[i] = DataSet(n_points[i], values[0]).compress(100);
		Y[i] = DataSet(1, values[2]);
		Z[i] = DataSet(1, values[3]);
		THETA[i] = DataSet(n_points[i], values[8]).compress(100);

		cout << files[i] << " " << atan2(Z[i][0].val(), Y[i][0].val()) << " " << sqrt(Y[i][0]*Y[i][0]+Z[i][0]*Z[i][0]).val() << endl;

	}


	double* r = new double[n_files];
	double* phi = new double[n_files];
	double* period = new double[n_files];
	double* lim1 = new double[n_files];
	double* lim2 = new double[n_files];
	for(int i = 0; i < n_files; i++){
		Var amax = THETA[i].getMax();
		Var amin = THETA[i].getMin();
		Var delta = (amax - amin) / Var(96);
		bool found_max = false;
		bool found_min = false;
		int j1, j2;
		for(int j = 0; j < THETA[i].size(); j++){
			if(abs(THETA[i][j] - amax) <= delta && !found_max){ 
				found_max = true;
				j1 = j;
			}
			if(!found_max) continue;
			if(abs(THETA[i][j] - amin) <= delta && !found_min) found_min = true;
			if(!found_min) continue;
			if(abs(THETA[i][j] - amax) <= delta){ 
				j2 = j;
				break;
			}
		}
		period[i] = 2*abs(T[i][j2].val() - T[i][j1].val());
		cout << files[i] << " " << i << " " << period[i] << " " << j1 << " " << j2 << endl;
		lim1[i] = T[i][j1].val();
		lim2[i] = T[i][j2].val();

		r[i] = sqrt(Y[i][0]*Y[i][0]+Z[i][0]*Z[i][0]).val();
		phi[i] = atan2(Z[i][0].val(), Y[i][0].val());
	}






	//=====================================================================

	TCanvas* c1 = new TCanvas("c1", "", 1500, 1000);
	c1->cd();
	
	/*
	TGraph* graph_period = new TGraph(n_files, r, period);
	graph_period->SetTitle("Period of Rotation as Function of r_{0}, #phi_{0}=22.5, px_{0}=-10, a_{0}=30;r_{0};T");
	graph_period->Draw("ACP");

	c1->SaveAs("plotCpp.png");
	*/
	
	TGraph* g; char out[16];
	TLine* line1; TLine* line2;
	for(int i = 0; i < n_files; i++){
		sprintf(out, "plotCpp%s.png", files[i].substr(5).c_str());
		g = GetTGraph(T[i], THETA[i]);
		g->SetTitle(files[i].c_str());

		line1 = new TLine(lim1[i], 0, lim1[i], 3);
		line2 = new TLine(lim2[i], 0, lim2[i], 3);
		line1->SetLineColor(2);
		line1->SetLineWidth(2);
		line2->SetLineColor(2);
		line2->SetLineWidth(2);

		g->Draw("ALP");
		line1->Draw("SAME");
		line2->Draw("SAME");
		c1->SaveAs(out);
	}

	delete[] r;
	delete[] phi;
	delete[] period;
	delete[] lim1;
	delete[] lim2;
	delete c1;
	delete g;
	delete line1;
	delete line2;	

}

int main(){

	//run_simulations(100,  1500,  1500000);
	//run_simulations( 50,  5000, 1500000);
	run_simulations( 30, 1000,  200000);

/*
	int n_files;
	//string* files = list_dir("Data_px-10_phi0_a30/", &n_files);
	//string* files = list_dir("Data_px-10_phi5_a30/", &n_files);
	//string* files = list_dir("Data_px-10_phi10_a30/", &n_files);
	//string* files = list_dir("Data_px-10_phi15_a30/", &n_files);
	//string* files = list_dir("Data_px-10_phi22.5_a30/", &n_files);
	//string* files = list_dir("Data_px-10_phi45_a30/", &n_files);
	string* files = list_dir("Data/", &n_files);

	
	int* n_cols = new int[n_files];
	int* n_points = new int[n_files];
	double** values;

	DataSet* T = new DataSet[n_files];
	DataSet* X = new DataSet[n_files];
	DataSet* Y = new DataSet[n_files];
	DataSet* Z = new DataSet[n_files];
	DataSet* PX = new DataSet[n_files];
	DataSet* PY = new DataSet[n_files];
	DataSet* PZ = new DataSet[n_files];
	DataSet* GAMMA = new DataSet[n_files];
	DataSet* THETA = new DataSet[n_files];
	DataSet* PTHETA = new DataSet[n_files];

	for(int i = 0; i < n_files; i++){
		values = ReadFile(files[i], &n_cols[i], &n_points[i], false, false);
		T[i] = DataSet(n_points[i], values[0]).compress(100);
		Y[i] = DataSet(1, values[2]);
		Z[i] = DataSet(1, values[3]);
		THETA[i] = DataSet(n_points[i], values[8]).compress(100);

		cout << files[i] << " " << atan2(Z[i][0].val(), Y[i][0].val()) << " " << sqrt(Y[i][0]*Y[i][0]+Z[i][0]*Z[i][0]).val() << endl;

	}


	double* r = new double[n_files];
	double* phi = new double[n_files];
	double* period = new double[n_files];
	double* lim1 = new double[n_files];
	double* lim2 = new double[n_files];
	for(int i = 0; i < n_files; i++){
		Var amax = THETA[i].getMax();
		Var amin = THETA[i].getMin();
		Var delta = (amax - amin) / Var(96);
		bool found_max = false;
		bool found_min = false;
		int j1, j2;
		for(int j = 0; j < THETA[i].size(); j++){
			if(abs(THETA[i][j] - amax) <= delta && !found_max){ 
				found_max = true;
				j1 = j;
			}
			if(!found_max) continue;
			if(abs(THETA[i][j] - amin) <= delta && !found_min) found_min = true;
			if(!found_min) continue;
			if(abs(THETA[i][j] - amax) <= delta){ 
				j2 = j;
				break;
			}
		}
		period[i] = 2*abs(T[i][j2].val() - T[i][j1].val());
		cout << files[i] << " " << i << " " << period[i] << " " << j1 << " " << j2 << endl;
		lim1[i] = T[i][j1].val();
		lim2[i] = T[i][j2].val();

		r[i] = sqrt(Y[i][0]*Y[i][0]+Z[i][0]*Z[i][0]).val();
		phi[i] = atan2(Z[i][0].val(), Y[i][0].val());
	}






	//=====================================================================

	TCanvas* c1 = new TCanvas("c1", "", 1500, 1000);
	c1->cd();
	
	
	TGraph* graph_period = new TGraph(n_files, r, period);
	graph_period->SetTitle("Period of Rotation as Function of r_{0}, #phi_{0}=22.5, px_{0}=-10, a_{0}=30;r_{0};T");
	graph_period->Draw("ACP");

	c1->SaveAs("plotCpp.png");
	
	
	TGraph* g; char out[16];
	TLine* line1; TLine* line2;
	for(int i = 0; i < n_files; i++){
		sprintf(out, "plotCpp%s.png", files[i].substr(5).c_str());
		g = GetTGraph(T[i], THETA[i]);
		g->SetTitle(files[i].c_str());

		line1 = new TLine(lim1[i], 0, lim1[i], 3);
		line2 = new TLine(lim2[i], 0, lim2[i], 3);
		line1->SetLineColor(2);
		line1->SetLineWidth(2);
		line2->SetLineColor(2);
		line2->SetLineWidth(2);

		g->Draw("ALP");
		line1->Draw("SAME");
		line2->Draw("SAME");
		c1->SaveAs(out);
	}
	
*/
}