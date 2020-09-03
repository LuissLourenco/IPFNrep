#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <chrono>
#include <ctime> 

#include "theory3.c"
#include "DataAnalysis.cpp"
#include "Read.cpp"

#include "TCanvas.h"
#include "TGraph2D.h"
#include "TApplication.h"

using namespace std;

int main(int argc, char **argv){

	clock_t c1 = clock();

	int process;

	double px0 = -2000;
	double kdamp = 0;
	double T = 120;
	int N = 100000;
	int pri = 100;
	int wave_type = 3;
	double tfwhm = 50;
	double stable = 0;
	double Eo = 1;
	double delta = 0;
	double w0 = 5;
	double lambda = 1;

	int l, p;
	if(argc == 6){
		sscanf(argv[1], "%i", &process);
		sscanf(argv[2], "%i", &l);
		sscanf(argv[3], "%i", &p);
		sscanf(argv[4], "%lf", &px0);
		sscanf(argv[5], "%lf", &kdamp);
	}else{
		cout << "ARGUMENT ERROR!" << argc << endl;
		return 1;
	}

	time_t start_time = chrono::system_clock::to_time_t(chrono::system_clock::now());
	cout << "Simulating l = " << l << " and p = " << p << " with px0 = " << px0;
	cout << " starting " << ctime(&start_time);

    int n_cols, n_points;
	double** values;

	double y_min = -20;
	double y_max = 20;
	double dy = 0.5;
	double y = y_min;

	double z_min = y_min;
	double z_max = y_max;
	double dz = dy;
	double z = z_min;

	FILE *input, *output;
	
	
	DataSet X(1);
	DataSet Y(1);
	DataSet Z(1);
	
	
	output = fopen(("Output_wv3_p"+to_string(p)+"l"+to_string(l)+"_px"+to_string(px0)+"_kdamp"+to_string(kdamp)+".txt").c_str(),"w");
	//output = fopen("teste","w");
	fprintf(output, "y\tz\tx_f\ty_f\tz_f\tpx_f\tpy_f\tpz_f\tp_y_max\tp_z_max\tE_f\tL_x_max\tL_x_f\t"); 
	fprintf(output, "px0=%.10e|", px0);	//p01
	fprintf(output, "kdamp=%.10e|", kdamp);	//kdamp
	fprintf(output, "T=%.10e|", T);	//T
	fprintf(output, "N=%i|", N);	//N
	fprintf(output, "pri=%i|", pri);	//pri
	fprintf(output, "tfwhm=%.10e|", tfwhm);	//tfwhm
	fprintf(output, "stable=%.10e|", stable);	//stable
	fprintf(output, "Eo=%.10e|", Eo);	//Eo
	fprintf(output, "delta=%.10e|", delta);	//delta
	fprintf(output, "w0=%.10e|", w0);	//w0
	fprintf(output, "lambda=%.10e|", lambda);	//lambda
	fprintf(output, "l=%i|", l);	//l
	fprintf(output, "p=%i", p);	//p

	int progress = 0;
	while(y <= y_max){
		z = z_min;
		while(z <= z_max){
			input = fopen(("InputToBatch"+to_string(process)+".txt").c_str(),"w");
			fprintf(input, "trash|x01|x02|x03|p01|p02|p03|kdamp|T|N|pri|dx|wave_type|tfwhm|stable|E0|delta|w0|lambda|n|eta|l|p\n"); 
			fprintf(input, "%.10e\n", 0.);	//x01
			fprintf(input, "%.10e\n", y);	//x02
			fprintf(input, "%.10e\n", z);	//x03
			fprintf(input, "%.10e\n", px0);	//p01
			fprintf(input, "%.10e\n", 0.);	//p02
			fprintf(input, "%.10e\n", 0.);	//p03
			fprintf(input, "%.10e\n", kdamp);	//kdamp
			fprintf(input, "%.10e\n", T);	//T
			fprintf(input, "%i\n", N);	//N
			fprintf(input, "%i\n", pri);	//pri
			fprintf(input, "%.10e\n", 5E-3);	//dx
			fprintf(input, "%i\n", wave_type);	//wave_type
			fprintf(input, "%.10e\n", tfwhm);	//tfwhm
			fprintf(input, "%.10e\n", stable);	//stable
			fprintf(input, "%.10e\n", Eo);	//Eo
			fprintf(input, "%.10e\n", delta);	//delta
			fprintf(input, "%.10e\n", w0);	//w0
			fprintf(input, "%.10e\n", lambda);	//lambda
			fprintf(input, "%.10e\n", 1.);	//n
			fprintf(input, "%.10e\n", 1.);	//eta
			fprintf(input, "%i\n", l);	//l
			fprintf(input, "%i", p);	//p
			fclose(input);

			run_theory3(process);

			values = ReadFile(("Out"+to_string(process)+".txt").c_str(), &n_cols, &n_points, false, false);
			// <Out3.txt> t x y z px py pz gamma

			fprintf(output, "\n");
			fprintf(output, "%.10e\t", y);	// y
			fprintf(output, "%.10e\t", z);	// z
			fprintf(output, "%.10e\t", values[1][n_points-1]);	// x_f
			fprintf(output, "%.10e\t", values[2][n_points-1]);	// y_f
			fprintf(output, "%.10e\t", values[3][n_points-1]);	// z_f
			fprintf(output, "%.10e\t", values[4][n_points-1]);	// px_f
			fprintf(output, "%.10e\t", values[5][n_points-1]);	// py_f
			fprintf(output, "%.10e\t", values[6][n_points-1]);	// pz_f

			DataSet PPERP = DataSet(n_points, values[5])*DataSet(n_points, values[5])+DataSet(n_points, values[6])*DataSet(n_points, values[6]);
			int pperpmaxi = PPERP.getMaxI();

			fprintf(output, "%.10e\t", DataSet(n_points, values[5])[pperpmaxi].val());	// p_y_max
			fprintf(output, "%.10e\t", DataSet(n_points, values[6])[pperpmaxi].val());	// p_z_max
			fprintf(output, "%.10e\t", sqrt(1 + values[4][n_points-1]*values[4][n_points-1]
											+ values[5][n_points-1]*values[5][n_points-1]
											+ values[6][n_points-1]*values[6][n_points-1]));	// E_f
			DataSet L_X = abs(DataSet(n_points, values[2]) * DataSet(n_points, values[6]) - DataSet(n_points, values[3]) * DataSet(n_points, values[5]));
			fprintf(output, "%.10e\t", L_X.getMax().val());	// L_x_max
			fprintf(output, "%.10e", L_X[-1].val());	// L_x_f

			
			X = X.concat(DataSet(n_points, values[1]));
			Y = Y.concat(DataSet(n_points, values[2]));
			Z = Z.concat(DataSet(n_points, values[3]));
			
			
			progress++;
			printf("\rSTEP %i of %i  |  %.5lf %%  |  RUN = %.3lf s  |  TIME LEFT = %.2lf min      ", progress, 
								(int)(((y_max-y_min)/dy+1)*((z_max-z_min)/dz+1)), 
								progress/(((y_max-y_min)/dy+1)*((z_max-z_min)/dz+1))*100, 
								(double)(clock()-c1)/(double)CLOCKS_PER_SEC,  
								((((y_max-y_min)/dy+1)*((z_max-z_min)/dz+1))-progress)*(double)(clock()-c1)/(double)CLOCKS_PER_SEC/60);
			fflush(stdout);
			c1 = clock();

			z += dz;
		}
		y += dy;
	}

	fclose(output);

	cout << endl << "Saved file <Output_wv3_p" << p << "l" << l << "_px" << to_string(px0) << "_kdamp" << kdamp << ".txt>" << endl;

	
    TApplication* MyRootApp;
	MyRootApp = new TApplication("MyRootApp", NULL, NULL);

	TCanvas* ca1 = new TCanvas("ca1", "", 1500, 1000);
	TGraph2D* graph = GetTGraph2D(X, Y, Z);
	graph->SetMarkerStyle(8);
	graph->SetMarkerSize(0.1);
	graph->SetMarkerColor(kRed);
	ca1->cd();
	graph->Draw("P");
	//ca1->SaveAs("Laguerre_Trajectories.png");

	MyRootApp->Run();
	

	return 0;

}

