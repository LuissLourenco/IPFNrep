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

	double px0 = 0;
	double kdamp = 0;
	double T = 200;
	int N = 100000;
	int pri = 100;
	double tfwhm = 50;
	double stable = 0;
	double Eo = 1;
	double delta = 0;
	double w0 = 10;
	double lambda = 1;

	int l, p;
	if(argc == 3){
		sscanf(argv[1], "%i", &l);
		sscanf(argv[2], "%i", &p);
	}else{
		l = 0;
		p = 0;
	}

	time_t start_time = chrono::system_clock::to_time_t(chrono::system_clock::now());
	cout << "Simulating l = " << l << " and p = " << p;
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

	//DataSet X(1);
	//DataSet Y(1);
	//DataSet Z(1);

	output = fopen(("Output_p"+to_string(p)+"l"+to_string(l)+"_px"+to_string(px0)+".txt").c_str(),"w");
	fprintf(output, "y\tz\tx_f\ty_f\tz_f\tpx_f\tpy_f\tpz_f\tp_y_max\tE_f\t"); 
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
			input = fopen("InputToBatch.txt","w");
			fprintf(input, "trash|x01|x02|x03|p01|p02|p03|kdamp|T|N|pri|dx|dy|wave_type|tfwhm|stable|E0|delta|w0|lambda|n|eta|l|p\n"); 
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
			fprintf(input, "%.10e\n", 5E-3);	//dy
			fprintf(input, "%i\n", 2);	//wave_type
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

			run_theory3();

			values = ReadFile("Out3.txt", &n_cols, &n_points, false, false);
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
			fprintf(output, "%.10e\t", (abs(DataSet(n_points, values[5]))).getMax().val());	// p_perp_max
			fprintf(output, "%.10e", sqrt(1 + values[4][n_points-1]*values[4][n_points-1]
											+ values[5][n_points-1]*values[5][n_points-1]
											+ values[6][n_points-1]*values[6][n_points-1]));	// E_f

			//X.append(Var(y));
			//Y.append(Var(z));
			//Z.append(Var(values[2][n_points-1]-y));
			
			progress++;

			printf("\rSTEP %i of %i \t  %.5lf %%  \t RUN = %.3lf s \t TIME LEFT = %.2lf min      ", progress, 
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

	cout << endl << "Saved file <Output_p" << p << "l" << l << "_px" << to_string(px0) << ".txt>" << endl;

	/*
    TApplication* MyRootApp;
	MyRootApp = new TApplication("MyRootApp", NULL, NULL);

	TCanvas* c1 = new TCanvas("c1", "", 1500, 1000);
	TGraph2D* graph = GetTGraph2D(X, Y, Z);
	c1->cd();
	graph->Draw("TRI");
	c1->SaveAs("Laguerre_Trajectories.png");

	MyRootApp->Run();
	*/
	return 0;

}

