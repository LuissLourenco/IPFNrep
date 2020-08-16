#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>

#include "theory3.c"
#include "DataAnalysis.cpp"
#include "Read.cpp"

#include "TCanvas.h"
#include "TGraph2D.h"
#include "TApplication.h"

using namespace std;

int main(){

	int n_cols, n_points;
	double** values;

	double r_min = 0;
	double r_max = 20;
	double dr = 1;
	double r = r_min;

	double phi_min = 0;
	double phi_max = 2*M_PI;
	double dphi = M_PI/16;
	double phi = phi_min;

	FILE *file;

	DataSet X(1);
	DataSet Y(1);
	DataSet Z(1);

	while(phi <= phi_max){
		r = r_min;
		while(r <= r_max){
			file = fopen("InputToBatch.txt","w");
			fprintf(file, "trash|x01|x02|x03|p01|p02|p03|kdamp|T|N|pri|dx|dy|wave_type|tfwhm|stable|E0|delta|w0|lambda|n|eta|l|p\n"); 
			fprintf(file, "%.10e\n", 0.);	//x01
			fprintf(file, "%.10e\n", r*cos(phi));	//x02
			fprintf(file, "%.10e\n", r*sin(phi));	//x03
			fprintf(file, "%.10e\n", 0.);	//p01
			fprintf(file, "%.10e\n", 0.);	//p02
			fprintf(file, "%.10e\n", 0.);	//p03
			fprintf(file, "%.10e\n", 0.);	//kdamp
			fprintf(file, "%.10e\n", 100.);	//T
			fprintf(file, "%i\n", 10000);	//N
			fprintf(file, "%i\n", 10);	//pri
			fprintf(file, "%.10e\n", 5E-3);	//dx
			fprintf(file, "%.10e\n", 5E-3);	//dy
			fprintf(file, "%i\n", 2);	//wave_t
			fprintf(file, "%.10e\n", 50.);	//tfwhm
			fprintf(file, "%.10e\n", 0.);	//stable
			fprintf(file, "%.10e\n", 1.);	//Eo
			fprintf(file, "%.10e\n", 0.);	//delta
			fprintf(file, "%.10e\n", 10.);	//w0
			fprintf(file, "%.10e\n", 1.);	//lambda
			fprintf(file, "%.10e\n", 1.);	//n
			fprintf(file, "%.10e\n", 1.);	//eta
			fprintf(file, "%i\n", 4);	//l
			fprintf(file, "%i\n", 4);	//p
			fclose(file);

			run_theory3();

			values = ReadFile("Out3.txt", &n_cols, &n_points, false, false);

			X = X.concat(DataSet(n_points, values[1], 0));
			Y = Y.concat(DataSet(n_points, values[2], 0));
			Z = Z.concat(DataSet(n_points, values[3], 0));

			cout << "phi = " << phi << "\t r = " << r << endl;

			r += dr;
		}
		phi += dphi;
	}

    TApplication* MyRootApp;
	MyRootApp = new TApplication("MyRootApp", NULL, NULL);

	TCanvas* c1 = new TCanvas("c1", "", 1500, 1000);
	TGraph2D* graph = GetTGraph2D(X, Y, Z);
	c1->cd();
	graph->Draw("P");
	c1->SaveAs("Laguerre_Trajectories.png");

	MyRootApp->Run();

	return 0;

}
