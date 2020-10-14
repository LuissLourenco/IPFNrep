#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <chrono>
#include <ctime> 

#include "theory4.c"
#include "DataAnalysis.cpp"

#include "TCanvas.h"
#include "TGraph2D.h"
#include "TApplication.h"
#include "TPolyLine3D.h"

using namespace std;



void input(double x01, double x02, double x03, double p01, double p02, double p03, double kdamp, double T, int N, int pri, double dx, int wave_type, double tfwhm, double stable, double Eo, double delta, double w0, double lambda, double n, double eta, int l, int p){
	FILE* itb = fopen("../../rk4_cooling/InputToBatch.txt","w");
	fprintf(itb, "trash|x01|x02|x03|p01|p02|p03|kdamp|T|N|pri|dx|wave_type|tfwhm|stable|E0|delta|w0|lambda|n|eta|l|p\n"); 
	fprintf(itb, "%.10e\n", x01);	//x01
	fprintf(itb, "%.10e\n", x02);	//x02
	fprintf(itb, "%.10e\n", x03);	//x03
	fprintf(itb, "%.10e\n", p01);	//p01
	fprintf(itb, "%.10e\n", p02);	//p02
	fprintf(itb, "%.10e\n", p03);	//p03
	fprintf(itb, "%.10e\n", kdamp);	//kdamp
	fprintf(itb, "%.10e\n", T);	//T
	fprintf(itb, "%i\n", N);	//N
	fprintf(itb, "%i\n", pri);	//pri
	fprintf(itb, "%.10e\n", dx);	//dx
	fprintf(itb, "%i\n", wave_type);	//wave_type
	fprintf(itb, "%.10e\n", tfwhm);	//tfwhm
	fprintf(itb, "%.10e\n", stable);	//stable
	fprintf(itb, "%.10e\n", Eo);	//Eo
	fprintf(itb, "%.10e\n", delta);	//delta
	fprintf(itb, "%.10e\n", w0);	//w0
	fprintf(itb, "%.10e\n", lambda);	//lambda
	fprintf(itb, "%.10e\n", n);	//n
	fprintf(itb, "%.10e\n", eta);	//eta
	fprintf(itb, "%i\n", l);	//l
	fprintf(itb, "%i", p);	//p
	fclose(itb);
}


int main(int argc, char **argv){

	
	
	double a0=10;

	int p=0;
	int l=1;

	double r0=2;
	double phi0=180*M_PI/180;

	double p0=-10;

	double dt=2e-3;
	double T=5e3;
	int N=T/dt;

	

	int num=7;
	for(int i=num; i<=num; i++){

		input(0,r0*cos(phi0),r0*sin(phi0),  p0,0,0,   1.18e-8, T, N, 10, 5e-3, 3, 50, T*10, a0, -1, 5, 1, 1, 1, l, p);
		
		run_theory3(i);
	}
	cout<<endl<<"gotovo"<<endl;





	

	return 0;
}

