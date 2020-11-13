#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <chrono>
#include <ctime> 

#include "../src/theory3.c"

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
	double Eo = 30;
	double delta = 0;
	double w0 = 5;
	double lambda = 1;
	int l = 0;
	int p = 0;

	double r, phi, y, z;


	if(argc == 18){
		sscanf(argv[1], "%i", &process);
		sscanf(argv[2], "%lf", &r);
		sscanf(argv[3], "%lf", &phi);
		sscanf(argv[4], "%lf", &px0);
		sscanf(argv[5], "%lf", &kdamp);
		sscanf(argv[6], "%lf", &T);
		sscanf(argv[7], "%i", &N);
		sscanf(argv[8], "%i", &pri);
		sscanf(argv[9], "%i", &wave_type);
		sscanf(argv[10], "%lf", &tfwhm);
		sscanf(argv[11], "%lf", &stable);
		sscanf(argv[12], "%lf", &Eo);
		sscanf(argv[13], "%lf", &delta);
		sscanf(argv[14], "%lf", &w0);
		sscanf(argv[15], "%lf", &lambda);
		sscanf(argv[16], "%i", &l);
		sscanf(argv[17], "%i", &p);
	}else{
		cout << "ARGUMENT ERROR!" << argc << endl;
		return 1;
	}

	y = r * cos(phi/180*M_PI);
	z = r * sin(phi/180*M_PI);

	FILE *input;

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

	return 0;

}

