#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1

#include <math.h>
#include <iostream>

#include "TVectorField.cpp"
#include "DataAnalysis.cpp"

using namespace std;

int main(){

	int n_points, n_cols;
	double** values = ReadFile("Output_wv3_p3l3_px-0.000000.txt", &n_cols, &n_points, true);
	//y	z	x_f	y_f	z_f	px_f	py_f	pz_f	p_y_max	E_f	
	DataSet Y(n_points, values[0]);
	DataSet Z(n_points, values[1]);
	DataSet YF(n_points, values[3]);
	DataSet ZF(n_points, values[4]);

	TCanvas* c1 = new TCanvas("c", "", 1000, 1000);

	c1->cd();

	TVectorField* v = new TVectorField(Y, Z, YF, ZF);

	//v->SetLimits(-10, -10, 10, 10);

	v->Draw();


	
	c1->SaveAs("test.png");

	return 0;

}
