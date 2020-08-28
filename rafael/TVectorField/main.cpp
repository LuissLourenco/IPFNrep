#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1

#include <math.h>
#include <iostream>

#include "TVectorField.cpp"
#include "DataAnalysis.cpp"

using namespace std;

int main(){

	int n_points, n_cols;
	double** values = ReadFile("Output_wv3_p2l3_px-0.000000.txt", &n_cols, &n_points, true);
	//y	z	x_f	y_f	z_f	px_f	py_f	pz_f	p_y_max p_z_max	E_f	L_x_max L_x_f
	DataSet Y(n_points, values[0]);
	DataSet Z(n_points, values[1]);
	
	DataSet DY = DataSet(n_points, values[10]);
	DataSet DZ = DataSet(n_points, values[10]);

	TCanvas* c1 = new TCanvas("c", "", 2000, 2000);

	c1->cd();

	TVectorField* v = new TVectorField(Y, Z, DY, DZ);

	//v->SetLimits(-10, -10, 10, 10);

	v->Draw("F");


	
	c1->SaveAs("test.png");

	return 0;

}
