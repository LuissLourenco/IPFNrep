#include <iostream>
#include <cstdio>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <ctime>

#include "ODEpoint.h"
#include "ODEsolver.h"
#include "DataAnalysis.h"

#define XCANVAS 900*2
#define YCANVAS 1950*2

using namespace std;

int main(){

	double E0 = 1;
	double B0 = E0;
	double omega = 2*M_PI/1;
	double k = 1./omega;

	double t_max = 100;
	double dt = 1E-3;

	vector<double> y0;
	y0.push_back(0);
	y0.push_back(0);
	y0.push_back(0);
	y0.push_back(0);

	ODEpoint p0(0, y0);

	TFormula dxdt("dxdt", "x[2]");
	TFormula dydt("dydt", "x[3]");
	TFormula dvxdt("dvxdt", "-x[3]*[0]*cos([1]*x[4] - [2]*x[0])");
	TFormula dvydt("dvydt", "(-[3] + x[2]*[0])*cos([1]*x[4] - [2]*x[0])");
	dvxdt.SetParameters(B0, omega, k);
	dvydt.SetParameters(B0, omega, k, E0);

	vector<TFormula> form;
	form.push_back(dxdt);
	form.push_back(dydt);
	form.push_back(dvxdt);
	form.push_back(dvydt);
	
	ODEsolver solver(form);
	vector<ODEpoint> Solution = solver.RK4solver(p0, 0, t_max, dt);

	double *t_arr = new double[Solution.size()];
	double *x_arr = new double[Solution.size()];
	double *y_arr = new double[Solution.size()];
	double *vx_arr = new double[Solution.size()];
	double *vy_arr = new double[Solution.size()];

	for(int i=0; i<Solution.size(); i++){
		t_arr[i] = Solution[i].Get_Time();
		x_arr[i] = Solution[i].Get_Var_vec()[0];
		y_arr[i] = Solution[i].Get_Var_vec()[1];
		vx_arr[i] = Solution[i].Get_Var_vec()[2];
		vy_arr[i] = Solution[i].Get_Var_vec()[3];
	}

	DataSet VX(Solution.size(), vx_arr, 0);
	double v_drift = 0;//VX.getMean().val();

	for(int i = 0; i < VX.size(); i++) v_drift += VX[i].val() * dt / t_max;

	cout << "V_DRIFT = " << v_drift << endl;

	for(int i=0; i<Solution.size(); i++){
		x_arr[i] -= v_drift*t_arr[i];
	}

	//==========================================================

	TCanvas *c = new TCanvas("c0", "", XCANVAS, YCANVAS);

	TGraph* graph = new TGraph(Solution.size(), x_arr, y_arr);
	graph->SetMarkerSize(0.1);
	graph->SetMarkerColor(kRed);
	graph->SetLineColor(kRed);

	c->cd();
	graph->Draw("AC");
	c->SaveAs("Orbit.png");
	
	return 0;
}
