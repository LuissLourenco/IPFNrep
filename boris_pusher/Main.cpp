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
#include "Vec.h"

#define XCANVAS 1000
#define YCANVAS 1000

using namespace std;

int main(){

	double dt = 0.000001;
	int n_iter = 2000000;

	double m = 1;
	double q = 1;
	double omega = 10;
	double k = 1./omega;

	double E0y = 1;
	double B0y = 1;
	double phase = M_PI/2;

	double E0z = 1;
	double B0z = 1;

	double *t;
	Vec *pos, *vel;
	t = new double[n_iter];
	pos = new Vec[n_iter];
	vel = new Vec[n_iter];

	TF2* Ey = new TF2("Ey", "[0]*cos([1]*x[0]-[2]*x[1]+[3])", -100, 100, -100, 100);
	TF2* By = new TF2("By", "[0]*cos([1]*x[0]-[2]*x[1]+[3])", -100, 100, -100, 100);
	Ey->SetParameters(E0y, omega, k, phase);
	By->SetParameters(B0y, omega, k, 0);

	TF2* Ez = new TF2("Ez", "[0]*cos([1]*x[0]-[2]*x[1]+[3])", -100, 100, -100, 100);
	TF2* Bz = new TF2("Bz", "[0]*cos([1]*x[0]-[2]*x[1]+[3])", -100, 100, -100, 100);
	Ez->SetParameters(E0z, omega, k, 0);
	Bz->SetParameters(B0z, omega, k, phase);

	Vec EVec, BVec, S, U, H, UL;
	double ql = dt * q / (2 * m);

	pos[0] = Vec(3, 0., 0., 0.);
	vel[0] = Vec(3, 0., 0., 0.);
	t[0] = 0;

	for(int i = 1; i < n_iter; i++){
		t[i] = i*dt;
		EVec = Vec(3, 0., Ey->Eval(t[i], pos[i-1][0]), Ez->Eval(t[i], pos[i-1][0]));
		BVec = Vec(3, 0., By->Eval(t[i], pos[i-1][0]), Bz->Eval(t[i], pos[i-1][0]));
		H = BVec * ql;
		S = H * (2 / (1 + H.dot(H)));
		U = vel[i-1] + EVec * ql;
		UL = U + (U + (U % H)) % S;
		vel[i] = UL + EVec * ql;
		pos[i] = pos[i-1] + vel[i] * dt;
		printProgress((double)i/n_iter);
	}
	cout << endl;

	int step = 1000;

	double *x, *y, *z, *vx, *vy, *vz;
	x = new double[n_iter/step];
	y = new double[n_iter/step];
	z = new double[n_iter/step];
	vx = new double[n_iter/step];
	vy = new double[n_iter/step];
	vz = new double[n_iter/step];
	for(int i = 0; i < n_iter/step; i++){
		x[i] = pos[i*step][0];
		y[i] = pos[i*step][1];
		z[i] = pos[i*step][2];
		vx[i] = vel[i*step][0];
		vy[i] = vel[i*step][1];
		vz[i] = vel[i*step][2];
		printProgress((double)i/n_iter*step);
	}
	cout << endl;

	//==========================================================
	
	TApplication *myapp=new TApplication("myapp",0,0);

	TCanvas *c = new TCanvas("c0", "", XCANVAS, YCANVAS);

	TGraph2D* graph = new TGraph2D(n_iter/step, x, y, z);
	graph->SetMarkerStyle(20);
	//graph->SetMarkerSize(1);
	//graph->SetMarkerColor(kRed);
	//graph->SetLineColor(kRed);

	c->cd();
	graph->Draw("P0");

	myapp->Run();


	c->SaveAs("Plot.png");
	

	return 0;
}
