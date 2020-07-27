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

	double dt = 0.00001;
	int n_iter = 2000000;

	double c = 1;
	double m = 1;
	double q = 1;
	double omega = 1;
	double kx = 1./omega;
	double ky = 0;
	double kz = 0;

	double a0 = 0.5;
	double delta = 1;//1./sqrt(2);

	double *t;
	Vec *pos, *vel;
	t = new double[n_iter];
	pos = new Vec[n_iter];
	vel = new Vec[n_iter];

	TF2* Phi = new TF2("phi", "[omega]*x[0]-[kx]*x[1]-[ky]*x[2]-[kz]*x[3]", -1000, 1000, -1000, 1000);
	Phi->SetParameter("omega", omega); 
	Phi->SetParameter("kx", kx);
	Phi->SetParameter("ky", ky);
	Phi->SetParameter("kz", kz);

	TF1* Ax = new TF1("Ax", "0", -1000, 1000);
	TF1* Ay = new TF1("Ay", "[a0]*[delta]*sin(x)", -1000, 1000);
	TF1* Az = new TF1("Az", "[a0]*sqrt(1-[delta]*[delta])*cos(x)", -1000, 1000);
	Ay->SetParameter("a0", a0);
	Ay->SetParameter("delta", delta);
	Az->SetParameter("a0", a0);
	Az->SetParameter("delta", delta);

	Vec EVec, BVec, S, U, H, UL;
	double phi, bx_aux, by_aux, bz_aux, gamma, ql;

	t[0] = 0;
	pos[0] = Vec(3, 0., 0., 1.);
	phi = Phi->Eval(t[0], pos[0][0], pos[0][1], pos[0][2]);
	vel[0] = Vec(3, -a0*a0/(4+a0*a0), Ay->Eval(phi), Az->Eval(phi));
	
	for(int i = 1; i < n_iter; i++){

		gamma = sqrt(1 + (vel[i-1][0] * vel[i-1][0] + vel[i-1][1] * vel[i-1][1] + vel[i-1][2] * vel[i-1][2]) / c / c);
		ql = dt * q / (2 * m * c * gamma);

		t[i] = i*dt;
		phi = Phi->Eval(t[i], pos[i-1][0], pos[i-1][1], pos[i-1][2]);

		EVec = Vec(3, -Ax->Derivative(phi)*omega, -Ay->Derivative(phi)*omega, -Az->Derivative(phi)*omega);
		bx_aux = Az->Derivative(phi)*(-ky) - Ay->Derivative(phi)*(-kz);
		by_aux = Ax->Derivative(phi)*(-kz) - Az->Derivative(phi)*(-kx);
		bz_aux = Ax->Derivative(phi)*(-ky) - Ay->Derivative(phi)*(-kx);
		BVec = Vec(3, bx_aux, by_aux, bz_aux);

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

	TCanvas *canvas = new TCanvas("c0", "", XCANVAS, YCANVAS);
	
	//TGraph* graph = new TGraph(n_iter/step, x, y);

	TGraph2D* graph = new TGraph2D(n_iter/step, x, y, z);
	graph->SetMarkerStyle(20);
	//graph->SetMarkerSize(1);
	//graph->SetMarkerColor(kRed);
	//graph->SetLineColor(kRed);
	graph->SetTitle("Boris;x;y;z");

	canvas->cd();
	graph->Draw("P0");

	myapp->Run();


	canvas->SaveAs("Plot.png");
	

	return 0;
}
