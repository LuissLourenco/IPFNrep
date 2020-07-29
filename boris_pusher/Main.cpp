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

template <typename F>
double Der(F&& f, double x){
	double h = 1E-4;
	return (f(x+h)-f(x))/h; 
};

int main(){

	FILE *fo;
 	fo = fopen("BorisOut.txt","w");

	double dt = 0.004;
	int n_iter = 10000;

	int step = 1; //numero de passos a saltar aquando da impressao

	double c = 1;
	double m = 1;
	double q = 1;
	double omega = 1;
	double kx = 1./omega;
	double ky = 0;
	double kz = 0;

	double a0 = 0.5;
	double delta = 1./sqrt(2);

	double t_now;
	Vec pos, vel;

	double tfwhm = 20;
	auto Poly = [](double x){
		return 10*x*x*x-15*x*x*x*x+6*x*x*x*x*x;
	};
	auto Envelope = [&](double phi){
		if(phi > 2 * tfwhm) return 0.;
		if(phi > tfwhm) return Poly((-phi+2*tfwhm)/tfwhm);
		if(phi > 0) return Poly(phi / tfwhm);
		else return 0.;
 	};
	

	auto Phi = [&](double x, double y, double z, double t){
		return omega*t - kx*x - ky*y - kz*z;
	};

	auto Ax = [&](double phi){
		return 0*Envelope(phi);
	};
	auto Ay = [&](double phi){
		return a0*delta*sin(phi)*Envelope(phi);
	};
	auto Az = [&](double phi){
		return a0*sqrt(1-delta*delta)*cos(phi)*Envelope(phi);
	};


	Vec EVec, BVec, S, U, H, UL;
	double phi, bx_aux, by_aux, bz_aux, gamma, ql;

	t_now = 0;
	pos = Vec(3, 0., 0., 1.);
	phi = Phi(pos[0], pos[1], pos[2], t_now);
	vel = Vec(3, -Ax(phi), -Ay(phi), -Az(phi));

	double *t, *x, *y, *z, *vx, *vy, *vz, *g;
	t = new double[n_iter/step];
	x = new double[n_iter/step];
	y = new double[n_iter/step];
	z = new double[n_iter/step];
	vx = new double[n_iter/step];
	vy = new double[n_iter/step];
	vz = new double[n_iter/step];
	g = new double[n_iter/step];

	fprintf(fo, "%.14e \t ", t_now); 
	fprintf(fo, "%.14e \t ", pos[0]); 
	fprintf(fo, "%.14e \t ", pos[1]); 
	fprintf(fo, "%.14e \t ", pos[2]); 
	fprintf(fo, "%.14e \t ", vel[0]); 
	fprintf(fo, "%.14e \t ", vel[1]); 
	fprintf(fo, "%.14e \t ", vel[2]); 
	fprintf(fo, "%.14e", sqrt(1 + (vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]) / c / c)); 
	
	for(int i = 1; i < n_iter; i++){

		t_now = i*dt;
		phi = Phi(pos[0], pos[1], pos[2], t_now);

		EVec = Vec(3, -Der(Ax, phi)*omega, -Der(Ay, phi)*omega, -Der(Az, phi)*omega);
		bx_aux = Der(Az, phi)*(-ky) - Der(Ay, phi)*(-kz);
		by_aux = Der(Ax, phi)*(-kz) - Der(Az, phi)*(-kx);
		bz_aux = Der(Ax, phi)*(-ky) - Der(Ay, phi)*(-kx);
		BVec = Vec(3, bx_aux, by_aux, bz_aux);

		gamma = sqrt(1 + (vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]) / c / c);
		ql = dt * q / (2 * m * c * gamma);

		H = BVec * ql;
		S = H * (2 / (1 + H.dot(H)));
		U = vel + EVec * ql;
		UL = U + (U + (U % H)) % S;

		gamma = sqrt(1 + (vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]) / c / c);
		ql = dt * q / (2 * m * c * gamma);

		vel = UL + EVec * ql;
		pos = pos + vel * dt;

		if(i % step == 0){
			t [i / step - 1] = t_now;
			x [i / step - 1] = pos[0];
			y [i / step - 1] = pos[1];
			z [i / step - 1] = pos[2];
			vx[i / step - 1] = vel[0];
			vy[i / step - 1] = vel[1];
			vz[i / step - 1] = vel[2];
			g [i / step - 1] = gamma;
			fprintf(fo, "\n%.14e \t ", t [i / step - 1]); 
			fprintf(fo, "%.14e \t ", x [i / step - 1]); 
			fprintf(fo, "%.14e \t ", y [i / step - 1]); 
			fprintf(fo, "%.14e \t ", z [i / step - 1]); 
			fprintf(fo, "%.14e \t ", vx[i / step - 1]); 
			fprintf(fo, "%.14e \t ", vy[i / step - 1]); 
			fprintf(fo, "%.14e \t ", vz[i / step - 1]); 
			fprintf(fo, "%.14e",  g [i / step - 1]); 
		}

		printProgress((double)i/n_iter);

	}

	fclose(fo);


	//==========================================================
	
	TApplication *myapp=new TApplication("myapp",0,0);

	TCanvas *canvas = new TCanvas("c0", "", XCANVAS, YCANVAS);
	
	//TGraph* graph = new TGraph(n_iter/step, x, y);

	TGraph2D* graph = new TGraph2D(n_iter/step, x, y, z);
	graph->SetMarkerStyle(20);
	graph->SetMarkerSize(0.5);
	//graph->SetMarkerColor(kRed);
	//graph->SetLineColor(kRed);
	graph->SetTitle("Boris;x;y;z");

	canvas->cd();
	graph->Draw("P0");

	myapp->Run();


	canvas->SaveAs("Plot.png");
	

	return 0;
}
