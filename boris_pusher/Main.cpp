#include <iostream>
#include <cstdio>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <ctime>

#include "DataAnalysis.h"
#include "Vec.h"

using namespace std;

template <typename F>
double Der(F&& f, double x){
	double h = 1E-5;
	return (f(x+h/2.)-f(x-h/2.))/h; 
};

int main(){

	FILE *fo;
 	fo = fopen("BorisOut.txt","w");

	double dt = 1e-3;
	double t_max = 150;
	int n_iter = t_max / dt;

	int step = 1; //numero de passos a saltar aquando da impressao
	double t_stabilize = 0; //tempo a partir do qual comeca a guardar dados

	double c = 1;
	double m = 1;
	double q = 1;
	double omega = 1;
	double kx = 1./omega;
	double ky = 0;
	double kz = 0;

	double a0 = 1;
	double delta = 1./sqrt(2);

	double t_now;
	Vec pos, vel;

	// Funcoes para criar o Envelope
	double rise = 50;
	double stable = 30;
	auto Poly = [](double x){
		return 10*x*x*x-15*x*x*x*x+6*x*x*x*x*x;
	};
	auto Envelope = [&](double phi){
		if(phi > 2 * rise + stable) return 0.;
		if(phi > rise + stable) return Poly((-phi+2*rise+stable)/rise);
		if(phi > rise) return 1.;
		if(phi > 0) return Poly(phi / rise);
		else return 0.;
 	};
	
 	// Funcoes da Onda Eletromagnetica
	auto Phi = [&](double x, double y, double z, double t){
		return omega*t - kx*x - ky*y - kz*z;
	};

	auto Ax = [&](double phi){
		return 0*Envelope(phi);
	};
	auto Ay = [&](double phi){
		return a0*delta*cos(phi)*Envelope(phi);
	};
	auto Az = [&](double phi){
		return a0*sqrt(1-delta*delta)*sin(phi)*Envelope(phi);
	};


	Vec EVec, BVec, S, U, H, UL;
	double phi, bx_aux, by_aux, bz_aux, gamma, ql;

	t_now = 0;
	pos = Vec(3, 0., 0., 1.);
	phi = Phi(pos[0], pos[1], pos[2], t_now);
	//vel = Vec(3, 1.016*a0*a0/4, 0., 0.); // Average Rest Frame (?)
	vel = Vec(3, 0., 0., 0.);

	if(t_now >= t_stabilize){
		fprintf(fo, "%.14e \t ", t_now); 
		fprintf(fo, "%.14e \t ", pos[0]); 
		fprintf(fo, "%.14e \t ", pos[1]); 
		fprintf(fo, "%.14e \t ", pos[2]); 
		fprintf(fo, "%.14e \t ", vel[0]); 
		fprintf(fo, "%.14e \t ", vel[1]); 
		fprintf(fo, "%.14e \t ", vel[2]); 
		fprintf(fo, "%.14e\n", sqrt(1 + (vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]) / c / c)); 
	}

	for(int i = 1; i < n_iter; i++){

		t_now = i*dt;
		phi = Phi(pos[0], pos[1], pos[2], t_now);

		// Calculo dos campos
		EVec = Vec(3, -Der(Ax, phi)*omega, -Der(Ay, phi)*omega, -Der(Az, phi)*omega);
		bx_aux = Der(Az, phi)*(-ky) - Der(Ay, phi)*(-kz);
		by_aux = Der(Ax, phi)*(-kz) - Der(Az, phi)*(-kx);
		bz_aux = Der(Ax, phi)*(-ky) - Der(Ay, phi)*(-kx);
		BVec = Vec(3, bx_aux, by_aux, bz_aux);

		// Constantes
		gamma = sqrt(1 + (vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]) / c / c);
		ql = dt * q / (2 * m * c * gamma);

		// Avancar metade do efeito de E e rodar por B
		H = BVec * ql;
		S = H * (2. / (1. + H.dot(H)));
		U = vel + EVec * ql;
		UL = U + (U + (U % H)) % S;

		// Atualizar as constantes
		gamma = sqrt(1 + (vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]) / c / c);
		ql = dt * q / (2. * m * c * gamma);

		// Aplicar a outra metade do efeito de E e atualizar posicoes
		vel = UL + EVec * ql;
		pos = pos + vel * dt;

		// Print e guardar tudo
		if(i % step == 0 && t_now >= t_stabilize){
			fprintf(fo, "%.14e \t ", t_now); 
			fprintf(fo, "%.14e \t ", pos[0]); 
			fprintf(fo, "%.14e \t ", pos[1]); 
			fprintf(fo, "%.14e \t ", pos[2]); 
			fprintf(fo, "%.14e \t ", vel[0]); 
			fprintf(fo, "%.14e \t ", vel[1]); 
			fprintf(fo, "%.14e \t ", vel[2]); 
			fprintf(fo, "%.14e\n",  gamma); 
		}

		printProgress((double)i/n_iter);

	}
	cout << "Закончилось" << endl;

	fclose(fo);

	return 0;
}
