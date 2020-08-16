
//#include "laguerre_gaussian.cpp"

#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>

#include <TF2.h>
#include <TCanvas.h>

using namespace std;



complex<double> imu(0,1);

double n=1;
double eta=1;
double w0=10;
double lambda=1;
double w =1;
double Eo=1;


complex<double> u_aux(double x, double y, double z, double l, double p){

	complex<double> res(1,0);

	double r = sqrt(y*y+z*z);
	double phi = atan2(z,y);
	double kg = 2*M_PI*n/lambda;
	double zr = M_PI*w0*w0*n/lambda;
	double wz = w0*sqrt(1+x*x/zr/zr);

	res *= Eo*w0/wz;
	res *= pow(r*sqrt(2)/wz, abs(l));
	res *= assoc_laguerre(abs(p), abs(l), 2*r*r/wz/wz);
	res *= exp(-r*r/wz/wz);
	res *= exp(-imu*kg*r*r*x/2./(x*x+zr*zr));
	res *= exp(-imu*l*phi);
	res *= exp(imu*(2*p+abs(l)+1)*atan(x/zr));

	return res;
}

double nao_delta(int a, int b){
	if(a!=b) return 1;
	else return 0;
}


complex<double> magnetic_field_luis(double x, double y, double z, double l, double p, double t){
	//\vec{B} = B \hat{z}
	complex<double> res;

	double r = sqrt(y*y+z*z);
	double phi = atan2(z,y);
	double kg = 2*M_PI*n/lambda;
	double zr = M_PI*w0*w0*n/lambda;
	double wz = w0*sqrt(1+x*x/zr/zr);

	complex<double> mood, mood1, mood2;

	mood1 = -x/(x*x+zr*zr)-abs(l)*x/(x*x+zr*zr)*nao_delta(l,0)+2*r*r*x*zr*zr/w0/w0/(x*x+zr*zr)/(x*x+zr*zr)-imu*kg*r*r*(0.5-x*x/(x*x+zr*zr))/(x*x+zr*zr)+imu*(2.*p+abs(l)+1.)*zr/(x*x+zr*zr)-imu*kg;
	mood1 *= u_aux(x,y,z,p,l);
	mood2 = u_aux(x,y,z,p-1,l+1)*4.*r*r*x*zr*zr/w0/w0/(x*x+zr*zr)/(x*x+zr*zr)*nao_delta(p,0);
	mood = mood1+mood2;

	res=mood*exp(imu*(w*t-kg*x));

	return res;
}

double fun1(double*x,double*par){
	return magnetic_field_luis(par[0],x[0],x[1],par[1],par[2],par[3]).real();
}


int main(){

	auto f1 = new TF2("f1",fun1,-50,50,-50,50,4);
	f1->SetParameter(0,0); //x
	f1->SetParameter(1,100); //l
	f1->SetParameter(2,100); //p
	f1->SetParameter(3,0); //t
	f1->SetNpx(1000);
	f1->SetNpy(1000);

	auto c1 = new TCanvas("c1", "", 1200, 1200);
	f1->Draw("colz");

	c1->SaveAs("Plot.png");



	delete c1;
	delete f1;


	return 0;
}
