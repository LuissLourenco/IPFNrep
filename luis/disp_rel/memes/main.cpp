#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>

#include <TF2.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TStyle.h>

using namespace std;

double n=1;
double eta=1;
double w0=10;
double lambda=1;
double Eo=1;
double zr = M_PI*w0*w0*n/lambda;
double p=1,l=1;

double w=1;

double* k(x,y,z){
	double r = sqrt(x*x+y*y);
	double phi = atan2(y,x);

	double res[3];
	res[0] -
}

double Ax(double x, double y, double z, double t){
	double r = sqrt(x*x+y*y);
	double phi = atan2(y,x);
	double wz = w0*sqrt(1+z*z/zr/zr);
	double res = 1;
	res *= Eo*w0/wz * pow(r*sqrt(2.)/wz, l) * assoc_laguerre(abs(p),abs(l), 2.*r*r/wz/wz) * exp(-r*r/wz/wz);
	res *= cos(w*t-k*z-k*r*r*z/2./(z*z+zr*zr)-l*phi+(2.*p+abs(l)+1.)*atan(z/zr));
	return res;
}

double Ex(double x, double y, double z, double t){
	double h=1e-5;
	return -(Ax(x,y,z,t+h/2.)-Ax(x,y,z,t-h/2.))/h;
}
double By(double x, double y, double z, double t){
	double h=1e-5;
	return (Ax(x,y,z+h/2.,t)-Ax(x,y,z-h/2.,t))/h;
}
double Bz(double x, double y, double z, double t){
	double h=1e-5;
	return -(Ax(x,y+h/2.,z,t)-Ax(x,y-h/2.,z,t))/h;
}



int main(int argc, char** argv){



	return 0;
}