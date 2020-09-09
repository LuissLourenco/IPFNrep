
#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <unistd.h>

#include "ram.h"

#include <TF2.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TApplication.h>
#include <TPad.h>
#include <TPadPainter.h>
#include <TImage.h>
#include <TROOT.h>
#include <TPaletteAxis.h>
#include <TH1.h>

using namespace std;

complex<double> imu(0,1);

char trash[128];
double x01,x02,x03,p01,p02,p03, kdamp, T;
long long int N; int pri;
double dx; int wave_type;
double tfwhm, stable, Eo, delta, w0, lambda, n, eta;
int lr,pr;
double k,kg,zr,w,Bo;

double der_coef4[5] = {1./12., -2./3., 0., 2./3., -1./12.};

double p,l,t;





//FIELDS FROM RK4--------------------------------------------------------------

double Envelope(double x, double t){return 1;}

double Efx(double x, double y, double z){  //Ex interpolation to (x,y,z)
 if(wave_type == 0) return 0; 
 if(wave_type == 1) return 0;
 if(wave_type == 2) return 0;
 if(wave_type == 3){ //return 0;
 	double r = sqrt(y*y+z*z);
	double phi = atan2(z,y);
 	double arg = w*t-k*x-l*phi;
 	double amp = Eo * exp(-r*r/w0/w0); amp *= (l!=0) ? pow(r*sqrt(2.)/w0 , abs(l))*assoc_laguerre(p, abs(l), 2.*r*r/w0/w0) : (p!=0) ? assoc_laguerre(p, abs(l), 2.*r*r/w0/w0) : 1;
 	double res = 2.*r/w0/w0*amp;
 	if(l!=0) res-= abs(l)*amp/r;
 	if(p!=0){double amp2 = Eo*exp(-r*r/w0/w0)*assoc_laguerre(p-1, abs(l)+1, 2.*r*r/w0/w0); amp2 *= (l!=0) ? pow(r*sqrt(2.)/w0 , abs(l)) : 1; res+=amp2*4*r/w0/w0;}
 	res *= cos(arg)*cos(phi);
 	if(l!=0) res+=l*amp/r*sin(arg)*sin(phi);
 	return res*Envelope(x,t);
 }
 else return 0.;
}

double DerEfx(double x, double y, double z){ //Ex time derivative at (x,y,z)
 if(wave_type == 0) return 0; 
 if(wave_type == 1) return 0;
 if(wave_type == 2) return 0;
 if(wave_type == 3){ 
 	double r = sqrt(y*y+z*z);
	double phi = atan2(z,y);
 	double arg = w*t-k*x-l*phi;
 	double amp = Eo * exp(-r*r/w0/w0); amp *= (l!=0) ? pow(r*sqrt(2.)/w0 , abs(l))*assoc_laguerre(p, abs(l), 2.*r*r/w0/w0) : (p!=0) ? assoc_laguerre(p, abs(l), 2.*r*r/w0/w0) : 1;
 	double res = 2.*r/w0/w0*amp;
 	if(l!=0) res-= abs(l)*amp/r;
 	if(p!=0){double amp2 = Eo*exp(-r*r/w0/w0)*assoc_laguerre(p-1, abs(l)+1, 2.*r*r/w0/w0); amp2 *= (l!=0) ? pow(r*sqrt(2.)/w0 , abs(l)) : 1; res+=amp2*4*r/w0/w0;}
 	res *= -w*sin(arg)*cos(phi);
 	if(l!=0) res+=w*l*amp/r*cos(arg)*sin(phi);
 	return res*Envelope(x,t);
 }
 else return 0.;
}

double Efy(double x, double y, double z){  //Ey interpolation to (x,y,z)
 if(wave_type == 0) return delta*w*Eo*sin(w*t-k*x)*Envelope(x,t);
 if(wave_type == 1){
 	double r = sqrt(y*y+z*z);
 	double wz = w0 * sqrt(1 + x*x/(zr*zr));
 	double R_1 = x / (x*x + zr*zr);
 	double psi = atan(x/zr); 
 	return w *Eo*w0/wz * exp(-r*r/(wz*wz)) * sin(w*t - kg*x - kg*r*r*R_1/2. + psi) * Envelope(x, t);
 }
 if(wave_type == 2){
 	double r = sqrt(y*y+z*z);
	double phi = atan2(z,y);
	double wz = w0*sqrt(1+x*x/zr/zr);
	double R_1 = x / (x*x + zr*zr);
	double N = abs(l) + 2*p;
	double arg = w*t-k*x-k*r*r*R_1/2.-l*phi+(N+1)*atan(x/zr);
 	double amp = Eo*w0*wz*exp(-r*r/wz/wz); amp *= (l!=0)? pow(r*sqrt(2.)/wz, abs(l))*assoc_laguerre(abs(p), abs(l), 2.*r*r/wz/wz) : (p!=0)? assoc_laguerre(abs(p), abs(l), 2.*r*r/wz/wz) : 1;
 	return w*amp*sin(arg)*Envelope(x,t);
 }
 if(wave_type==3){ //return 0;
 	double r = sqrt(y*y+z*z);
	double phi = atan2(z,y);
	double arg = w*t-k*x-l*phi;
	double amp=1;
	amp *= Eo;
	amp *= exp(-r*r/w0/w0);
	if(l!=0) amp *= pow(r*sqrt(2.)/w0 , abs(l));
	if (l!=0 || p!=0) amp *= assoc_laguerre(abs(p), abs(l), 2.*r*r/w0/w0); 	
 	double res = w*amp*sin(arg); 
 	return res*Envelope(x,t);
 } 
 else return 0.;
}

double DerEfy(double x, double y, double z){ //Ey time derivative at (x,y,z)
 if(wave_type == 0) return delta*w*w*Eo*cos(w*t-k*x)*Envelope(x,t); 
 if(wave_type == 1){
 	double r = sqrt(y*y+z*z);
 	double wz = w0 * sqrt(1 + x*x/(zr*zr));
 	double R_1 = x / (x*x + zr*zr);
 	double psi = atan(x/zr);
 	return w*w * Eo * w0/wz * exp(-r*r/(wz*wz)) * cos(w*t - kg*x - kg*r*r*R_1/2. + psi) * Envelope(x, t);
 }
 if(wave_type == 2){
 	double r = sqrt(y*y+z*z);
	double phi = atan2(z,y);
	double wz = w0*sqrt(1+x*x/zr/zr);
	double R_1 = x / (x*x + zr*zr);
	double N = abs(l) + 2*p;
	double arg = w*t-k*x-k*r*r*R_1/2.-l*phi+(N+1)*atan(x/zr);
 	double amp = Eo*w0*wz*exp(-r*r/wz/wz); amp *= (l!=0)? pow(r*sqrt(2.)/wz, abs(l))*assoc_laguerre(abs(p), abs(l), 2.*r*r/wz/wz) : (p!=0)? assoc_laguerre(abs(p), abs(l), 2.*r*r/wz/wz) : 1;
 	return w*w*amp*cos(arg)*Envelope(x,t);
 }
 if(wave_type==3){ 
 	double r = sqrt(y*y+z*z);
	double phi = atan2(z,y);
	double arg = w*t-k*x-l*phi;
	double amp=1;
	amp *= Eo;
	amp *= exp(-r*r/w0/w0);
	if(l!=0) amp *= pow(r*sqrt(2.)/w0 , abs(l));
	if (l!=0 || p!=0) amp *= assoc_laguerre(abs(p), abs(l), 2.*r*r/w0/w0); 	
 	double res = w*w*amp*cos(arg); 
 	return res*Envelope(x,t);
 } 
 else return 0.;
}

double Efz(double x, double y, double z){  //Ez interpolation to (x,y,z)
 if(wave_type == 0) return -sqrt(1-delta*delta)*w*Eo*cos(w*t-k*x)*Envelope(x,t); 
 if(wave_type == 1) return 0;
 if(wave_type == 2) return 0;
 if(wave_type == 3) return 0;
 else return 0.;
}

double DerEfz(double x, double y, double z){ //Ez time derivative at (x,y,z)
 if(wave_type == 0) return sqrt(1-delta*delta)*w*w*Eo*sin(w*t-k*x)*Envelope(x,t);
 if(wave_type == 1) return 0;
 if(wave_type == 2) return 0;
 if(wave_type == 3) return 0;
 else return 0.;
}

double Bfx( double x, double y, double z){  //Bx interpolation to (x,y,z)
 if(wave_type == 0) return 0; 
 if(wave_type == 1) return 0;
 if(wave_type == 2){
 	double r = sqrt(y*y+z*z);
	double phi = atan2(z,y);
	double wz = w0*sqrt(1+x*x/zr/zr);
	double R_1 = x / (x*x + zr*zr);
	double N = abs(l) + 2*p;
	double arg = w*t-k*x-k*r*r*R_1/2.-l*phi+(N+1)*atan(x/zr);
 	double amp = Eo*w0*wz*exp(-r*r/wz/wz); amp *= (l!=0)? pow(r*sqrt(2.)/wz, abs(l))*assoc_laguerre(abs(p), abs(l), 2.*r*r/wz/wz) : (p!=0)? assoc_laguerre(abs(p), abs(l), 2.*r*r/wz/wz) : 1;
 	double f1 = (2*r/wz/wz-abs(l)/r)*amp;
 	if(p!=0){double amp2 = Eo*w0/wz*exp(-r*r/wz/wz)*assoc_laguerre(p-1, abs(l)+1, 2.*r*r/wz/wz); amp2 *= (l!=0) ? pow(r*sqrt(2.)/wz , abs(l)) : 1; f1+=amp2*4*r/wz/wz;}
 	f1*=sin(phi);
 	double f2 = -(l/r*cos(phi)+k*r*R_1*sin(phi))*amp;
 	return (f1*cos(arg)+f2*sin(arg))*Envelope(x,t);
 }
 if(wave_type == 3){ //return 0;
 	double r = sqrt(y*y+z*z);
	double phi = atan2(z,y);
 	double arg = w*t-k*x-l*phi;
 	double amp = Eo * exp(-r*r/w0/w0); amp *= (l!=0) ? pow(r*sqrt(2.)/w0 , abs(l))*assoc_laguerre(p, abs(l), 2.*r*r/w0/w0) : (p!=0) ? assoc_laguerre(p, abs(l), 2.*r*r/w0/w0) : 1;
 	double res = 2.*r/w0/w0*amp;
 	if(l!=0) res-= abs(l)*amp/r;
 	if(p!=0){double amp2 = Eo*exp(-r*r/w0/w0)*assoc_laguerre(p-1, abs(l)+1, 2.*r*r/w0/w0); amp2 *= (l!=0) ? pow(r*sqrt(2.)/w0 , abs(l)) : 1; res+=amp2*4*r/w0/w0;}
 	res *= cos(arg)*sin(phi);
 	if(l!=0) res-=l*amp/r*sin(arg)*cos(phi);
 	return res*Envelope(x,t);
 }

 else return 0.;
}

double DerBfx(double x, double y, double z){ //Bx time derivative at (x,y,z)
 if(wave_type == 0) return 0; 
 if(wave_type == 1) return 0;
 if(wave_type == 2){
 	double r = sqrt(y*y+z*z);
	double phi = atan2(z,y);
	double wz = w0*sqrt(1+x*x/zr/zr);
	double R_1 = x / (x*x + zr*zr);
	double N = abs(l) + 2*p;
	double arg = w*t-k*x-k*r*r*R_1/2.-l*phi+(N+1)*atan(x/zr);
 	double amp = Eo*w0*wz*exp(-r*r/wz/wz); amp *= (l!=0)? pow(r*sqrt(2.)/wz, abs(l))*assoc_laguerre(abs(p), abs(l), 2.*r*r/wz/wz) : (p!=0)? assoc_laguerre(abs(p), abs(l), 2.*r*r/wz/wz) : 1;
 	double f1 = (2*r/wz/wz-abs(l)/r)*amp;
 	if(p!=0){double amp2 = Eo*w0/wz*exp(-r*r/wz/wz)*assoc_laguerre(p-1, abs(l)+1, 2.*r*r/wz/wz); amp2 *= (l!=0) ? pow(r*sqrt(2.)/wz , abs(l)) : 1; f1+=amp2*4*r/wz/wz;}
 	f1*=sin(phi);
 	double f2 = -(l/r*cos(phi)+k*r*R_1*sin(phi))*amp;
 	return w*(-f1*sin(arg)+f2*cos(arg))*Envelope(x,t);
 }
 if(wave_type == 3){ 
 	double r = sqrt(y*y+z*z);
	double phi = atan2(z,y);
 	double arg = w*t-k*x-l*phi;
 	double amp = Eo * exp(-r*r/w0/w0); amp *= (l!=0) ? pow(r*sqrt(2.)/w0 , abs(l))*assoc_laguerre(p, abs(l), 2.*r*r/w0/w0) : (p!=0) ? assoc_laguerre(p, abs(l), 2.*r*r/w0/w0) : 1;
 	double res = 2.*r/w0/w0*amp;
 	if(l!=0) res-= abs(l)*amp/r;
 	if(p!=0){double amp2 = Eo*exp(-r*r/w0/w0)*assoc_laguerre(p-1, abs(l)+1, 2.*r*r/w0/w0); amp2 *= (l!=0) ? pow(r*sqrt(2.)/w0 , abs(l)) : 1; res+=amp2*4*r/w0/w0;}
 	res *= -w*sin(arg)*sin(phi);
 	if(l!=0) res-=w*l*amp/r*cos(arg)*cos(phi);
 	return res*Envelope(x,t);
 }
 else return 0.;
}

double Bfy(double x, double y, double z){  //By interpolation to (x,y,z)
 if(wave_type == 0) return sqrt(1-delta*delta)*k*Bo*cos(w*t-k*x)*Envelope(x,t);
 if(wave_type == 1) return 0;
 if(wave_type == 2) return 0;
 if(wave_type == 3) return 0;
 else return 0.;
}

double DerBfy(double x, double y, double z){ //By time derivative at (x,y,z)
 if(wave_type == 0) return -sqrt(1-delta*delta)*w*k*Bo*sin(w*t-k*x)*Envelope(x,t);
 if(wave_type == 1) return 0;
 if(wave_type == 2) return 0;
 if(wave_type == 3) return 0;
 else return 0.;
}

double Bfz(double x, double y, double z){  //Bz interpolation to (x,y,z)
 if(wave_type == 0) return delta*k*Bo*sin(w*t-k*x)*Envelope(x,t);
 if(wave_type == 1){
 	double r = sqrt(y*y+z*z);
 	double wz = w0 * sqrt(1 + x*x/(zr*zr));
 	double R_1 = x / (x*x + zr*zr);
 	double psi = atan(x/zr);
 	return kg*Eo / eta * w0/wz * exp(-r*r/(wz*wz)) * sin(w*t - kg*x - kg*r*r*R_1/2. + psi) * Envelope(x, t);
 }
 if(wave_type == 2){
 	double r = sqrt(y*y+z*z);
	double phi = atan2(z,y);
	double wz = w0*sqrt(1+x*x/zr/zr);
	double R_1 = x / (x*x + zr*zr);
	double N = abs(l) + 2*p;
	double arg = w*t-k*x-k*r*r*R_1/2.-l*phi+(N+1)*atan(x/zr);
 	double amp = Eo*w0*wz*exp(-r*r/wz/wz); amp *= (l!=0)? pow(r*sqrt(2.)/wz, abs(l))*assoc_laguerre(abs(p), abs(l), 2.*r*r/wz/wz) : (p!=0)? assoc_laguerre(abs(p), abs(l), 2.*r*r/wz/wz) : 1;
 	double f1 = (2*r*r/wz/wz-1-abs(l))*amp;
 	if(p!=0){double amp2 = Eo*w0/wz*exp(-r*r/wz/wz)*assoc_laguerre(p-1, abs(l)+1, 2.*r*r/wz/wz); amp2 *= (l!=0) ? pow(r*sqrt(2.)/wz , abs(l)) : 1; f1+=amp2*2*r*r/wz/wz;}
 	f1*=w0*w0*x/zr/zr/wz/wz;
 	double f2 = (k-k*r*r/2*(x*x-zr*zr)/(x*x+zr*zr)/(x*x+zr*zr)-(N+1)*zr/(x*x+zr*zr))*amp;
 	return (f1*cos(arg)+f2*sin(arg))*Envelope(x,t);
 }
 if(wave_type==3){ //eturn 0;
 	double r = sqrt(y*y+z*z);
	double phi = atan2(z,y);
	double arg = w*t-k*x-l*phi;
	double amp=1;
	amp *= Eo;
	amp *= exp(-r*r/w0/w0);
	if(l!=0) amp *= pow(r*sqrt(2.)/w0 , abs(l));
	if (l!=0 || p!=0) amp *= assoc_laguerre(abs(p), abs(l), 2.*r*r/w0/w0); 	
 	double res = k*amp*sin(arg); 
 	return res*Envelope(x,t);
 } 

 else return 0.;
}

double DerBfz(double x, double y, double z){ //Bz time derivative at (x,y,z)
 if(wave_type == 0) return delta*w*k*Bo*cos(w*t-k*x)*Envelope(x,t); 
 if(wave_type == 1){
 	double r = sqrt(y*y+z*z);
 	double wz = w0 * sqrt(1 + x*x/(zr*zr));
 	double R_1 = x / (x*x + zr*zr);
 	double psi = atan(x/zr);
 	return w*kg * Eo / eta * w0/wz * exp(-r*r/(wz*wz)) * cos(w*t - kg*x - kg*r*r*R_1/2. + psi) * Envelope(x, t);
 }
 if(wave_type == 2){
 	double r = sqrt(y*y+z*z);
	double phi = atan2(z,y);
	double wz = w0*sqrt(1+x*x/zr/zr);
	double R_1 = x / (x*x + zr*zr);
	double N = abs(l) + 2*p;
	double arg = w*t-k*x-k*r*r*R_1/2.-l*phi+(N+1)*atan(x/zr);
 	double amp = Eo*w0*wz*exp(-r*r/wz/wz); amp *= (l!=0)? pow(r*sqrt(2.)/wz, abs(l))*assoc_laguerre(abs(p), abs(l), 2.*r*r/wz/wz) : (p!=0)? assoc_laguerre(abs(p), abs(l), 2.*r*r/wz/wz) : 1;
 	double f1 = (2*r*r/wz/wz-1-abs(l))*amp;
 	if(p!=0){double amp2 = Eo*w0/wz*exp(-r*r/wz/wz)*assoc_laguerre(p-1, abs(l)+1, 2.*r*r/wz/wz); amp2 *= (l!=0) ? pow(r*sqrt(2.)/wz , abs(l)) : 1; f1+=amp2*2*r*r/wz/wz;}
 	f1*=w0*w0*x/zr/zr/wz/wz;
 	double f2 = (k-k*r*r/2*(x*x-zr*zr)/(x*x+zr*zr)/(x*x+zr*zr)-(N+1)*zr/(x*x+zr*zr))*amp;
 	return w*(-f1*sin(arg)+f2*cos(arg))*Envelope(x,t);
 }
 if(wave_type==3){
 	double r = sqrt(y*y+z*z);
	double phi = atan2(z,y);
	double arg = w*t-k*x-l*phi;
	double amp=1;
	amp *= Eo;
	amp *= exp(-r*r/w0/w0);
	if(l!=0) amp *= pow(r*sqrt(2.)/w0 , abs(l));
	if (l!=0 || p!=0) amp *= assoc_laguerre(abs(p), abs(l), 2.*r*r/w0/w0); 	
 	double res = k*w*amp*cos(arg); 
 	return res*Envelope(x,t);
 } 
 else return 0.;
}

//----------------------------------------------------------------------------


//CHECK TIME DERIVATIVES
double dt_Efx(double x, double y, double z){
	double h=1e-5;
	t+=h/2; double a1=Efx(x,y,z);
	t-=h; double a2=Efx(x,y,z);
	t+=h/2; double res= (a1-a2)/h;
	return res;
}
double dt_Efy(double x, double y, double z){
	double h=1e-5;
	t+=h/2; double a1=Efy(x,y,z);
	t-=h; double a2=Efy(x,y,z);
	t+=h/2; double res= (a1-a2)/h;
	return res;
}
double dt_Efz(double x, double y, double z){
	double h=1e-5;
	t+=h/2; double a1=Efz(x,y,z);
	t-=h; double a2=Efz(x,y,z);
	t+=h/2; double res= (a1-a2)/h;
	return res;
}
double dt_Bfx(double x, double y, double z){
	double h=1e-5;
	t+=h/2; double a1=Bfx(x,y,z);
	t-=h; double a2=Bfx(x,y,z);
	t+=h/2; double res= (a1-a2)/h;
	return res;
}
double dt_Bfy(double x, double y, double z){
	double h=1e-5;
	t+=h/2; double a1=Bfy(x,y,z);
	t-=h; double a2=Bfy(x,y,z);
	t+=h/2; double res= (a1-a2)/h;
	return res;
}
double dt_Bfz(double x, double y, double z){
	double h=1e-5;
	t+=h/2; double a1=Bfz(x,y,z);
	t-=h; double a2=Bfz(x,y,z);
	t+=h/2; double res= (a1-a2)/h;
	return res;
}



//CHECK DIV,ROT---------------------------------------------------------

double divEf(double x, double y, double z){
	double h=1e-5;
	double a1 = (Efx(x+h/2,y,z)-Efx(x-h/2,y,z))/h;
	double a2 = (Efy(x,y+h/2,z)-Efy(x,y-h/2,z))/h;
	double a3 = (Efz(x,y,z+h/2)-Efz(x,y,z-h/2))/h;
	return a1+a2+a3;
}
double divBf(double x, double y, double z){
	double h=1e-5;
	double a1 = (Bfx(x+h/2,y,z)-Bfx(x-h/2,y,z))/h;
	double a2 = (Bfy(x,y+h/2,z)-Bfy(x,y-h/2,z))/h;
	double a3 = (Bfz(x,y,z+h/2)-Bfz(x,y,z-h/2))/h;
	return a1+a2+a3;
}
double* curlEf(double x, double y, double z){
	double h=1e-5;
	double* res = new double[3];
	res[0] = (Efz(x,y+h/2,z)-Efz(x,y-h/2,z))/h - (Efy(x,y,z+h/2)-Efy(x,y,z-h/2))/h;
	res[1] = (Efx(x,y,z+h/2)-Efx(x,y,z-h/2))/h - (Efz(x+h/2,y,z)-Efz(x-h/2,y,z))/h;
	res[2] = (Efy(x+h/2,y,z)-Efy(x-h/2,y,z))/h - (Efx(x,y+h/2,z)-Efx(x,y-h/2,z))/h;
	return res;
}
double* curlBf(double x, double y, double z){
	double h=1e-5;
	double* res = new double[3];
	res[0] = (Bfz(x,y+h/2,z)-Bfz(x,y-h/2,z))/h - (Bfy(x,y,z+h/2)-Bfy(x,y,z-h/2))/h;
	res[1] = (Bfx(x,y,z+h/2)-Bfx(x,y,z-h/2))/h - (Bfz(x+h/2,y,z)-Bfz(x-h/2,y,z))/h;
	res[2] = (Bfy(x+h/2,y,z)-Bfy(x-h/2,y,z))/h - (Bfx(x,y+h/2,z)-Bfx(x,y-h/2,z))/h;
	return res;
}

//--------------------------------------------------------------------------


double teste(double*x,double*par){


	//CHECK DIV,ROT
	//return divEf(par[0],x[0],x[1]);
	//return divBf(par[0],x[0],x[1]);
	//check curlE
	//double a1,a2,a3;
	//a1 = pow(curlEf(par[0],x[0],x[1])[0]+DerBfx(par[0],x[0],x[1]),2);
	//a2 = pow(curlEf(par[0],x[0],x[1])[1]+DerBfy(par[0],x[0],x[1]),2);
	//a3 = pow(curlEf(par[0],x[0],x[1])[2]+DerBfz(par[0],x[0],x[1]),2);
	//return sqrt(a1+a2+a3);
	//check curlB
	//double a1,a2,a3;
	//a1 = pow( curlBf(par[0],x[0],x[1])[0]-DerEfx(par[0],x[0],x[1]) ,2);
	//a2 = pow( curlBf(par[0],x[0],x[1])[1]-DerEfy(par[0],x[0],x[1]) ,2);
	//a3 = pow( curlBf(par[0],x[0],x[1])[2]-DerEfz(par[0],x[0],x[1]) ,2);
	//return sqrt(a1+a2+a3);

	//CHECK TIME DERIVATIVES
	//double b1 = DerEfx(par[0],x[0],x[1])-dt_Efx(par[0],x[0],x[1]);
	//double b2 = DerEfy(par[0],x[0],x[1])-dt_Efy(par[0],x[0],x[1]);
	//ouble b3 = DerEfz(par[0],x[0],x[1])-dt_Efz(par[0],x[0],x[1]);
	//double b4 = DerBfx(par[0],x[0],x[1])-dt_Bfx(par[0],x[0],x[1]);
	//double b5 = DerBfy(par[0],x[0],x[1])-dt_Bfy(par[0],x[0],x[1]);
	//double b6 = DerBfz(par[0],x[0],x[1])-dt_Bfz(par[0],x[0],x[1]);
	//return sqrt(b1*b1+b2*b2+b3*b3+b4*b4+b5*b5+b6*b6);

	//CHECK SOME COMPONENT
	return Efy(par[0],x[0],x[1]);
}















int main(int argc, char** argv){

	FILE*foo;
	foo=fopen("../../rk4_cooling/InputToBatch.txt","r");
	fscanf(foo,"%s %lf %lf %lf %lf %lf %lf %lf %lf %lli %i %lf %i %lf %lf %lf %lf %lf %lf %lf %lf %i %i", 
				trash, &x01, &x02, &x03, &p01, &p02, &p03, &kdamp, &T, &N, &pri, &dx, 
				&wave_type, &tfwhm, &stable, &Eo, &delta, &w0, &lambda, &n, &eta, &lr, &pr);
	fclose(foo);

	
	k=1;
	kg = 2*M_PI*n/lambda;
	if(wave_type != 0) k = kg;
	zr = M_PI*w0*w0*n/lambda;
	w=k;


	cout<<"WAVE_TYPE = "<<wave_type<<endl;
	cout<<"Ao = "<<Eo<<endl;
	cout<<"zr = "<<zr<<endl;
	cout<<"delta = "<<delta<<endl;
	cout<<"w,k = "<<w<<endl;



	double range=1;

	int n_l=6; n_l=1;
	int n_p=4; n_p=1;
	int side=600;

	auto f1 = new TF2**[n_p];
	auto t1 = new TLatex**[n_p];

	TApplication* MyRootApp;
	MyRootApp = new TApplication("MyRootApp", NULL, NULL);

	auto c1 = new TCanvas("c1", "", n_l*side, n_p*side);
	c1->Divide(n_l,n_p);
	c1->SetRightMargin(0.0);
	c1->SetLeftMargin(0.0);
	c1->SetBottomMargin(0.0);
	c1->SetTopMargin(0.0);


	gStyle->SetPalette(kBird);


	double MINI=0, MAXI=0;

	cout<<endl;
	t = 0;
	double ti=0; double dt=0.05; double tempo=50;
	while(ti<tempo){
	for(int pi=0; pi<n_p; pi++){
		f1[pi] = new TF2*[n_l];
		t1[pi] = new TLatex*[n_l];
		for(int li=0; li<n_l; li++){
			f1[pi][li] = new TF2("",teste,-20,20,-20,20,1);
			//p=pi; l=li;
			p=pr; l=lr;
			t=ti;
			f1[pi][li]->SetParameter(0,0); //x
			range = max(abs(MAXI), abs(MINI));
			f1[pi][li]->SetMinimum(-range); 
			f1[pi][li]->SetMaximum(range); 
			f1[pi][li]->SetNpx(500);
			f1[pi][li]->SetNpy(500);
			//cout<<f1[pi][li]->GetMaximum()<<"\t"<<f1[pi][li]->GetMinimum()<<endl;
			c1->cd(1+li+pi*n_l);
			
			//c1->cd(1+li+pi*n_l)->SetRightMargin(0.0);
			//c1->cd(1+li+pi*n_l)->SetLeftMargin(0.0);
			//c1->cd(1+li+pi*n_l)->SetBottomMargin(0.0);
			//c1->cd(1+li+pi*n_l)->SetTopMargin(0.0);
			
			f1[pi][li]->Draw("colz");
				
			//t1[pi][li] = new TLatex(-10,16,("#font[132]{p = "+to_string((int)p)+" | l = "+to_string((int)l)+"}").c_str());
			//t1[pi][li]->SetTextSize(0.12);
			t1[pi][li] = new TLatex(-7,16,("#font[132]{p = "+to_string((int)p)+" | l = "+to_string((int)l)+"}").c_str());
			t1[pi][li]->SetTextSize(0.07);
			t1[pi][li]->Draw("SAME");
			
			//c1->cd(1+li+pi*n_l)->Update();
			gPad->Update();

			//c1->SaveAs("plot.png");

			double maxi = f1[0][0]->GetMaximum(), mini=f1[0][0]->GetMinimum();
			cout<<"\33[2K\r";
			cout<<"TIME: "<<t<<"\t\t"<<"pl="<<p<<l<<"\t\t( "<<MINI<<" , "<<MAXI<<" )";
			if(max(abs(mini), abs(maxi))>1e-7) cout<<"\t\t***NOT ZERO***";
			if(maxi>MAXI) MAXI = maxi;
			if(mini<MINI) MINI = mini;
			cout<<flush;

		}
	}
	
	c1->Update();
	ti+=dt;
	if(situacao()>97){ cout<<endl<<endl<<"DEATH BY RAM"<<endl<<endl; return 1;}
	}

	cout<<endl<<endl;

	return 0;
}

