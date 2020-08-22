
//#include "laguerre_gaussian.cpp"

#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>

#include <TF2.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TApplication.h>

using namespace std;

complex<double> imu(0,1);


char trash[128];
	double x01,x02,x03,p01,p02,p03, kdamp, T;
	long long int N; int pri;
	double dx, dy; int wave_type;
	double tfwhm, stable, Eo, delta, w0, lambda, n, eta;
	int lr,pr;
	double k,kg,zr,w,Bo;


double der_coef4[5] = {1./12., -2./3., 0., 2./3., -1./12.};


double Envelope(double x, double t){return 1;}

double Efx(double x, double y, double z, double p, double l, double t){  //Ex interpolation to (x,y,z)
 if(wave_type == 0) return 0; 
 if(wave_type == 1) return 0;
 if(wave_type == 2) return 0;
 if(wave_type == 3) return 0;
 else return 0.;
}

double DerEfx(double x, double y, double z, double p, double l, double t){ //Ex time derivative at (x,y,z)
 if(wave_type == 0) return 0; 
 if(wave_type == 1) return 0;
 if(wave_type == 2) return 0;
 if(wave_type == 3) return 0;
 else return 0.;
}

double Efy(double x, double y, double z, double p, double l, double t){  //Ey interpolation to (x,y,z)
 if(wave_type == 0) return delta*w*Eo*sin(w*t-k*x)*Envelope(x,t);
 if(wave_type == 1){
 	double r = sqrt(y*y+z*z);
 	double wz = w0 * sqrt(1 + x*x/(zr*zr));
 	double R_1 = x / (x*x + zr*zr);
 	double psi = atan(x/zr); 
 	psi=0;
 	return Eo * w0/wz * exp(-r*r/(wz*wz)) * sin(w*t - kg*x - kg*r*r*R_1/2 + psi) * Envelope(x, t);
 }
 if(wave_type == 2){
 	double r = sqrt(y*y+z*z);
	double phi = atan2(z,y);
	double wz = w0*sqrt(1+x*x/zr/zr);

 	double res=1;

 	res *= w;

 	res *= Eo*w0/wz;
 	res *= pow(r*sqrt(2.)/wz, abs(l));
 	res *= assoc_laguerre(abs(p), abs(l), 2.*r*r/wz/wz);
 	res *= exp(-r*r/wz/wz);

 	double arg = w*t-k*x-k*r*r*x/2./(x*x+zr*zr)-l*phi+(2.*(double)p+(double)abs(l)+1.)*atan(x/zr);
 	res *= sin(arg);

 	return res*Envelope(x,t);
 }
 if(wave_type==3){
 	double r = sqrt(y*y+z*z);
	double phi = atan2(z,y);
 	double amp = Eo * exp(-r*r/w0/w0) * pow(r*sqrt(2.)/w0 , abs(l)) * assoc_laguerre(abs(p), abs(l), 2.*r*r/w0/w0);
 	double arg = w*t-k*x-l*phi;
 	double res = w*amp*sin(arg); 
 	return res*Envelope(x,t);
 } 
 else return 0.;
}

double DerEfy(double x, double y, double z, double p, double l, double t){ //Ey time derivative at (x,y,z)
 if(wave_type == 0) return delta*w*w*Eo*cos(w*t-k*x)*Envelope(x,t); 
 if(wave_type == 1){
 	double r = sqrt(y*y+z*z);
 	double wz = w0 * sqrt(1 + x*x/(zr*zr));
 	double R_1 = x / (x*x + zr*zr);
 	double psi = atan(x/zr);
 	psi=0;
 	return w * Eo * w0/wz * exp(-r*r/(wz*wz)) * cos(w*t - kg*x - kg*r*r*R_1/2 + psi) * Envelope(x, t);
 }
 if(wave_type == 2){
 	double r = sqrt(y*y+z*z);
	double phi = atan2(z,y);
	double wz = w0*sqrt(1+x*x/zr/zr);

 	double res=1;

 	res *= w*w;

 	res *= Eo*w0/wz;
 	res *= pow(r*sqrt(2.)/wz, abs(l));
 	res *= assoc_laguerre(abs(p), abs(l), 2.*r*r/wz/wz);
 	res *= exp(-r*r/wz/wz);

 	double arg = w*t-k*x-k*r*r*x/2./(x*x+zr*zr)-l*phi+(2.*(double)p+(double)abs(l)+1.)*atan(x/zr);
 	res *= cos(arg);

 	return res*Envelope(x,t);
 }
 if(wave_type==3){
 	double r = sqrt(y*y+z*z);
	double phi = atan2(z,y);
 	double amp = Eo * exp(-r*r/w0/w0) * pow(r*sqrt(2.)/w0 , abs(l)) * assoc_laguerre(abs(p), abs(l), 2.*r*r/w0/w0);
 	double arg = w*t-k*x-l*phi;
 	double res = w*w*amp*cos(arg); 
 	return res*Envelope(x,t);;
 } 
 else return 0.;
}

double Efz(double x, double y, double z, double p, double l, double t){  //Ez interpolation to (x,y,z)
 if(wave_type == 0) return -sqrt(1-delta*delta)*w*Eo*cos(w*t-k*x)*Envelope(x,t); 
 if(wave_type == 1) return 0;
 if(wave_type == 2) return 0;
 if(wave_type == 3) return 0;
 else return 0.;
}

double DerEfz(double x, double y, double z, double p, double l, double t){ //Ez time derivative at (x,y,z)
 if(wave_type == 0) return sqrt(1-delta*delta)*w*w*Eo*sin(w*t-k*x)*Envelope(x,t);
 if(wave_type == 1) return 0;
 if(wave_type == 2) return 0;
 if(wave_type == 3) return 0;
 else return 0.;
}

double Bfx( double x, double y, double z, double p, double l, double t){  //Bx interpolation to (x,y,z)
 if(wave_type == 0) return 0; 
 if(wave_type == 1) return 0;
 if(wave_type == 2) return 0;
 if(wave_type == 3){return 0;
 	double r = sqrt(y*y+z*z);
	double phi = atan2(z,y);
 	double amp = Eo * exp(-r*r/w0/w0) * pow(r*sqrt(2.)/w0 , abs(l)) * assoc_laguerre(abs(p), abs(l), 2.*r*r/w0/w0);
 	double arg = w*t-k*x-l*phi;

 	double res = 2.*r*r/w0/w0 * cos(arg)*sin(phi);
 	if(p!=0){
 		double amp3 = Eo * exp(-r*r/w0/w0) * pow(r*sqrt(2.)/w0 , abs(l)) * assoc_laguerre(abs(p)-1, abs(l)+1, 2.*r*r/w0/w0); //com Laguerre_(p-1)_(l+1)
 		amp += 2.*amp3;
 	}
 	res *= amp;
 	if(l!=0){
 		double amp2 = Eo * exp(-r*r/w0/w0) * pow(sqrt(2.)/w0 , abs(l))*pow(r, abs(l)-1) * assoc_laguerre(abs(p), abs(l), 2.*r*r/w0/w0); //com r^(|l|-1)
 		res += amp2*(-abs(l)*cos(arg)*sin(phi) -l*sin(arg)*cos(phi));
 	}

 	return res*Envelope(x,t);
 }

 else return 0.;
}

double DerBfx(double x, double y, double z, double p, double l, double t){ //Bx time derivative at (x,y,z)
 if(wave_type == 0) return 0; 
 if(wave_type == 1) return 0;
 if(wave_type == 2) return 0;
 if(wave_type == 3){return 0;
 	double r = sqrt(y*y+z*z);
	double phi = atan2(z,y);
 	double amp = Eo * exp(-r*r/w0/w0) * pow(r*sqrt(2.)/w0 , abs(l)) * assoc_laguerre(abs(p), abs(l), 2.*r*r/w0/w0);
 	double arg = w*t-k*x-l*phi;

 	double res = 2.*r*r/w0/w0 * (-w*sin(arg))*sin(phi);
 	if(p!=0){
 		double amp3 = Eo * exp(-r*r/w0/w0) * pow(r*sqrt(2.)/w0 , abs(l)) * assoc_laguerre(abs(p)-1, abs(l)+1, 2.*r*r/w0/w0); //com Laguerre_(p-1)_(l+1)
 		amp += 2.*amp3;
 	}
 	res *= amp;
 	if(l!=0){
 		double amp2 = Eo * exp(-r*r/w0/w0) * pow(sqrt(2.)/w0 , abs(l))*pow(r, abs(l)-1) * assoc_laguerre(abs(p), abs(l), 2.*r*r/w0/w0); //com r^(|l|-1)
 		res += amp2*(-abs(l)*(-w*sin(arg))*sin(phi) -l*(w*cos(arg))*cos(phi));
 	}

 	return res*Envelope(x,t);
 }
 else return 0.;
}

double Bfy(double x, double y, double z, double p, double l, double t){  //By interpolation to (x,y,z)
 if(wave_type == 0) return sqrt(1-delta*delta)*k*Bo*cos(w*t-k*x)*Envelope(x,t);
 if(wave_type == 1) return 0;
 if(wave_type == 2) return 0;
 if(wave_type == 3) return 0;
 else return 0.;
}

double DerBfy(double x, double y, double z, double p, double l, double t){ //By time derivative at (x,y,z)
 if(wave_type == 0) return -sqrt(1-delta*delta)*w*k*Bo*sin(w*t-k*x)*Envelope(x,t);
 if(wave_type == 1) return 0;
 if(wave_type == 2) return 0;
 if(wave_type == 3) return 0;
 else return 0.;
}

double Bfz(double x, double y, double z, double p, double l, double t){  //Bz interpolation to (x,y,z)
 if(wave_type == 0) return delta*k*Bo*sin(w*t-k*x)*Envelope(x,t);
 if(wave_type == 1){
 	double r = sqrt(y*y+z*z);
 	double wz = w0 * sqrt(1 + x*x/(zr*zr));
 	double R_1 = x / (x*x + zr*zr);
 	double psi = atan(x/zr);
 	psi=0;
 	return Eo / eta * w0/wz * exp(-r*r/(wz*wz)) * sin(w*t - kg*x - kg*r*r*R_1/2 + psi) * Envelope(x, t);
 }
 if(wave_type == 2){
 	double r = sqrt(y*y+z*z);
	double phi = atan2(z,y);
	double wz = w0*sqrt(1+x*x/zr/zr);

	double res;
	double A1_re=0; double A1_im=0; double A2_re=0; double A2_im=0; double A3_re=0; double A3_im=0;
 	
 	double aux1 = Eo*w0/wz *  pow(r*sqrt(2.)/wz, abs(l)) * assoc_laguerre(abs(p), abs(l), 2.*r*r/wz/wz) * exp(-r*r/wz/wz);
 	double arg = -k*r*r*x/2./(x*x+zr*zr)-l*phi+(2.*(double)p+(double)abs(l)+1.)*atan(x/zr);
 	
 	double A1 = -w0*w0*x/zr/zr/wz/wz*aux1;
 	if(l!=0) A1 *= (1.+(double)abs(l));
 	A1_re = A1 * cos(arg);
 	A1_im = A1 * sin(arg);

 	if(p!=0){
 		double aux2 = Eo*w0/wz *  pow(r*sqrt(2.)/wz, abs(l)) * assoc_laguerre(abs(p)-1, abs(l)+1, 2.*r*r/wz/wz) * exp(-r*r/wz/wz);
 		double A2 = 4*r*r*w0*w0*x/zr/zr/wz/wz/wz/wz * aux2;
 		A2_re = A2 * cos(arg);
 		A2_im = A2 * sin(arg);
 	}

 	double A3 = aux1*(-k*r*r/2.*(zr*zr-x*x)/(x*x+zr*zr)/(x*x+zr*zr)+(2.*(double)p+(double)abs(l)+1.)*zr/(x*x+zr*zr)-k);
 	A3_re = -A3*sin(arg);
 	A3_im = A3*cos(arg);

	res = cos(w*t-k*x)*(A1_re+A2_re+A3_re) - sin(w*t-k*x)*(A1_im+A2_im+A3_im);
 	return res*Envelope(x,t);
 }
 if(wave_type==3){
 	double r = sqrt(y*y+z*z);
	double phi = atan2(z,y);
 	double amp = Eo * exp(-r*r/w0/w0) * pow(r*sqrt(2.)/w0 , abs(l)) * assoc_laguerre(abs(p), abs(l), 2.*r*r/w0/w0);
 	double arg = w*t-k*x-l*phi;
 	double res = k*amp*sin(arg); 
 	return res*Envelope(x,t);
 } 

 else return 0.;
}

double DerBfz(double x, double y, double z, double p, double l, double t){ //Bz time derivative at (x,y,z)
 if(wave_type == 0) return delta*w*k*Bo*cos(w*t-k*x)*Envelope(x,t); 
 if(wave_type == 1){
 	double r = sqrt(y*y);
 	double wz = w0 * sqrt(1 + x*x/(zr*zr));
 	double R_1 = x / (x*x + zr*zr);
 	double psi = atan(x/zr);
 	psi=0;
 	return w * Eo / eta * w0/wz * exp(-r*r/(wz*wz)) * cos(w*t - kg*x - kg*r*r*R_1/2 + psi) * Envelope(x, t);
 }
 if(wave_type == 2){
 	double r = sqrt(y*y+z*z);
	double phi = atan2(z,y);
	double wz = w0*sqrt(1+x*x/zr/zr);

	double res;
	double A1_re=0; double A1_im=0; double A2_re=0; double A2_im=0; double A3_re=0; double A3_im=0;
 	
 	double aux1 = Eo*w0/wz *  pow(r*sqrt(2.)/wz, abs(l)) * assoc_laguerre(abs(p), abs(l), 2.*r*r/wz/wz) * exp(-r*r/wz/wz);
 	double arg = -k*r*r*x/2./(x*x+zr*zr)-l*phi+(2.*(double)p+(double)abs(l)+1.)*atan(x/zr);
 	
 	double A1 = -w0*w0*x/zr/zr/wz/wz*aux1;
 	if(l!=0) A1 *= (1.+(double)abs(l));
 	A1_re = A1 * cos(arg);
 	A1_im = A1 * sin(arg);

 	if(p!=0){
 		double aux2 = Eo*w0/wz *  pow(r*sqrt(2.)/wz, abs(l)) * assoc_laguerre(abs(p)-1, abs(l)+1, 2.*r*r/wz/wz) * exp(-r*r/wz/wz);
 		double A2 = 4*r*r*w0*w0*x/zr/zr/wz/wz/wz/wz * aux2;
 		A2_re = A2 * cos(arg);
 		A2_im = A2 * sin(arg);
 	}

 	double A3 = aux1*(-k*r*r/2.*(zr*zr-x*x)/(x*x+zr*zr)/(x*x+zr*zr)+(2.*(double)p+(double)abs(l)+1.)*zr/(x*x+zr*zr)-k);
 	A3_re = -A3*sin(arg);
 	A3_im = A3*cos(arg);

	res = -w*sin(w*t-k*x)*(A1_re+A2_re+A3_re) - w*cos(w*t-k*x)*(A1_im+A2_im+A3_im);
 	return res*Envelope(x,t);
 }
 if(wave_type==3){
 	double r = sqrt(y*y+z*z);
	double phi = atan2(z,y);
 	double amp = Eo * exp(-r*r/w0/w0) * pow(r*sqrt(2.)/w0 , abs(l)) * assoc_laguerre(abs(p), abs(l), 2.*r*r/w0/w0);
 	double arg = w*t-k*x-l*phi;
 	double res = k*w*amp*cos(arg); 
 	return res*Envelope(x,t);
 } 
 else return 0.;
}







complex<double> upl(double r, double phi, double z, double p, double l){

	complex<double> res(1,0);
	double wz = w0*sqrt(1+z*z/zr/zr);

	res *= Eo*w0/wz;
	res *= pow(r*sqrt(2.)/wz, abs(l));
	res *= assoc_laguerre((int)abs(p), (int)abs(l), 2.*r*r/wz/wz);
	res *= exp(-r*r/wz/wz);
	double aux = -kg*r*r*z/2./(z*z+zr*zr)-l*phi+(2.*p+abs(l)+1.)*atan(z/zr);
	res *= exp(imu*aux);

	return res;
}
complex<double> Alg(double r, double phi, double z, double p, double l, double t){
	return upl(r,phi,z,p,l)*exp(imu*(w*t-kg*z));
}

complex<double> electric_field(double x, double y, double z, double p, double l, double t){
	complex<double> res;
	double r = sqrt(y*y+z*z);
	double phi = atan2(z,y);
	res = -imu*w*upl(r,phi,x,p,l)*exp(imu*(w*t-kg*x));
	return res;
}

complex<double> magnetic_field(double x, double y, double z, double p, double l, double t){
	double r = sqrt(y*y+z*z);
	double phi = atan2(z,y);
	double h=1e-5;
 	complex<double> res=0;
 	for(int i=0; i<5; i++) res += 1./h * der_coef4[i] * Alg(r,phi,x+(double)(i-2)*h,p,l,t);   
 	return res;
}

complex<double> magnetic_field_luis(double x, double y, double z, double p, double l, double t){
	
	complex<double> res(0,0);

	double r = sqrt(y*y+z*z);
	double phi = atan2(z,y);
	double wz = w0*sqrt(1+x*x/zr/zr);

	//adicionar mood = du/dz-iku
	res += -w0*w0*x/zr/zr/wz/wz * upl(r,phi,x,p,l);
	if((int)l != 0) res *= 1.+abs(l);
	if((int)p != 0){
		complex<double> upl_gay(1,0);
		upl_gay *= Eo*w0/wz;
		upl_gay *= pow(r*sqrt(2.)/wz, abs(l));
		upl_gay *= assoc_laguerre((int)abs(p)-1, (int)abs(l)+1, 2.*r*r/wz/wz); //<<<===
		upl_gay *= exp(-r*r/wz/wz);
		double aux = -kg*r*r*x/2./(x*x+zr*zr)-l*phi+(2.*p+abs(l)+1.)*atan(x/zr);
		upl_gay *= exp(imu*aux);

		res += 4*r*r*w0*w0*x/zr/zr/wz/wz/wz/wz*upl_gay;
	}
	double aux = -kg*r*r/2.*(zr*zr-x*x)/(x*x+zr*zr)/(x*x+zr*zr)+(2.*p+abs(l)+1.)*zr/(x*x+zr*zr)-kg;
	res += imu*aux*upl(r,phi,x,p,l);

	res *= exp(imu*(w*t-kg*x));

	return res;
}


double electric_field_trig(double x, double y, double z, double p, double l, double t){
	double r = sqrt(y*y+z*z);
	double phi = atan2(z,y);
	double wz = w0*sqrt(1+x*x/zr/zr);

 	double res=1;

 	res *= w;
 	res *= Eo*w0/wz;
 	res *= pow(r*sqrt(2.)/wz, abs(l));
 	res *= assoc_laguerre(abs(p), abs(l), 2.*r*r/wz/wz);
 	res *= exp(-r*r/wz/wz);

 	double arg = w*t-k*x-k*r*r*x/2./(x*x+zr*zr)-l*phi+(2.*(double)p+(double)abs(l)+1.)*atan(x/zr);
 	res *= sin(arg);

 	return res;
}
double magnetic_field_trig(double x, double y, double z, double p, double l, double t){
	double r = sqrt(y*y+z*z);
	double phi = atan2(z,y);
	double wz = w0*sqrt(1+x*x/zr/zr);

	double res;
	double A1_re=0; double A1_im=0; double A2_re=0; double A2_im=0; double A3_re=0; double A3_im=0;
 	
 	double aux1 = Eo*w0/wz *  pow(r*sqrt(2.)/wz, abs(l)) * assoc_laguerre(abs(p), abs(l), 2.*r*r/wz/wz) * exp(-r*r/wz/wz);
 	double arg = -k*r*r*x/2./(x*x+zr*zr)-l*phi+(2.*(double)p+(double)abs(l)+1.)*atan(x/zr);
 	
 	double A1 = -w0*w0*x/zr/zr/wz/wz*aux1;
 	if(l!=0) A1 *= (1.+(double)abs(l));
 	A1_re = A1 * cos(arg);
 	A1_im = A1 * sin(arg);

 	if(p!=0){
 		double aux2 = Eo*w0/wz *  pow(r*sqrt(2.)/wz, abs(l)) * assoc_laguerre(abs(p)-1, abs(l)+1, 2.*r*r/wz/wz) * exp(-r*r/wz);
 		double A2 = 4*r*r*w0*w0*x/zr/zr/wz/wz/wz/wz * aux2;
 		A2_re = A2 * cos(arg);
 		A2_im = A2 * sin(arg);
 	}

 	double A3 = aux1*(-k*r*r/2.*(zr*zr-x*x)/(x*x+zr*zr)/(x*x+zr*zr)+(2.*(double)p+(double)abs(l)+1.)*zr/(x*x+zr*zr)-k);
 	A3_re = -A3*sin(arg);
 	A3_im = A3*cos(arg);

	res = cos(w*t-k*x)*(A1_re+A2_re+A3_re) - sin(w*t-k*x)*(A1_im+A2_im+A3_im);
 	return res;
}



double funBluis(double*x,double*par){
	return magnetic_field_luis(par[0],x[0],x[1],par[1],par[2],par[3]).real();
}

double funB(double*x,double*par){
	return magnetic_field(par[0],x[0],x[1],par[1],par[2],par[3]).real();
}

double funE(double*x,double*par){
	return electric_field(par[0],x[0],x[1],par[1],par[2],par[3]).real();
}

double funEtrig(double*x,double*par){
	return electric_field_trig(par[0],x[0],x[1],par[1],par[2],par[3]);
}

double funBtrig(double*x,double*par){
	return magnetic_field_trig(par[0],x[0],x[1],par[1],par[2],par[3]);
}

double teste(double*x,double*par){
	wave_type=1;
	double a1 = Bfz(par[0],x[0],x[1],par[1],par[2],par[3]);
	wave_type=2;
	double a2 = Bfz(par[0],x[0],x[1],par[1],par[2],par[3]);
	return a1-a2;
}


double teste1(double*x,double*par){
	wave_type=1;
	return Bfz(par[0],x[0],0.,par[1],par[2],par[3]);
}
double teste2(double*x,double*par){
	wave_type=3;
	return Bfz(par[0],x[0],0.,par[1],par[2],par[3]);
}



int main(int argc, char** argv){

	FILE*foo;
	foo=fopen("../../rk4_cooling/InputToBatch.txt","r");
	fscanf(foo,"%s %lf %lf %lf %lf %lf %lf %lf %lf %lli %i %lf %lf %i %lf %lf %lf %lf %lf %lf %lf %lf %i %i", 
				trash, &x01, &x02, &x03, &p01, &p02, &p03, &kdamp, &T, &N, &pri, &dx, &dy, 
				&wave_type, &tfwhm, &stable, &Eo, &delta, &w0, &lambda, &n, &eta, &lr, &pr);
	fclose(foo);

	
	k=1;
	kg = 2*M_PI*n/lambda;
	if(wave_type != 0) k = kg;
	zr = M_PI*w0*w0*n/lambda;
	w=1;


	cout<<"WAVE_TYPE="<<wave_type<<endl;
	cout<<"Eo="<<Eo<<endl;
	cout<<"zr="<<zr<<endl;
	cout<<"delta="<<delta<<endl;




	TApplication* MyRootApp;
	MyRootApp = new TApplication("MyRootApp", NULL, NULL);

	auto c1 = new TCanvas("c1", "", 500,500);

	auto f1 = new TF1*[2];
	f1[0] = new TF1("f1", teste1, -20,20,4);
	f1[1] = new TF1("f2", teste2, -20,20,4);

	for(int i=0; i<2; i++){
		f1[i]->SetParameter(0,0);
		f1[i]->SetParameter(2,0);
		f1[i]->SetParameter(1,0);
	}

	f1[0]->SetLineColor(kRed);
	f1[1]->SetLineColor(kBlue);

	f1[0]->SetMinimum(-25);
	f1[0]->SetMaximum(25);


	double t=0, dt=0.01, tempo=100;
	while(t<tempo){

		f1[0]->SetParameter(3,t);
		f1[1]->SetParameter(3,t);

		f1[0]->Draw();	
		f1[1]->Draw("same");
			
		c1->Update();

		t+=dt;
	}

	c1->SaveAs("Plot.png");


	return 0;
	
}
