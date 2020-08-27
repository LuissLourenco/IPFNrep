
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
 	return w*Eo * w0/wz * exp(-r*r/(wz*wz)) * sin(w*t - kg*x - kg*r*r*R_1/2. + psi) * Envelope(x, t);
 }
 if(wave_type == 2){
 	double r = sqrt(y*y+z*z);
	double phi = atan2(z,y);
	double wz = w0*sqrt(1+x*x/zr/zr);
	double R_1 = x / (x*x + zr*zr);
	double N = abs(l) + 2*p;
	double arg = w*t-k*x-k*r*r*R_1/2.-l*phi+(N+1)*atan(x/zr);
 	double res=1;

 	res *= w;
 	res *= Eo*w0/wz;
 	if(l!=0) res *= pow(r*sqrt(2.)/wz, abs(l));
 	if(l!=0 && p!=0) res *= assoc_laguerre(abs(p), abs(l), 2.*r*r/wz/wz);
 	res *= exp(-r*r/wz/wz); 	
 	res *= sin(arg);
 	return res*Envelope(x,t);
 }
 if(wave_type==3){
 	double r = sqrt(y*y+z*z);
	double phi = atan2(z,y);
	double arg = w*t-k*x-l*phi;
	double amp=1;
	amp *= Eo;
	amp *= exp(-r*r/w0/w0);
	if(l!=0) amp *= pow(r*sqrt(2.)/w0 , abs(l));
	if (l!=0 && p!=0) amp *= assoc_laguerre(abs(p), abs(l), 2.*r*r/w0/w0); 	
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
 	return w*w * Eo * w0/wz * exp(-r*r/(wz*wz)) * cos(w*t - kg*x - kg*r*r*R_1/2. + psi) * Envelope(x, t);
 }
 if(wave_type == 2){
 	double r = sqrt(y*y+z*z);
	double phi = atan2(z,y);
	double wz = w0*sqrt(1+x*x/zr/zr);
	double R_1 = x / (x*x + zr*zr);
	double N = abs(l) + 2*p;
	double arg = w*t-k*x-k*r*r*R_1/2.-l*phi+(N+1)*atan(x/zr);
 	double res=1;

 	res *= w*w;
 	res *= Eo*w0/wz;
 	if(l!=0) res *= pow(r*sqrt(2.)/wz, abs(l));
 	if(l!=0 && p!=0) res *= assoc_laguerre(abs(p), abs(l), 2.*r*r/wz/wz);
 	res *= exp(-r*r/wz/wz); 	
 	res *= cos(arg);
 	return res*Envelope(x,t);
 }
 if(wave_type==3){
 	double r = sqrt(y*y+z*z);
	double phi = atan2(z,y);
	double arg = w*t-k*x-l*phi;
	double amp=1;
	amp *= Eo;
	amp *= exp(-r*r/w0/w0);
	if(l!=0) amp *= pow(r*sqrt(2.)/w0 , abs(l));
	if (l!=0 && p!=0) amp *= assoc_laguerre(abs(p), abs(l), 2.*r*r/w0/w0); 	
 	double res = w*w*amp*cos(arg); 
 	return res*Envelope(x,t);
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
 if(wave_type == 3){
 	double r = sqrt(y*y+z*z);
	double phi = atan2(z,y);
 	double arg = w*t-k*x-l*phi;
 	
 	double amp = Eo * exp(-r*r/w0/w0);
 	if(l!=0) amp*= pow(r*sqrt(2.)/w0 , abs(l));
 	if(l!=0 && p!=0) amp *= assoc_laguerre(abs(p), abs(l), 2.*r*r/w0/w0);

 	double res = 2.*r/w0/w0 * cos(arg)*sin(phi);
 	if(p!=0){
 		double amp3 = Eo * exp(-r*r/w0/w0) * pow(r*sqrt(2.)/w0 , abs(l)) * assoc_laguerre(abs(p)-1, abs(l)+1, 2.*r*r/w0/w0); //com Laguerre_(p-1)_(l+1)
 		amp += 2.*amp3;
 	}
 	res *= amp;
 	if(l!=0){
 		double amp2 = Eo * exp(-r*r/w0/w0) * pow(sqrt(2.)/w0 , abs(l)) * pow(r, abs(l)-1) * assoc_laguerre(abs(p), abs(l), 2.*r*r/w0/w0); //com r^(|l|-1)
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
 if(wave_type == 3){
 	double r = sqrt(y*y+z*z);
	double phi = atan2(z,y);
 	double arg = w*t-k*x-l*phi;
 	
 	double amp = Eo * exp(-r*r/w0/w0);
 	if(l!=0) amp*= pow(r*sqrt(2.)/w0 , abs(l));
 	if(l!=0 && p!=0) amp *= assoc_laguerre(abs(p), abs(l), 2.*r*r/w0/w0);

 	double res = -2.*r/w0/w0 * w*sin(arg)*sin(phi);
 	if(p!=0){
 		double amp3 = Eo * exp(-r*r/w0/w0) * pow(r*sqrt(2.)/w0 , abs(l)) * assoc_laguerre(abs(p)-1, abs(l)+1, 2.*r*r/w0/w0); //com Laguerre_(p-1)_(l+1)
 		amp += 2.*amp3;
 	}
 	res *= amp;
 	if(l!=0){
 		double amp2 = Eo * exp(-r*r/w0/w0) * pow(sqrt(2.)/w0 , abs(l)) * pow(r, abs(l)-1) * assoc_laguerre(abs(p), abs(l), 2.*r*r/w0/w0); //com r^(|l|-1)
 		res += amp2*w*(abs(l)*sin(arg)*sin(phi) -l*cos(arg)*cos(phi));
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
 	return kg*Eo / eta * w0/wz * exp(-r*r/(wz*wz)) * sin(w*t - kg*x - kg*r*r*R_1/2. + psi) * Envelope(x, t);
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
	double arg = w*t-k*x-l*phi;
	double amp=1;
	amp *= Eo;
	amp *= exp(-r*r/w0/w0);
	if(l!=0) amp *= pow(r*sqrt(2.)/w0 , abs(l));
	if (l!=0 && p!=0) amp *= assoc_laguerre(abs(p), abs(l), 2.*r*r/w0/w0); 	
 	double res = k*amp*sin(arg); 
 	return res*Envelope(x,t);
 } 

 else return 0.;
}

double DerBfz(double x, double y, double z, double p, double l, double t){ //Bz time derivative at (x,y,z)
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
	double arg = w*t-k*x-l*phi;
	double amp=1;
	amp *= Eo;
	amp *= exp(-r*r/w0/w0);
	if(l!=0) amp *= pow(r*sqrt(2.)/w0 , abs(l));
	if (l!=0 && p!=0) amp *= assoc_laguerre(abs(p), abs(l), 2.*r*r/w0/w0); 	
 	double res = k*w*amp*cos(arg); 
 	return res*Envelope(x,t);
 } 
 else return 0.;
}



double A3(double r, double phi, double x, double p, double l, double t){
	double arg = w*t-k*x-l*phi;
	return Eo*exp(-r*r/w0/w0) * pow(r*sqrt(2.)/w0, abs(l)) * assoc_laguerre(p,abs(l), 2*r*r/w0/w0) * cos(arg);
}
double Bfx3(double x, double y, double z, double p, double l, double t){
	double h=1e-5;
	double r = sqrt(y*y+z*z);
	double phi = atan2(z,y);
	double der1 = ((r+h/2)*A3(r+h/2,phi,x,p,l,t)-(r-h/2)*A3(r-h/2,phi,x,p,l,t))/h * (-sin(phi));
	double der2 = (A3(r,phi+h/2,x,p,l,t)*cos(phi+h/2)-A3(r,phi-h/2,x,p,l,t)*cos(phi-h/2))/h;
	return 1./r * (der1-der2);
}
double DerBfx3(double x, double y, double z, double p, double l, double t){
	double h=1e-5;
	return (Bfx3(x,y,z,p,l,t+h/2)-Bfx3(x,y,z,p,l,t-h/2))/h;
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
	wave_type=3;
	double a1 = Bfx(par[0],x[0],x[1],par[1],par[2],par[3]);
	double a2 = Bfx3(par[0],x[0],x[1],par[1],par[2],par[3]);
	return a1-a2;
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
	w=k;


	cout<<"WAVE_TYPE = "<<wave_type<<endl;
	cout<<"Ao = "<<Eo<<endl;
	cout<<"zr = "<<zr<<endl;
	cout<<"delta = "<<delta<<endl;
	cout<<"w,k = "<<w<<endl;



	double range=1;

	int n_l=6; n_l=1;
	int n_p=4; n_p=1;
	int side=500;

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


	double t=0; double dt=0.05; double tempo=1000;
	while(t<tempo){
	for(int p=0; p<n_p; p++){
		f1[p] = new TF2*[n_l];
		t1[p] = new TLatex*[n_l];
		for(int l=0; l<n_l; l++){
			f1[p][l] = new TF2("",teste,-20,20,-20,20,4);
			f1[p][l]->SetParameter(1,pr);
			f1[p][l]->SetParameter(2,lr);
			f1[p][l]->SetParameter(0,0); //x
			f1[p][l]->SetParameter(3,t); //t
			f1[p][l]->SetMinimum(-range); 
			f1[p][l]->SetMaximum(range); 
			f1[p][l]->SetNpx(500);
			f1[p][l]->SetNpy(500);
			cout<<f1[p][l]->GetMaximum()<<"\t"<<f1[p][l]->GetMinimum()<<endl;
			c1->cd(1+l+p*n_l);
			//c1->cd(1+l+p*n_l)->SetRightMargin(0.0);
			//c1->cd(1+l+p*n_l)->SetLeftMargin(0.0);
			//c1->cd(1+l+p*n_l)->SetBottomMargin(0.0);
			//c1->cd(1+l+p*n_l)->SetTopMargin(0.0);
			f1[p][l]->Draw("colz");
			t1[p][l] = new TLatex(-5,17,("#font[132]{p = "+to_string(p)+" | l = "+to_string(l)+"}").c_str());
			t1[p][l]->Draw("SAME");
		}
	}
	c1->Update();
	t+=dt;
	}


	c1->SaveAs("Plot.png");


	return 0;
}
