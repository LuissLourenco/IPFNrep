#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
#include<iostream>
#include<cmath>
#include<stdlib.h>
using namespace std;
double t, xgrid, ygrid,gam,Fx,Fy,Fz,pE,k,dx,dy,dz,w,Eo,Bo,delta,kdamp, tfwhm, stable;  
int pri, Ni, Nj;

int wave_type = 1; // 0-> Plane Wave; 1-> Gaussian Beam; 
					//2->Laguerre-Gaussian Beam; 
					//3->Laguerre-Gaussian Beam simplificado
double lambda, w0, n, eta, zr, kg;
int l,p;

bool run_kdamp = true;

double Poly(double x){
 	return 10*x*x*x-15*x*x*x*x+6*x*x*x*x*x;
}

double Envelope(double x,double t){
	double res,x1,x2; 
	x1=(x-t+2*tfwhm+stable)/tfwhm; 
	x2=(t-x)/tfwhm;
	if (x<t-2*tfwhm-stable) res=0;
	else if (x<t-tfwhm-stable) res=Poly(x1);
	else if (x<t-tfwhm) res = 1;
	else if (x<t) res=Poly(x2);
	else res=0; 
	return res;
}



/*
double Ax(double phi){return 0;}
double Ay(double phi){return delta*Eo*sin(phi);}
double Az(double phi){return sqrt(1-delta*delta)*Eo*cos(phi);}

double Phi(double x, double y, double z, double t){
		return w*t - k*x;
	};

double Der(double(*f)(double), double x){
	double h = 1e-5;
	return (f(x+h/2.)-f(x-h/2.))/h; 
};

double Der2(double(*f)(double), double x){
	double h = 1e-5;
	return (f(x+h)-2*f(x)+f(x-h))/h/h;
}

double Efx(double x, double y, double z){return 0;}
double Efy(double x, double y, double z){return -w*Der(Ay,Phi(x,y,z,t)) *Envelope(x,t);}
double Efz(double x, double y, double z){return -w*Der(Az,Phi(x,y,z,t)) *Envelope(x,t);}

double DerEfx(double x, double y, double z){return 0;}
double DerEfy(double x, double y, double z){return -w*w*Der2(Ay,Phi(x,y,z,t)) *Envelope(x,t);}
double DerEfz(double x, double y, double z){return -w*w*Der2(Az,Phi(x,y,z,t)) *Envelope(x,t);}

double Bfx(double x, double y, double z){return 0;}
double Bfy(double x, double y, double z){return k*Der(Az,Phi(x,y,z,t)) *Envelope(x,t);}
double Bfz(double x, double y, double z){return -k*Der(Ay,Phi(x,y,z,t)) *Envelope(x,t);}

double DerBfx(double x, double y, double z){return 0;}
double DerBfy(double x, double y, double z){return w*k*Der2(Az,Phi(x,y,z,t)) *Envelope(x,t);}
double DerBfz(double x, double y, double z){return -w*k*Der2(Ay,Phi(x,y,z,t)) *Envelope(x,t);}


*/





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
 	double res=1;

 	res *= w;
 	res *= Eo*w0/wz;
 	if(l!=0) res *= pow(r*sqrt(2.)/wz, abs(l));
 	if(l!=0 && p!=0) res *= assoc_laguerre(abs(p), abs(l), 2.*r*r/wz/wz);
 	res *= exp(-r*r/wz/wz); 	
 	res *= sin(arg);
 	return res*Envelope(x,t);
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
 if(wave_type == 2) return 0;
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
 if(wave_type == 2) return 0;
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
	if (l!=0 || p!=0) amp *= assoc_laguerre(abs(p), abs(l), 2.*r*r/w0/w0); 	
 	double res = k*w*amp*cos(arg); 
 	return res*Envelope(x,t);
 } 
 else return 0.;
}





double fun1( double px,double py,double pz,double x,double y, double z){   //  dpx/dt
	if(run_kdamp){
		double res,A,B,C,D;
		//A=gam*DerEfx(x,y,z)+px*(Efx(x+dx,y,z)-Efx(x-dx,y,z))/(2*dx)+py*(Efx(x,y+dy,z)-Efx(x,y-dy,z))/(2*dy)+py*DerBfz(x,y,z)-pz*DerBfy(x,y,z);
		//B=py*px*(Bfz(x+dx,y,z)-Bfz(x-dx,y,z))/(2*dx*gam)+py*py*(Bfz(x,y+dy,z)-Bfz(x,y-dy,z))/(2*dy*gam)-pz*px*(Bfy(x+dx,y,z)-Bfy(x-dx,y,z))/(2*dx*gam)-pz*py*(Bfy(x,y+dy,z)-Bfy(x,y-dy,z))/(2*dy*gam);
		
		A  = gam*DerEfx(x,y,z);
		A += px*(Efx(x+dx,y,z)-Efx(x-dx,y,z))/(2*dx);
		A += py*(Efx(x,y+dy,z)-Efx(x,y-dy,z))/(2*dy);
		A += pz*(Efx(x,y,z+dz)-Efx(x,y,z-dz))/(2*dz);
		A += py*DerBfz(x,y,z)-pz*DerBfy(x,y,z);

		B  = py*px*(Bfz(x+dx,y,z)-Bfz(x-dx,y,z))/(2*dx*gam);
		B += py*py*(Bfz(x,y+dy,z)-Bfz(x,y-dy,z))/(2*dy*gam);
		B += py*pz*(Bfz(x,y,z+dz)-Bfz(x,y,z-dz))/(2*dz*gam);
		B -= pz*px*(Bfy(x+dx,y,z)-Bfy(x-dx,y,z))/(2*dx*gam);
		B -= pz*py*(Bfy(x,y+dy,z)-Bfy(x,y-dy,z))/(2*dy*gam);
		B -= pz*pz*(Bfy(x,y,z+dz)-Bfy(x,y,z-dz))/(2*dz*gam);

		C=Efy(x,y,z)*Bfz(x,y,z)-Efz(x,y,z)*Bfy(x,y,z)+(Bfy(x,y,z)*Bfx(x,y,z)*py-Bfy(x,y,z)*Bfy(x,y,z)*px-Bfz(x,y,z)*Bfz(x,y,z)*px+Bfz(x,y,z)*Bfx(x,y,z)*pz)/gam;
		D=Efx(x,y,z)*pE/gam+px*pE*pE/gam-gam*px*(Fx*Fx+Fy*Fy+Fz*Fz);
		res=Fx+kdamp*(A+B+C+D);
		return res;
	}else return Fx;
}

double fun2( double px,double py,double pz,double x,double y, double z){   //  dpy/dt
	if(run_kdamp){
		double res,A,B,C,D;
		//A=gam*DerEfy(x,y,z)+px*(Efy(x+dx,y,z)-Efy(x-dx,y,z))/(2*dx)+py*(Efy(x,y+dy,z)-Efy(x,y-dy,z))/(2*dy)+pz*DerBfx(x,y,z)-px*DerBfz(x,y,z);
		//B=pz*px*(Bfx(x+dx,y,z)-Bfx(x-dx,y,z))/(2*dx*gam)+pz*py*(Bfx(x,y+dy,z)-Bfx(x,y-dy,z))/(2*dy*gam)-px*px*(Bfz(x+dx,y,z)-Bfz(x-dx,y,z))/(2*dx*gam)-px*py*(Bfz(x,y+dy,z)-Bfz(x,y-dy,z))/(2*dy*gam);
		
		A  = gam*DerEfy(x,y,z);
		A += px*(Efy(x+dx,y,z)-Efy(x-dx,y,z))/(2*dx);
		A += py*(Efy(x,y+dy,z)-Efy(x,y-dy,z))/(2*dy);
		A += pz*(Efy(x,y,z+dz)-Efy(x,y,z-dz))/(2*dz);
		A += pz*DerBfx(x,y,z)-px*DerBfz(x,y,z);

		B  = pz*px*(Bfx(x+dx,y,z)-Bfx(x-dx,y,z))/(2*dx*gam);
		B += pz*py*(Bfx(x,y+dy,z)-Bfx(x,y-dy,z))/(2*dy*gam);
		B += pz*pz*(Bfx(x,y,z+dz)-Bfx(x,y,z-dz))/(2*dz*gam);
		B -= px*px*(Bfz(x+dx,y,z)-Bfz(x-dx,y,z))/(2*dx*gam);
		B -= px*py*(Bfz(x,y+dy,z)-Bfz(x,y-dy,z))/(2*dy*gam);
		B -= px*pz*(Bfz(x,y,z+dz)-Bfz(x,y,z-dz))/(2*dz*gam);

		C=Efz(x,y,z)*Bfx(x,y,z)-Efx(x,y,z)*Bfz(x,y,z)+(Bfz(x,y,z)*Bfy(x,y,z)*pz-Bfz(x,y,z)*Bfz(x,y,z)*py-Bfx(x,y,z)*Bfx(x,y,z)*py+Bfx(x,y,z)*Bfy(x,y,z)*px)/gam;
		D=Efy(x,y,z)*pE/gam+py*pE*pE/gam-gam*py*(Fx*Fx+Fy*Fy+Fz*Fz);
		res=Fy+kdamp*(A+B+C+D);
		return res;
	}else return Fy;
}

double fun3( double px,double py,double pz,double x,double y, double z){   //  dpz/dt
	if(run_kdamp){
		double res,A,B,C,D;
		//A=gam*DerEfz(x,y,z)+px*(Efz(x+dx,y,z)-Efz(x-dx,y,z))/(2*dx)+px*DerBfy(x,y,z)-py*DerBfx(x,y,z);
		//B=px*px*(Bfy(x+dx,y,z)-Bfy(x-dx,y,z))/(2*dx*gam)+px*py*(Bfy(x,y+dy,z)-Bfy(x,y-dy,z))/(2*dy*gam)-py*px*(Bfx(x+dx,y,z)-Bfx(x-dx,y,z))/(2*dx*gam)-py*py*(Bfx(x,y+dy,z)-Bfx(x,y-dy,z))/(2*dy*gam);

		A  = gam*DerEfz(x,y,z);
		A += px*(Efz(x+dx,y,z)-Efz(x-dx,y,z))/(2*dx);
		A += py*(Efz(x,y+dy,z)-Efz(x,y-dy,z))/(2*dy);
		A += pz*(Efz(x,y,z+dz)-Efz(x,y,z-dz))/(2*dz);
		A += px*DerBfy(x,y,z)-py*DerBfx(x,y,z);

		B  = px*px*(Bfy(x+dx,y,z)-Bfy(x-dx,y,z))/(2*dx*gam);
		B += px*py*(Bfy(x,y+dy,z)-Bfy(x,y-dy,z))/(2*dy*gam);
		B += px*pz*(Bfy(x,y,z+dz)-Bfy(x,y,z-dz))/(2*dz*gam);
		B -= py*px*(Bfx(x+dx,y,z)-Bfx(x-dx,y,z))/(2*dx*gam);
		B -= py*py*(Bfx(x,y+dy,z)-Bfx(x,y-dy,z))/(2*dy*gam);
		B -= py*pz*(Bfx(x,y,z+dz)-Bfx(x,y,z-dz))/(2*dz*gam);

		C=Efx(x,y,z)*Bfy(x,y,z)-Efy(x,y,z)*Bfx(x,y,z)+(Bfx(x,y,z)*Bfz(x,y,z)*px-Bfx(x,y,z)*Bfx(x,y,z)*pz-Bfy(x,y,z)*Bfy(x,y,z)*pz+Bfy(x,y,z)*Bfz(x,y,z)*py)/gam;
		D=Efz(x,y,z)*pE/gam+pz*pE*pE/gam-gam*pz*(Fx*Fx+Fy*Fy+Fz*Fz);
		res=Fz+kdamp*(A+B+C+D);
		return res;
	}else return Fz;
}


double fun4(double px,double py, double pz, double x, double y, double z)   // dx/dt
{double res;
 res=px/sqrt(1+px*px+py*py+pz*pz);
 return res;
}

double fun5(double px,double py, double pz, double x, double y, double z)    //  dy/dt
{double res;
 res=py/sqrt(1+px*px+py*py+pz*pz);
 return res;
}

double fun6(double px,double py, double pz, double x, double y, double z)    //  dz/dt
{double res;
 res=pz/sqrt(1+px*px+py*py+pz*pz);
 return res;
}


//procedure RK does runge-kutta integration of the model defined with functions above (fun1,2,3,4,5,6)

double RK(double T, long long int N, double p01, double p02,double p03, double x01, double x02, double x03) 
{double k11, k12, k13, k14,k15, k16, k21, k22, k23, k24, k25, k26, k31, k32, k33, k34, k35, k36, k41, k42, k43, k44, k45, k46;
 double w1, w2, w3, w4, w5, w6, h;
 double px, py, pz, x, y, z;
 long long int i; 
 FILE *fo;
 
 
 fo=fopen("Out3.txt","w");
 //fprintf(fo,"t	 x		 y      z	   px		 py		  pz	 gamma\n");
 h=T/N; 
 
 w1=p01; w2=p02; w3=p03; w4=x01; w5=x02; w6=x03;  gam=sqrt(1+p01*p01+p02*p02+p03*p03); t=0;
 

 fprintf(fo, "%.14e    	%.14e 	 %.14e  	%.14e    %.14e	 %.14e	%.14e		%.14e\n",t, w4, w5, w6, w1, w2, w3, gam); 
 for(i=1; i<N; i++)
 	{
 	 {
 	 t=(i-1)*h;
 	
 	 px=w1; py=w2; pz=w3; x=w4; y=w5; z=w6;
 	 gam=sqrt(1+px*px+py*py+pz*pz); pE=px*Efx(x,y,z)+py*Efy(x,y,z)+pz*Efz(x,y,z); 
 	 Fx=Efx(x,y,z)+(py*Bfz(x,y,z)-pz*Bfy(x,y,z))/gam; 	 
 	 Fy=Efy(x,y,z)+(pz*Bfx(x,y,z)-px*Bfz(x,y,z))/gam; 
 	 Fz=Efz(x,y,z)+(px*Bfy(x,y,z)-py*Bfx(x,y,z))/gam; 
 	 
 	 k11=h*fun1(px,py,pz,x,y,z); k12=h*fun2(px,py,pz,x,y,z); k13=h*fun3(px,py,pz,x,y,z); 
 	 k14=h*fun4(px,py,pz,x,y,z); k15=h*fun5(px,py,pz,x,y,z); k16=h*fun6(px,py,pz,x,y,z);
 	 
 	//------------------------------------------------------------------------
 	 
 	 t=(i-1)*h+h/2;
 
 	 
 	 px=w1+k11*1/2; py=w2+k12*1/2; pz=w3+k13*1/2; x=w4+k14*1/2; y=w5+k15*1/2; z=w6+k16*1/2;
 	 
 	 gam=sqrt(1+px*px+py*py+pz*pz); pE=px*Efx(x,y,z)+py*Efy(x,y,z)+pz*Efz(x,y,z); 
 	 Fx=Efx(x,y,z)+(py*Bfz(x,y,z)-pz*Bfy(x,y,z))/gam; 	 
 	 Fy=Efy(x,y,z)+(pz*Bfx(x,y,z)-px*Bfz(x,y,z))/gam; 
 	 Fz=Efz(x,y,z)+(px*Bfy(x,y,z)-py*Bfx(x,y,z))/gam; 
 	  
	 k21=h*fun1(px,py,pz,x,y,z); 
 	 k22=h*fun2(px,py,pz,x,y,z);
 	 k23=h*fun3(px,py,pz,x,y,z);
 	 k24=h*fun4(px,py,pz,x,y,z);
 	 k25=h*fun5(px,py,pz,x,y,z);
 	 k26=h*fun6(px,py,pz,x,y,z);

 	 //-----------------------------------------------------------------------
 	 
 	 px=w1+k21*1/2; py=w2+k22*1/2; pz=w3+k23*1/2; x=w4+k24*1/2; y=w5+k25*1/2; z=w6+k26*1/2;
 	
 	 gam=sqrt(1+px*px+py*py+pz*pz); pE=px*Efx(x,y,z)+py*Efy(x,y,z)+pz*Efz(x,y,z); 
 	 Fx=Efx(x,y,z)+(py*Bfz(x,y,z)-pz*Bfy(x,y,z))/gam; 	 
 	 Fy=Efy(x,y,z)+(pz*Bfx(x,y,z)-px*Bfz(x,y,z))/gam; 
 	 Fz=Efz(x,y,z)+(px*Bfy(x,y,z)-py*Bfx(x,y,z))/gam; 
 	 
 	 k31=h*fun1(px,py,pz,x,y,z); 
 	 k32=h*fun2(px,py,pz,x,y,z);
 	 k33=h*fun3(px,py,pz,x,y,z);
 	 k34=h*fun4(px,py,pz,x,y,z);
 	 k35=h*fun5(px,py,pz,x,y,z);
 	 k36=h*fun6(px,py,pz,x,y,z);

 	 //---------------------------------------------------------------------------
 	  
 	 t=i*h;
 	 
 	 px=w1+k31; py=w2+k32; pz=w3+k33; x=w4+k34; y=w5+k35; z=w6+k36;
 	 
 	 gam=sqrt(1+px*px+py*py+pz*pz); pE=px*Efx(x,y,z)+py*Efy(x,y,z)+pz*Efz(x,y,z); 
 	 Fx=Efx(x,y,z)+(py*Bfz(x,y,z)-pz*Bfy(x,y,z))/gam; 	 
 	 Fy=Efy(x,y,z)+(pz*Bfx(x,y,z)-px*Bfz(x,y,z))/gam; 
 	 Fz=Efz(x,y,z)+(px*Bfy(x,y,z)-py*Bfx(x,y,z))/gam; 
 	 
 	 k41=h*fun1(px,py,pz,x,y,z);  
 	 k42=h*fun2(px,py,pz,x,y,z);
 	 k43=h*fun3(px,py,pz,x,y,z);
 	 k44=h*fun4(px,py,pz,x,y,z);
 	 k45=h*fun5(px,py,pz,x,y,z);
 	 k46=h*fun6(px,py,pz,x,y,z);

 	 
 	 w1=w1+(k11+2*k21+2*k31+k41)/6; 
 	 w2=w2+(k12+2*k22+2*k32+k42)/6; 
 	 w3=w3+(k13+2*k23+2*k33+k43)/6; 
 	 w4=w4+(k14+2*k24+2*k34+k44)/6; 
 	 w5=w5+(k15+2*k25+2*k35+k45)/6;
 	 w6=w6+(k16+2*k26+2*k36+k46)/6;
 	 
 
 	 //---------------------------------------------------------------------------
 	 
 	 if ((i % pri) == 0)
 	 {gam=sqrt(1+w1*w1+w2*w2+w3*w3); 
 	  fprintf(fo, "%.14e    	%.14e 	 %.14e   %.14e 	%.14e	 %.14e	%.14e		%.14e\n",t, w4, w5, w6, w1, w2, w3, gam); 
 	  // prints t, x, y, z, px, py, pz, gamma
	  }
	 }
 	}
 	fclose(fo);
 	return 0;
}


int main(){ 

	double p01, p02, p03, x01, x02, x03;
  	double T; 
	long long int N; 
	FILE *foo;

	char trash[128];
 
	foo=fopen("InputToBatch.txt","r");
	fscanf(foo,"%s %lf %lf %lf %lf %lf %lf %lf %lf %lli %i %lf %i %lf %lf %lf %lf %lf %lf %lf %lf %i %i", 
				trash, &x01, &x02, &x03, &p01, &p02, &p03, &kdamp, &T, &N, &pri, &dx, 
				&wave_type, &tfwhm, &stable, &Eo, &delta, &w0, &lambda, &n, &eta, &l, &p);
	fclose(foo);

	//trash|x01|x02|x03|p01|p02|p03|kdamp|T|N|pri|dx|wave_type|tfwhm|stable|E0|delta|w0|lambda|n|eta|l|p

	dy = dx; dz = dx;
 	k = 1;
 	Bo = Eo; 
 	if(wave_type>=1){
 		kg = 2. * M_PI * n / lambda;
 		k = kg;
 	}
 	w=k;
 	zr = M_PI*w0*w0*n/lambda;

 	RK(T, N, p01, p02, p03, x01, x02, x03);

 	printf("gotovo\n");

 	return 1;

 }
 

 /*
 gcc -O2 trial.c -o trial
./trial  
*/