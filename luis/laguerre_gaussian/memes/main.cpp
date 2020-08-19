
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

using namespace std;



complex<double> imu(0,1);

double n=1;
double eta=1;
double w0=10;
double lambda=1;
double w =1;
double Eo=1;
double kg = 2*M_PI*n/lambda;
double k = kg;
double zr = M_PI*w0*w0*n/lambda;

double der_coef4[5] = {1./12., -2./3., 0., 2./3., -1./12.};






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
	return funB(x,par)-funBtrig(x,par);
}


int main(int argc, char** argv){

	int n_l=6;
	int n_p=4;
	int side=500;

	auto f1 = new TF2**[n_p];
	auto t1 = new TLatex**[n_p];

	auto c1 = new TCanvas("c1", "", n_l*side, n_p*side);
	c1->Divide(n_l,n_p);
	c1->SetRightMargin(0.0);
	c1->SetLeftMargin(0.0);
	c1->SetBottomMargin(0.0);
	c1->SetTopMargin(0.0);

	gStyle->SetPalette(kBird);

	for(int p=0; p<n_p; p++){
		f1[p] = new TF2*[n_l];
		t1[p] = new TLatex*[n_l];
		for(int l=0; l<n_l; l++){
			f1[p][l] = new TF2("",teste,-20,20,-20,20,4);
			f1[p][l]->SetParameter(1,p);
			f1[p][l]->SetParameter(2,l);
			f1[p][l]->SetParameter(0,0); //x
			f1[p][l]->SetParameter(3,M_PI/2./w); //t
			f1[p][l]->SetNpx(500);
			f1[p][l]->SetNpy(500);
			cout<<f1[p][l]->GetMaximum()<<"\t"<<f1[p][l]->GetMinimum()<<endl;
			c1->cd(1+l+p*n_l);
			c1->cd(1+l+p*n_l)->SetRightMargin(0.0);
			c1->cd(1+l+p*n_l)->SetLeftMargin(0.0);
			c1->cd(1+l+p*n_l)->SetBottomMargin(0.0);
			c1->cd(1+l+p*n_l)->SetTopMargin(0.0);
			f1[p][l]->Draw("colz");
			t1[p][l] = new TLatex(-5,17,("#font[132]{p = "+to_string(p)+" | l = "+to_string(l)+"}").c_str());
			t1[p][l]->Draw("SAME");
		}
	}
	


	c1->SaveAs("Plot.png");


	return 0;
}
