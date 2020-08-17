
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
double w0=5;
double lambda=1;
double w =1;
complex<double> Eo=imu;
double kg = 2*M_PI*n/lambda;


complex<double> u_aux(double x, double y, double z, double l, double p){

	complex<double> res(1,0);

	double r = sqrt(y*y+z*z);
	double phi = atan2(z,y);
	double kg = 2.*M_PI*n/lambda;
	double zr = M_PI*w0*w0*n/lambda;
	double wz = w0*sqrt(1+x*x/zr/zr);

	res *= Eo*w0/wz;
	res *= pow(r*sqrt(2.)/wz, abs(l));
	res *= assoc_laguerre(abs(p), abs(l), 2*r*r/wz/wz);
	res *= exp(-r*r/wz/wz);
	res *= exp(-imu*kg*r*r*x/2./(x*x+zr*zr));
	res *= exp(-imu*l*phi);
	res *= exp(imu*(2.*p+abs(l)+1.)*atan(x/zr));

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
	double zr = M_PI*w0*w0*n/lambda;
	double wz = w0*sqrt(1+x*x/zr/zr);

	complex<double> mood, mood1(0,0), mood2(0,0);

	mood1 = -x/(x*x+zr*zr)+2*r*r*x*zr*zr/w0/w0/(x*x+zr*zr)/(x*x+zr*zr)-imu*kg*r*r*(0.5-x*x/(x*x+zr*zr))/(x*x+zr*zr)+imu*(2.*p+abs(l)+1.)*zr/(x*x+zr*zr)-imu*kg;
	if((int)l!=0) mood1 += -fabs(l)*x/(x*x+zr*zr);
	mood1 *= u_aux(x,y,z,l,p);
	if (p!=0) mood2 = u_aux(x,y,z,l+1,p-1)*4.*r*r*x*zr*zr/w0/w0/(x*x+zr*zr)/(x*x+zr*zr);
	mood = mood1+mood2;

	res=mood*exp(imu*(w*t-kg*x));

	return res;
}

complex<double> electric_field(double x, double y, double z, double l, double p, double t){
	complex<double> res;
	double r = sqrt(y*y+z*z);
	double phi = atan2(z,y);
	res = -imu*w*u_aux(x,y,z,l,p)*exp(imu*(w*t-kg*x));
	return res;
}

complex<double> magnetic_field(double x, double y, double z, double l, double p, double t){
	complex<double> res;
	double h=1e-5;
	res = (u_aux(x+h/2.,y,z,l,p)-u_aux(x-h/2.,y,z,l,p))/h - imu*kg*u_aux(x,y,z,l,p);
	res *= exp(imu*(w*t-kg*x));
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

double teste(double*x,double*par){
	return funB(x,par)-funBluis(x,par);
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
			f1[p][l]->SetParameter(1,l);
			f1[p][l]->SetParameter(2,p);
			f1[p][l]->SetParameter(0,0); //x
			f1[p][l]->SetParameter(3,0); //t
			f1[p][l]->SetNpx(500);
			f1[p][l]->SetNpy(500);
			cout<<f1[p][l]->GetMaximum()<<endl;
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
