#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1

#include <math.h>
#include <iostream>

#include "TF2.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TEllipse.h"

using namespace std;

double Eo = 1;
double n = 1;
double w = 1;

double lambda = 1;
double w0 = 5;

double laguerre_mode(double* coord, double* par){
	double y = coord[0];
	double z = coord[1];
	double x = par[0];
	double t = par[1];

	int p = par[2];
	int l = par[3];

	double r = sqrt(y*y + z*z);
	double phi = atan2(z, y);

 	double zr = M_PI * w0*w0 * n / lambda;
 	double wz = w0 * sqrt(1 + x*x/(zr*zr));
 	double R_1 = x / (x*x + zr*zr);
 	double kg = 2 * M_PI * n / lambda;
 	double psi = atan(x/zr);
 
	// RETURNS
	double res = cos(w*t - kg*x - kg*r*r*R_1/2 - l*phi + (2*p+abs(l)+1)*psi);
	res *= Eo * w0/wz * exp(-r*r/(wz*wz));
	res *= pow(sqrt(2)*r/wz, abs(l));
	res *= assoc_laguerre(p, abs(l), 2*r*r/(wz*wz));

	return res;
}


double magnetic_field_rafa(double* coord, double* par){

	double h = 0.00001;

	double par1[4] = {par[0]  , par[1]+M_PI/2/w, par[2], par[3]};
	double par2[4] = {par[0]+h, par[1]+M_PI/2/w, par[2], par[3]};

	return (laguerre_mode(coord, par2) - laguerre_mode(coord, par1)) / h;

}


int main(){
	/*
	TF2* f1 = new TF2("f1", laguerre_mode, -20, 20, -20, 20, 4); 
	f1->SetTitle(";z;r;Gaussian Beam");
	///gStyle->SetPalette(kMint);
	f1->SetNpx(500);
	f1->SetNpy(500);

	f1->SetParameter(0, 0); // SET X COORDINATE - PROPAGATION 
	f1->SetParameter(1, 0); // SET TIME
	f1->SetParameter(2, 10); // SET "p" VARIABLE
	f1->SetParameter(3, 10); // SET "l" VARIABLE

	TCanvas* c1 = new TCanvas("c", "", 1000, 1000);

	c1->cd();

	f1->Draw("COL"); 
	c1->SaveAs("LaguerreGaussianWaveModes.png");
	*/


	int n_l = 6;
	int n_p = 4;
	int side = 400;

	TCanvas* c1 = new TCanvas("c", "", n_l*side, n_p*side);
	c1->Divide(n_l, n_p);
	c1->SetRightMargin(0.0);
	c1->SetLeftMargin(0.0);
	c1->SetBottomMargin(0.0);
	c1->SetTopMargin(0.0);

	gStyle->SetPalette(kBird);

	TF2*** fun = new TF2**[n_p];
	TLatex*** text = new TLatex**[n_p];
	TEllipse*** ellipse = new TEllipse**[n_p];
	for(int p = 0; p < n_p; p++){
		fun[p] = new TF2*[n_l];
		text[p] = new TLatex*[n_l];
		ellipse[p] = new TEllipse*[n_l];
		for(int l = 0; l < n_l; l++){
			fun[p][l] = new TF2(("f_"+to_string(p)+"_"+to_string(l)).c_str(), magnetic_field_rafa, -20, 20, -20, 20, 4); 
			fun[p][l]->SetTitle(";;;");
			fun[p][l]->SetNpx(500);
			fun[p][l]->SetNpy(500);

			fun[p][l]->SetParameter(0, 0); // SET X COORDINATE - PROPAGATION 
			fun[p][l]->SetParameter(1, 0); // SET TIME
			fun[p][l]->SetParameter(2, p); // SET "p" VARIABLE
			fun[p][l]->SetParameter(3, l); // SET "l" VARIABLE

			double range = fun[p][l]->GetMaximum();

			//fun[p][l]->SetMaximum(range);
			//fun[p][l]->SetMinimum(-range);	

			c1->cd(1+l+p*n_l);

			c1->cd(1+l+p*n_l)->SetRightMargin(0.0);
			c1->cd(1+l+p*n_l)->SetLeftMargin(0.0);
			c1->cd(1+l+p*n_l)->SetBottomMargin(0.0);
			c1->cd(1+l+p*n_l)->SetTopMargin(0.0);

			fun[p][l]->Draw("COL");	

			text[p][l] = new TLatex(-13, 15.5, ("#font[132]{p = "+to_string(p)+" | l = "+to_string(l)+"}").c_str());
			text[p][l]->SetTextSize(0.15);
			text[p][l]->Draw("SAME");

			ellipse[p][l] = new TEllipse(0, 0, w0);
			ellipse[p][l]->SetFillStyle(0);
			ellipse[p][l]->SetLineColor(2);
			ellipse[p][l]->SetLineWidth(2);
			ellipse[p][l]->Draw();
		}
	}

	c1->SaveAs("LaguerreGaussianWaveModes.png");



	return 0;

}