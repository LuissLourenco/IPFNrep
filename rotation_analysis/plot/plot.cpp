#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1

#include <math.h>
#include <iostream>

#include "../src/DataAnalysis.cpp"

#include "TF2.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TEllipse.h"
#include "TFrame.h"
#include "TPad.h"

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


int mood(string file, int p, int l, string bottom, string opt=""){

	int n_points, n_cols;
	double** values = ReadFile((file + ".txt").c_str(), &n_cols, &n_points, true);
	//y	z	phi	r	period	thetamin	thetamax
	DataSet Y(n_points, values[0]);
	DataSet Z(n_points, values[1]);
	DataSet PERIOD(n_points, values[4]);
	DataSet AMPLITUDE = abs(DataSet(n_points, values[5]) - DataSet(n_points, values[6]));

	if(opt == "quadrant"){
		Y = Y.concat(Y);
		Z = Z.concat(Z*Var(-1.));
		PERIOD = PERIOD.concat(PERIOD);
		AMPLITUDE = AMPLITUDE.concat(AMPLITUDE);
		Y = Y.concat(Y*Var(-1.));
		Z = Z.concat(Z);
		PERIOD = PERIOD.concat(PERIOD);
		AMPLITUDE = AMPLITUDE.concat(AMPLITUDE);
	}

	TCanvas* c1 = new TCanvas("c", "", 1500, 550);
	c1->SetRightMargin(0.0);
	c1->SetLeftMargin(0.0);
	c1->SetBottomMargin(0.0);
	c1->SetTopMargin(0.0);

	gStyle->SetPalette(kBird);
	gStyle->SetNumberContours(500);

	TPad* pad1 = new TPad("", "", 0, 1./11, 1./3, 1);
	TPad* pad2 = new TPad("", "", 1./3, 1./11, 2./3, 1);
	TPad* pad3 = new TPad("", "", 2./3, 1./11, 1, 1);
	TPad* pad4 = new TPad("", "", 0, 0, 1, 1./11);

	//=====================================

	TF2* fun = new TF2(("f_"+to_string(p)+"_"+to_string(l)).c_str(), laguerre_mode, -20, 20, -20, 20, 4); 
	fun->SetTitle(";;;");
	fun->SetNpx(500);
	fun->SetNpy(500);
	fun->SetParameter(0, 0); // SET X COORDINATE - PROPAGATION 
	fun->SetParameter(1, 0); // SET TIME
	fun->SetParameter(2, p); // SET "p" VARIABLE
	fun->SetParameter(3, l); // SET "l" VARIABLE

	pad1->cd();

	pad1->SetRightMargin( 0.14);
	pad1->SetLeftMargin(  0.06);
	pad1->SetBottomMargin(0.08);
	pad1->SetTopMargin(   0.08);

	TH1 *frame_wave = c1->cd(1)->DrawFrame(-3.9, -3.9, 3.9, 3.9);
	frame_wave->SetTitle(("Wave Shape ( p = " + to_string(p) + " | l = " + to_string(l) + " );y_{0};z_{0}").c_str());
	frame_wave->GetXaxis()->SetLabelSize(0.03);
	frame_wave->GetXaxis()->SetTitleSize(0.06);
	frame_wave->GetXaxis()->SetTitleOffset(0.5);
	frame_wave->GetYaxis()->SetLabelSize(0.03);
	frame_wave->GetYaxis()->SetTitleSize(0.06);
	frame_wave->GetYaxis()->SetTitleOffset(0.3);
	fun->Draw("same COL");

	pad1->SetGrid(1, 1);
	pad1->RedrawAxis("g");

	//=====================================

	TGraph2D* period = GetTGraph2D(Y, Z, PERIOD);
	period->SetTitle("Period od Oscilation;y_{0};z_{0}");
	period->SetNpx(500);
	period->SetNpy(500);
	
	pad2->cd();

	pad2->SetRightMargin( 0.14);
	pad2->SetLeftMargin(  0.06);
	pad2->SetBottomMargin(0.08);
	pad2->SetTopMargin(   0.08);

	TH1 *frame_period = c1->cd(2)->DrawFrame(-3.9, -3.9, 3.9, 3.9);
	frame_period->SetTitle("Period of Oscilation;y_{0};z_{0}");
	frame_period->GetXaxis()->SetLabelSize(0.03);
	frame_period->GetXaxis()->SetTitleSize(0.06);
	frame_period->GetXaxis()->SetTitleOffset(0.5);
	frame_period->GetYaxis()->SetLabelSize(0.03);
	frame_period->GetYaxis()->SetTitleSize(0.06);
	frame_period->GetYaxis()->SetTitleOffset(0.3);
	frame_period->GetZaxis()->SetMaxDigits(2);
	period->Draw("same COLZ");

	pad2->SetGrid(1, 1);
	pad2->RedrawAxis("g");

	//=====================================

	TGraph2D* amplitude = GetTGraph2D(Y, Z, AMPLITUDE);
	amplitude->SetTitle("p_{y}^{max};y_{0};z_{0}");
	amplitude->SetNpx(500);
	amplitude->SetNpy(500);

	pad3->cd();

	pad3->SetRightMargin( 0.14);
	pad3->SetLeftMargin(  0.06);
	pad3->SetBottomMargin(0.08);
	pad3->SetTopMargin(   0.08);

	TH1 *frame_amplitude = c1->cd(3)->DrawFrame(-3.9, -3.9, 3.9, 3.9);
	frame_amplitude->SetTitle("Amplitude of Oscilation, #sigma;z_{0}");
	frame_amplitude->GetXaxis()->SetLabelSize(0.03);
	frame_amplitude->GetXaxis()->SetTitleSize(0.06);
	frame_amplitude->GetXaxis()->SetTitleOffset(0.5);
	frame_amplitude->GetYaxis()->SetLabelSize(0.03);
	frame_amplitude->GetYaxis()->SetTitleSize(0.06);
	frame_amplitude->GetYaxis()->SetTitleOffset(0.3);
	frame_amplitude->GetZaxis()->SetMaxDigits(2);

	frame_amplitude->GetZaxis()->SetBinLabel(1,"-#pi");

	amplitude->Draw("same COLZ");

	pad3->SetGrid(1, 1);
	pad3->RedrawAxis("g");


	//=====================================


	pad4->cd();

	TLatex* bottom_string = new TLatex(0.01, 0.4, bottom.c_str());
	bottom_string->SetTextSize(0.4);
	bottom_string->Draw();


	//=====================================


	c1->cd(); 
	pad1->Draw();
	pad2->Draw();
	pad3->Draw();
	pad4->Draw();

	c1->SaveAs((file + ".png").c_str());


	delete c1;
	delete fun;
	delete period;
	delete amplitude;

	return 0;

}

int main(){

	/*
	mood(x1, x2, x3, x4, x5)
	x1 -> NAME OF THE TXT FILE WITHOUT THE EXTENSION
	x2 -> p OF THE WAVE
	x3 -> l OF THE WAVE
	x4 -> TEXT ON THE BOTTOM OF THE PLOTS
	x5 -> (OPTIONAL) IF THE TXT FILE HAS ONLY THE FIRST QUADRANT, USE "quadrant"
		TO REPLICATE TO OTHERS
	*/

	/*
	mood("../final_outputs/p0l1_px-10_a30_dt0.005" , 0, 1, "Laguerre Gaussian Mode ( p = 0 | l = 1 ) with Longitudinal B and E fields; Initial Longitudinal Momentum = -10; w0 = 5; a0 = 30; lambda = 1; kdamp = 0; Linear Polarization on y", "quadrant");
	mood("../final_outputs/p0l1_px-10_a50_dt0.003" , 0, 1, "Laguerre Gaussian Mode ( p = 0 | l = 1 ) with Longitudinal B and E fields; Initial Longitudinal Momentum = -10; w0 = 5; a0 = 50; lambda = 1; kdamp = 0; Linear Polarization on y", "quadrant");
	mood("../final_outputs/p0l1_px-10_a70_dt0.003" , 0, 1, "Laguerre Gaussian Mode ( p = 0 | l = 1 ) with Longitudinal B and E fields; Initial Longitudinal Momentum = -10; w0 = 5; a0 = 70; lambda = 1; kdamp = 0; Linear Polarization on y", "quadrant");
	mood("../final_outputs/p0l1_px-10_a100_dt0.001", 0, 1, "Laguerre Gaussian Mode ( p = 0 | l = 1 ) with Longitudinal B and E fields; Initial Longitudinal Momentum = -10; w0 = 5; a0 = 100; lambda = 1; kdamp = 0; Linear Polarization on y", "quadrant");
	*/	
	mood("../final_outputs/Data_Out", 0, 1, "pila");

	return 0;
}

