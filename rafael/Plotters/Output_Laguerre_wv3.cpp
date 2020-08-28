#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1

#include <math.h>
#include <iostream>

#include "DataAnalysis.cpp"

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


int mood(string file, int p, int l, string bottom){

	int n_points, n_cols;
	double** values = ReadFile((file + ".txt").c_str(), &n_cols, &n_points, true);
	//y	z	x_f	y_f	z_f	px_f	py_f	pz_f	p_y_max	p_z_max	E_f	L_x_max	L_x_f
	DataSet Y(n_points, values[0]);
	DataSet Z(n_points, values[1]);
	
	DataSet YF(n_points, values[3]);
	DataSet ZF(n_points, values[4]);

	DataSet DY = sqrt(YF*YF+ZF*ZF) - sqrt(Y*Y+Z*Z);
	DataSet PYM = atan2(YF, ZF) - atan2(Y, Z);

	TCanvas* c1 = new TCanvas("c", "", 1500, 1050);
	//c1->Divide(3, 1);
	c1->SetRightMargin(0.0);
	c1->SetLeftMargin(0.0);
	c1->SetBottomMargin(0.0);
	c1->SetTopMargin(0.0);

	gStyle->SetPalette(kBird);

	TPad* pad1 = new TPad("", "", 0, 11./21, 1./3, 1);
	TPad* pad2 = new TPad("", "", 1./3, 11./21, 2./3, 1);
	TPad* pad3 = new TPad("", "", 2./3, 11./21, 1, 1);
	TPad* pad4 = new TPad("", "", 0, 1./21, 1./3, 11./21);
	TPad* pad5 = new TPad("", "", 1./3, 1./21, 2./3, 11./21);
	TPad* pad6 = new TPad("", "", 2./3, 1./21, 1, 11./21);
	TPad* pad7 = new TPad("", "", 0, 0, 1, 1./21);

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

	TH1 *frame_wave = c1->cd(1)->DrawFrame(-19.9,-19.9,19.9,19.9);
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

	TGraph2D* disp = GetTGraph2D(Y, Z, DY);
	disp->SetTitle("#Delta y;y_{0};z_{0}");
	disp->SetNpx(500);
	disp->SetNpy(500);
	
	pad2->cd();

	pad2->SetRightMargin( 0.14);
	pad2->SetLeftMargin(  0.06);
	pad2->SetBottomMargin(0.08);
	pad2->SetTopMargin(   0.08);

	TH1 *frame_disp = c1->cd(2)->DrawFrame(-19.9,-19.9,19.9,19.9);
	frame_disp->SetTitle("#Delta r;y_{0};z_{0}");
	frame_disp->GetXaxis()->SetLabelSize(0.03);
	frame_disp->GetXaxis()->SetTitleSize(0.06);
	frame_disp->GetXaxis()->SetTitleOffset(0.5);
	frame_disp->GetYaxis()->SetLabelSize(0.03);
	frame_disp->GetYaxis()->SetTitleSize(0.06);
	frame_disp->GetYaxis()->SetTitleOffset(0.3);
	frame_disp->GetZaxis()->SetMaxDigits(2);
	disp->Draw("same COLZ");

	pad2->SetGrid(1, 1);
	pad2->RedrawAxis("g");

	//=====================================

	TGraph2D* pym = GetTGraph2D(Y, Z, PYM);
	pym->SetTitle("p_{y}^{max};y_{0};z_{0}");
	pym->SetNpx(500);
	pym->SetNpy(500);
	
	pad3->cd();

	pad3->SetRightMargin( 0.14);
	pad3->SetLeftMargin(  0.06);
	pad3->SetBottomMargin(0.08);
	pad3->SetTopMargin(   0.08);

	TH1 *frame_pym = c1->cd(3)->DrawFrame(-19.9,-19.9,19.9,19.9);
	frame_pym->SetTitle("#Delta #theta;y_{0};z_{0}");
	frame_pym->GetXaxis()->SetLabelSize(0.03);
	frame_pym->GetXaxis()->SetTitleSize(0.06);
	frame_pym->GetXaxis()->SetTitleOffset(0.5);
	frame_pym->GetYaxis()->SetLabelSize(0.03);
	frame_pym->GetYaxis()->SetTitleSize(0.06);
	frame_pym->GetYaxis()->SetTitleOffset(0.3);
	frame_pym->GetZaxis()->SetMaxDigits(2);
	pym->Draw("same COLZ");

	pad3->SetGrid(1, 1);
	pad3->RedrawAxis("g");


	//=====================================


	pad7->cd();

	TLatex* bottom_string = new TLatex(0.01, 0.4, bottom.c_str());
	bottom_string->SetTextSize(0.4);
	bottom_string->Draw();


	//=====================================


	c1->cd(); 
	pad1->Draw();
	pad2->Draw();
	pad3->Draw();
	pad4->Draw();
	pad5->Draw();
	pad6->Draw();
	pad7->Draw();

	c1->SaveAs((file + ".png").c_str());


	delete c1;
	delete fun;
	delete disp;
	delete pym;

	return 0;

}

int main(){

	mood("Output_Laguerre_wv3_data/Output_wv3_p0l0_px-0.000000"    , 0, 0, "Laguerre Gaussian Mode ( p = 0 | l = 0 ) without Longitudinal B field; Initial Longitudinal Momentum = 0; w0 = 5; E0 = 1; lambda = 1; Linear Polarization on y");
	/*mood("Output_Laguerre_wv3_data/Output_wv3_p0l1_px-0.000000"    , 0, 1, "Laguerre Gaussian Mode ( p = 0 | l = 1 ) without Longitudinal B field; Initial Longitudinal Momentum = 0; w0 = 5; E0 = 1; lambda = 1; Linear Polarization on y");
	mood("Output_Laguerre_wv3_data/Output_wv3_p0l2_px-0.000000"    , 0, 2, "Laguerre Gaussian Mode ( p = 0 | l = 2 ) without Longitudinal B field; Initial Longitudinal Momentum = 0; w0 = 5; E0 = 1; lambda = 1; Linear Polarization on y");
	mood("Output_Laguerre_wv3_data/Output_wv3_p0l3_px-0.000000"    , 0, 3, "Laguerre Gaussian Mode ( p = 0 | l = 3 ) without Longitudinal B field; Initial Longitudinal Momentum = 0; w0 = 5; E0 = 1; lambda = 1; Linear Polarization on y");
	mood("Output_Laguerre_wv3_data/Output_wv3_p1l0_px-0.000000"    , 1, 0, "Laguerre Gaussian Mode ( p = 1 | l = 0 ) without Longitudinal B field; Initial Longitudinal Momentum = 0; w0 = 5; E0 = 1; lambda = 1; Linear Polarization on y");
	mood("Output_Laguerre_wv3_data/Output_wv3_p1l1_px-0.000000"    , 1, 1, "Laguerre Gaussian Mode ( p = 1 | l = 1 ) without Longitudinal B field; Initial Longitudinal Momentum = 0; w0 = 5; E0 = 1; lambda = 1; Linear Polarization on y");
	mood("Output_Laguerre_wv3_data/Output_wv3_p1l2_px-0.000000"    , 1, 2, "Laguerre Gaussian Mode ( p = 1 | l = 2 ) without Longitudinal B field; Initial Longitudinal Momentum = 0; w0 = 5; E0 = 1; lambda = 1; Linear Polarization on y");
	mood("Output_Laguerre_wv3_data/Output_wv3_p1l3_px-0.000000"    , 1, 3, "Laguerre Gaussian Mode ( p = 1 | l = 3 ) without Longitudinal B field; Initial Longitudinal Momentum = 0; w0 = 5; E0 = 1; lambda = 1; Linear Polarization on y");
	mood("Output_Laguerre_wv3_data/Output_wv3_p2l0_px-0.000000"    , 2, 0, "Laguerre Gaussian Mode ( p = 2 | l = 0 ) without Longitudinal B field; Initial Longitudinal Momentum = 0; w0 = 5; E0 = 1; lambda = 1; Linear Polarization on y");
	mood("Output_Laguerre_wv3_data/Output_wv3_p2l1_px-0.000000"    , 2, 1, "Laguerre Gaussian Mode ( p = 2 | l = 1 ) without Longitudinal B field; Initial Longitudinal Momentum = 0; w0 = 5; E0 = 1; lambda = 1; Linear Polarization on y");
	mood("Output_Laguerre_wv3_data/Output_wv3_p2l2_px-0.000000"    , 2, 2, "Laguerre Gaussian Mode ( p = 2 | l = 2 ) without Longitudinal B field; Initial Longitudinal Momentum = 0; w0 = 5; E0 = 1; lambda = 1; Linear Polarization on y");
	mood("Output_Laguerre_wv3_data/Output_wv3_p2l3_px-0.000000"    , 2, 3, "Laguerre Gaussian Mode ( p = 2 | l = 3 ) without Longitudinal B field; Initial Longitudinal Momentum = 0; w0 = 5; E0 = 1; lambda = 1; Linear Polarization on y");
	mood("Output_Laguerre_wv3_data/Output_wv3_p3l0_px-0.000000"    , 3, 0, "Laguerre Gaussian Mode ( p = 3 | l = 0 ) without Longitudinal B field; Initial Longitudinal Momentum = 0; w0 = 5; E0 = 1; lambda = 1; Linear Polarization on y");
	mood("Output_Laguerre_wv3_data/Output_wv3_p3l1_px-0.000000"    , 3, 1, "Laguerre Gaussian Mode ( p = 3 | l = 1 ) without Longitudinal B field; Initial Longitudinal Momentum = 0; w0 = 5; E0 = 1; lambda = 1; Linear Polarization on y");
	mood("Output_Laguerre_wv3_data/Output_wv3_p3l2_px-0.000000"    , 3, 2, "Laguerre Gaussian Mode ( p = 3 | l = 2 ) without Longitudinal B field; Initial Longitudinal Momentum = 0; w0 = 5; E0 = 1; lambda = 1; Linear Polarization on y");
	mood("Output_Laguerre_wv3_data/Output_wv3_p3l3_px-0.000000"    , 3, 3, "Laguerre Gaussian Mode ( p = 3 | l = 3 ) without Longitudinal B field; Initial Longitudinal Momentum = 0; w0 = 5; E0 = 1; lambda = 1; Linear Polarization on y");
*/
}

