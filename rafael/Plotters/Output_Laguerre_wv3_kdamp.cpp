#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1

#include <math.h>
#include <iostream>

#include "DataAnalysis.cpp"
#include "TVectorField.cpp"

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
	
	DataSet PYMAX(n_points, values[8]);
	DataSet PZMAX(n_points, values[9]);

	DataSet E(n_points, values[10]);

	DataSet LXMAX(n_points, values[11]);
	DataSet LXF(n_points, values[12]);

	DataSet DY = DataSet(n_points, values[3]) - Y;
	DataSet DZ = DataSet(n_points, values[4]) - Z;


	TCanvas* c1 = new TCanvas("c", "", 4500, 3150);
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

	TGraph2D* graph2 = GetTGraph2D(Y, Z, PYMAX);
	graph2->SetTitle("Maximum y Momentum;y_{0};z_{0}");
	graph2->SetNpx(500);
	graph2->SetNpy(500);
	
	pad2->cd();

	pad2->SetRightMargin( 0.14);
	pad2->SetLeftMargin(  0.06);
	pad2->SetBottomMargin(0.08);
	pad2->SetTopMargin(   0.08);

	TH1 *frame_graph2 = c1->cd(2)->DrawFrame(-19.9,-19.9,19.9,19.9);
	frame_graph2->SetTitle("Maximum y Momentum;y_{0};z_{0}");
	frame_graph2->GetXaxis()->SetLabelSize(0.03);
	frame_graph2->GetXaxis()->SetTitleSize(0.06);
	frame_graph2->GetXaxis()->SetTitleOffset(0.5);
	frame_graph2->GetYaxis()->SetLabelSize(0.03);
	frame_graph2->GetYaxis()->SetTitleSize(0.06);
	frame_graph2->GetYaxis()->SetTitleOffset(0.3);
	frame_graph2->GetZaxis()->SetMaxDigits(2);
	graph2->Draw("same COLZ");

	pad2->SetGrid(1, 1);
	pad2->RedrawAxis("g");

	//=====================================

	TGraph2D* graph3 = GetTGraph2D(Y, Z, LXMAX);
	graph3->SetTitle("Maximum x Angular Momentum;y_{0};z_{0}");
	graph3->SetNpx(500);
	graph3->SetNpy(500);
	
	pad3->cd();

	pad3->SetRightMargin( 0.14);
	pad3->SetLeftMargin(  0.06);
	pad3->SetBottomMargin(0.08);
	pad3->SetTopMargin(   0.08);

	TH1 *frame_graph3 = c1->cd(3)->DrawFrame(-19.9,-19.9,19.9,19.9);
	frame_graph3->SetTitle("Maximum x Angular Momentum;y_{0};z_{0}");
	frame_graph3->GetXaxis()->SetLabelSize(0.03);
	frame_graph3->GetXaxis()->SetTitleSize(0.06);
	frame_graph3->GetXaxis()->SetTitleOffset(0.5);
	frame_graph3->GetYaxis()->SetLabelSize(0.03);
	frame_graph3->GetYaxis()->SetTitleSize(0.06);
	frame_graph3->GetYaxis()->SetTitleOffset(0.3);
	frame_graph3->GetZaxis()->SetMaxDigits(2);
	graph3->Draw("same COLZ");

	pad3->SetGrid(1, 1);
	pad3->RedrawAxis("g");


	//=====================================

	TVectorField* graph4 = new TVectorField(Y, Z, DY, DZ);
	//graph4->SetTitle("Maximum x Angular Momentum;y_{0};z_{0}");
	//graph4->SetNpx(500);
	//graph4->SetNpy(500);
	graph4->SetLimits(-19.9, -19.9, 19.9, 19.9);
	graph4->SetArrowColor(1);
	graph4->SetArrowSize(0.5, 0.002, 2);

	pad4->cd();

	pad4->SetRightMargin( 0.14);
	pad4->SetLeftMargin(  0.06);
	pad4->SetBottomMargin(0.08);
	pad4->SetTopMargin(   0.08);

	graph4->Draw("F");

	TH1 *frame_graph4 = graph4->GetTH1();
	frame_graph4->SetTitle("Displacement;y_{0};z_{0}");
	frame_graph4->SetTitleSize(100);
	//frame_graph4->SetTitleOffset(0.01);
	frame_graph4->GetXaxis()->SetLabelSize(0.03);
	frame_graph4->GetXaxis()->SetTitleSize(0.06);
	frame_graph4->GetXaxis()->SetTitleOffset(0.5);
	frame_graph4->GetYaxis()->SetLabelSize(0.03);
	frame_graph4->GetYaxis()->SetTitleSize(0.06);
	frame_graph4->GetYaxis()->SetTitleOffset(0.3);
	frame_graph4->GetZaxis()->SetMaxDigits(2);
	graph4->ReDraw("F");

	pad4->SetGrid(1, 1);
	pad4->RedrawAxis("g");


	//=====================================

	TGraph2D* graph5 = GetTGraph2D(Y, Z, PZMAX);
	graph5->SetTitle("Maximum z Momentum;y_{0};z_{0}");
	graph5->SetNpx(500);
	graph5->SetNpy(500);
	
	pad5->cd();

	pad5->SetRightMargin( 0.14);
	pad5->SetLeftMargin(  0.06);
	pad5->SetBottomMargin(0.08);
	pad5->SetTopMargin(   0.08);

	TH1 *frame_graph5 = c1->cd(3)->DrawFrame(-19.9,-19.9,19.9,19.9);
	frame_graph5->SetTitle("Maximum z Momentum;y_{0};z_{0}");
	frame_graph5->GetXaxis()->SetLabelSize(0.03);
	frame_graph5->GetXaxis()->SetTitleSize(0.06);
	frame_graph5->GetXaxis()->SetTitleOffset(0.5);
	frame_graph5->GetYaxis()->SetLabelSize(0.03);
	frame_graph5->GetYaxis()->SetTitleSize(0.06);
	frame_graph5->GetYaxis()->SetTitleOffset(0.3);
	frame_graph5->GetZaxis()->SetMaxDigits(2);
	graph5->Draw("same COLZ");

	pad5->SetGrid(1, 1);
	pad5->RedrawAxis("g");


	//=====================================

	TGraph2D* graph6 = GetTGraph2D(Y, Z, LXF);
	graph6->SetTitle("Final x Angular Momentum;y_{0};z_{0}");
	graph6->SetNpx(500);
	graph6->SetNpy(500);
	
	pad6->cd();

	pad6->SetRightMargin( 0.14);
	pad6->SetLeftMargin(  0.06);
	pad6->SetBottomMargin(0.08);
	pad6->SetTopMargin(   0.08);

	TH1 *frame_graph6 = c1->cd(3)->DrawFrame(-19.9,-19.9,19.9,19.9);
	frame_graph6->SetTitle("Final x Angular Momentum;y_{0};z_{0}");
	frame_graph6->GetXaxis()->SetLabelSize(0.03);
	frame_graph6->GetXaxis()->SetTitleSize(0.06);
	frame_graph6->GetXaxis()->SetTitleOffset(0.5);
	frame_graph6->GetYaxis()->SetLabelSize(0.03);
	frame_graph6->GetYaxis()->SetTitleSize(0.06);
	frame_graph6->GetYaxis()->SetTitleOffset(0.3);
	frame_graph6->GetZaxis()->SetMaxDigits(2);
	graph6->Draw("same COLZ");

	pad6->SetGrid(1, 1);
	pad6->RedrawAxis("g");


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
	delete graph2;
	delete graph3;

	return 0;

}

int main(){

	mood("Output_Laguerre_wv3_kdamp_data/Output_wv3_l1p0_px-2000.000000" , 0, 1, "Simplified Laguerre Gaussian Mode ( p = 0 | l = 1 ); Initial Longitudinal Momentum = -2000; w0 = 5; E0 = 1; lambda = 1; Linear Polarization on y");
	mood("Output_Laguerre_wv3_kdamp_data/Output_wv3_l2p0_px-2000.000000" , 0, 2, "Simplified Laguerre Gaussian Mode ( p = 0 | l = 2 ); Initial Longitudinal Momentum = -2000; w0 = 5; E0 = 1; lambda = 1; Linear Polarization on y");
	
}

