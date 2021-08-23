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
#include "TApplication.h"

using namespace std;

void ellipsis_period(){

	int n_points, n_cols;
	double** values = ReadFile("../final_outputs/Data_Out_Finner.txt", &n_cols, &n_points, true);
	// file phi	r	period	raio_max	px_mean	x_final	eta
	DataSet PER_F(n_points, values[3]);
	DataSet PXM_F(n_points, values[5]);
	DataSet ETA_F(n_points, values[7]);

	//ETA_F.round(0);

	values = ReadFile("../final_outputs/Logger_Finner", &n_cols, &n_points, true);
	// process	r	phi	px0	kdamp	T	N	pri	wave_type	tfwhm	stable	Eo	delta	w0	lambda	l	p
	DataSet PX0_F(n_points, values[3]);
	DataSet A0_F(n_points, values[11]);

	DataSet PER, PXM, ETA, PX0, A0, PARAM;
	for(int i = 0; i < PER_F.size(); i++){
		if(PER_F[i].val() < 0) continue;
		if(PX0_F[i] > A0_F[i]*Var(-0.2)+Var(0.5)) continue;
		if(PX0_F[i] > Var(-0.8)) continue;
		//if(abs(ETA_F[i] - Var(20)) > 1) continue;
		PER.append(PER_F[i]);
		PXM.append(PXM_F[i]);
		ETA.append(ETA_F[i]);
		PX0.append(PX0_F[i]);
		A0.append(A0_F[i]);
		PARAM.append(abs(PX0_F[i])*A0_F[i]);
		//PARAM.append(Var(A0.size()));
	}

	DataSet CUSTO = PER * PER * abs(PXM);

	//gStyle->SetPalette(kCMYK);
    gStyle->SetNumberContours(999);

	TCanvas* c1 = new TCanvas("c", "", 2000, 2000);

	//TGraph2D *g = GetTGraph2D(A0, PX0, log(ETA/CUSTO)+Var(13));
	//TGraph2D *g = GetTGraph2D(A0, PX0, CUSTO);
	TGraph2D *g = GetTGraph2D(A0, PX0, PER/ETA);
	g->SetNpx(500);
	g->SetNpy(500);
	g->SetTitle(";a_{0}; p_{x0}");


	
	c1->cd();

	c1->SetRightMargin( 0.14);
	c1->SetLeftMargin(  0.06);
	c1->SetBottomMargin(0.08);
	c1->SetTopMargin(   0.08);

	TH1 *frame_graph = c1->cd(3)->DrawFrame(A0.getMin().val()-0.5, PX0.getMin().val()-0.5, A0.getMax().val()+0.5, PX0.getMax().val()+0.5);
	frame_graph->GetXaxis()->SetLabelSize(0.03);
	frame_graph->GetXaxis()->SetTitleSize(0.05);
	frame_graph->GetXaxis()->SetTitleOffset(0.7);
	frame_graph->GetYaxis()->SetLabelSize(0.03);
	frame_graph->GetYaxis()->SetTitleSize(0.05);
	frame_graph->GetYaxis()->SetTitleOffset(0.4);
	frame_graph->GetZaxis()->SetMaxDigits(2);	
	frame_graph->SetTitle("Period of each Ellipsis;a_{0};p_{x0};T_{ell}(T_{las})");

	g->Draw("same COLZ");

	c1->SetGrid(1, 1);
	c1->RedrawAxis("g");

	auto fit_aux = [=](Double_t *x, Double_t *par) { // x[0] - v; x[1] - t
		//return par[0] + par[1] * x[0] + par[2] * x[1];
		return par[0] * atan2(x[1], x[0]) + par[1];
		//return par[0] * x[1]/x[0] + par[1];
	};
	
	TF2* fit = new TF2("fit2", fit_aux, A0.getMin().val()-0.5, A0.getMax().val()+0.5, PX0.getMin().val()-0.5, PX0.getMax().val()+0.5, 2);
	fit->SetParameters(1, 1, 1);
	g->Fit(fit, "R");

	DataSet AUX = PER/ETA;

	int pad_levels = 5;
	int n_levels = (AUX.getMax() - AUX.getMin()).val() + 2 * pad_levels;
	double *levels = new double[n_levels];
	for(int i = 0; i < n_levels; i++) 
		levels[i] = AUX.getMin().val() + i - pad_levels;
	fit->SetContour(n_levels, levels);
	fit->SetLineColor(kRed);

	fit->SetNpx(1000);
	fit->SetNpy(1000);
	fit->SetLineWidth(2);

	fit->Draw("SAME");

	char aux[128];
	sprintf(aux, "T_{ell} = %3.2lf #times atan #frac{p_{x0}}{a_{0}} + %3.2lf (T_{las})", fit->GetParameter(0), fit->GetParameter(1));
	TLatex* latex1 = new TLatex(2, -10, aux);
	latex1->SetTextSize(0.03);
	latex1->Draw();

	sprintf(aux, "#chi^{2}/NDF = %4.3lf", fit->GetChisquare() / fit->GetNDF());
	TLatex* latex2 = new TLatex(2, -9.5, aux);
	latex2->SetTextSize(0.03);
	latex2->Draw();

	c1->SaveAs("Ellipsis_Period.png");

	system("xdg-open Ellipsis_Period.png");

}

void cost(int eta_plot, double eta_range){

	int n_points, n_cols;
	double** values = ReadFile("../final_outputs/Data_Out_Finner.txt", &n_cols, &n_points, true);
	// file phi	r	period	raio_max	px_mean	x_final	eta
	DataSet PER_F(n_points, values[3]);
	DataSet PXM_F(n_points, values[5]);
	DataSet ETA_F(n_points, values[7]);

	//ETA_F.round(0);

	values = ReadFile("../final_outputs/Logger_Finner", &n_cols, &n_points, true);
	// process	r	phi	px0	kdamp	T	N	pri	wave_type	tfwhm	stable	Eo	delta	w0	lambda	l	p
	DataSet PX0_F(n_points, values[3]);
	DataSet A0_F(n_points, values[11]);

	DataSet PER, PXM, ETA, PX0, A0, PARAM;
	for(int i = 0; i < PER_F.size(); i++){
		if(PER_F[i].val() < 0) continue;
		if(PX0_F[i] > A0_F[i]*Var(-0.2)+Var(0.5)) continue;
		if(PX0_F[i] > Var(-0.8)) continue;
		if(abs(ETA_F[i] - Var(eta_plot)) > eta_range) continue;
		PER.append(PER_F[i]);
		PXM.append(PXM_F[i]);
		ETA.append(ETA_F[i]);
		PX0.append(PX0_F[i]);
		A0.append(A0_F[i]);
		PARAM.append(abs(PX0_F[i])*A0_F[i]);
		//PARAM.append(Var(A0.size()));
	}

	DataSet CUSTO = PER * PER * abs(PXM);

	//gStyle->SetPalette(kCMYK);
    gStyle->SetNumberContours(999);

    //TApplication *MyApp = new TApplication("MyApp", NULL, NULL);

	TCanvas* c1 = new TCanvas("c", "", 2000, 2000);

	//TGraph2D *g = GetTGraph2D(A0, PX0, log(ETA/CUSTO)+Var(13));
	//TGraph2D *g = GetTGraph2D(A0, PX0, CUSTO);
	TGraph2D *g = GetTGraph2D(A0, PX0, CUSTO);
	g->SetNpx(500);
	g->SetNpy(500);
	g->SetTitle(";a_{0}; p_{x0}");

	c1->cd();

	c1->SetRightMargin( 0.14);
	c1->SetLeftMargin(  0.06);
	c1->SetBottomMargin(0.08);
	c1->SetTopMargin(   0.08);

	//TH1 *frame_graph = c1->cd(3)->DrawFrame(A0.getMin().val()-0.5, PX0.getMin().val()-0.5, A0.getMax().val()+0.5, PX0.getMax().val()+0.5);
	TH1 *frame_graph = c1->cd(3)->DrawFrame(1.5, -10.5, 20.5, -0.3);
	frame_graph->GetXaxis()->SetLabelSize(0.03);
	frame_graph->GetXaxis()->SetTitleSize(0.05);
	frame_graph->GetXaxis()->SetTitleOffset(0.7);
	frame_graph->GetYaxis()->SetLabelSize(0.03);
	frame_graph->GetYaxis()->SetTitleSize(0.05);
	frame_graph->GetYaxis()->SetTitleOffset(0.4);
	frame_graph->GetZaxis()->SetMaxDigits(2);	
	frame_graph->SetTitle("Cost;a_{0};p_{x0}");
	
	g->Draw("same COLZ");

	c1->SetGrid(1, 1);
	c1->RedrawAxis("g");

	/*
	TGraph* g1 = GetTGraph(A0 / PX0, PER/ETA);
	g1->SetTitle("todos os freaking pontos;px0 #times a0; custo");
	g1->SetMarkerStyle(20);
	g1->SetMarkerSize(2);
	g1->Draw("APL");
	*/

	/*
	auto fit_aux = [=](Double_t *x, Double_t *par) { // x[0] - v; x[1] - t
		//return par[0] + par[1] * x[0] + par[2] * x[1];
		return par[0] * atan2(x[1], x[0]) + par[1];
	};
	
	TF2* fit = new TF2("fit2", fit_aux, A0.getMin().val()-0.5, A0.getMax().val()+0.5, PX0.getMin().val()-0.5, PX0.getMax().val()+0.5, 2);
	fit->SetParameters(1, 1, 1);
	g->Fit(fit, "R");

	DataSet AUX = PER/ETA;

	int pad_levels = 5;
	int n_levels = (AUX.getMax() - AUX.getMin()).val() + 2 * pad_levels;
	double *levels = new double[n_levels];
	for(int i = 0; i < n_levels; i++) 
		levels[i] = AUX.getMin().val() + i - pad_levels;
	fit->SetContour(n_levels, levels);
	fit->SetLineColor(kRed);

	fit->Draw("SAME");

	char aux[128];
	sprintf(aux, "T_{ell} = %3.2lf #times atan #frac{p_{x0}}{a_{0}} + %3.2lf (T_{las})", fit->GetParameter(0), fit->GetParameter(1));
	TLatex* latex1 = new TLatex(2, -10, aux);
	latex1->SetTextSize(0.03);
	latex1->Draw();

	sprintf(aux, "#chi^{2}/NDF = %4.3lf", fit->GetChisquare() / fit->GetNDF());
	TLatex* latex2 = new TLatex(2, -9.5, aux);
	latex2->SetTextSize(0.03);
	latex2->Draw();
	*/

	char title[64];
	sprintf(title, "Cost_Eta%02d.png", eta_plot);
	c1->SaveAs(title);

	//MyApp->Run();

}

int main(){
	
	//ellipsis_period();

	//cost(20, 1);
	//cost(30, 1);
	//cost(40, 1);
	//cost(50, 2);


	int n_points, n_cols;
	double** values = ReadFile("../final_outputs/Data_Out_Finner.txt", &n_cols, &n_points, true);
	// file phi	r	period	raio_max	px_mean	x_final	eta
	DataSet PER_F(n_points, values[3]);
	DataSet PXM_F(n_points, values[5]);
	DataSet ETA_F(n_points, values[7]);

	//ETA_F.round(0);

	values = ReadFile("../final_outputs/Logger_Finner", &n_cols, &n_points, true);
	// process	r	phi	px0	kdamp	T	N	pri	wave_type	tfwhm	stable	Eo	delta	w0	lambda	l	p
	DataSet PX0_F(n_points, values[3]);
	DataSet A0_F(n_points, values[11]);

	DataSet PER, PXM, ETA, PX0, A0, PARAM;
	for(int i = 0; i < PER_F.size(); i++){
		if(PER_F[i].val() < 0) continue;
		if(PX0_F[i] > A0_F[i]*Var(-0.2)+Var(0.5)) continue;
		if(PX0_F[i] > Var(-0.8)) continue;
		//if(abs(ETA_F[i] - Var(30)) > 1) continue;
		PER.append(PER_F[i]);
		PXM.append(PXM_F[i]);
		ETA.append(ETA_F[i]);
		PX0.append(PX0_F[i]);
		A0.append(A0_F[i]);
		PARAM.append(abs(PX0_F[i])*A0_F[i]);
		//PARAM.append(Var(A0.size()));
	}

	DataSet CUSTO = PER * PER * abs(PXM);

	//gStyle->SetPalette(kCMYK);
    gStyle->SetNumberContours(999);

    //TApplication *MyApp = new TApplication("MyApp", NULL, NULL);

	TCanvas* c1 = new TCanvas("c", "", 2000, 2000);

	//TGraph2D *g = GetTGraph2D(A0, PX0, log(ETA/CUSTO)+Var(13));
	//TGraph2D *g = GetTGraph2D(A0, PX0, CUSTO);
	TGraph2D *g = GetTGraph2D(A0, PX0, abs(PXM));
	g->SetNpx(500);
	g->SetNpy(500);

	c1->cd();

	c1->SetRightMargin( 0.14);
	c1->SetLeftMargin(  0.06);
	c1->SetBottomMargin(0.08);
	c1->SetTopMargin(   0.08);

	//TH1 *frame_graph = c1->cd(3)->DrawFrame(A0.getMin().val()-0.5, PX0.getMin().val()-0.5, A0.getMax().val()+0.5, PX0.getMax().val()+0.5);
	TH1 *frame_graph = c1->cd(3)->DrawFrame(1.5, -10.5, 20.5, -0.3);
	frame_graph->GetXaxis()->SetLabelSize(0.03);
	frame_graph->GetXaxis()->SetTitleSize(0.05);
	frame_graph->GetXaxis()->SetTitleOffset(0.7);
	frame_graph->GetYaxis()->SetLabelSize(0.03);
	frame_graph->GetYaxis()->SetTitleSize(0.05);
	frame_graph->GetYaxis()->SetTitleOffset(0.4);
	frame_graph->GetZaxis()->SetMaxDigits(2);	
	frame_graph->SetTitle("Average Longitudinal Momentum (m_{e}c);a_{0};p_{x0}");
	
	g->Draw("same COLZ");

	c1->SetGrid(1, 1);
	c1->RedrawAxis("g");

	/*
	TGraph* g1 = GetTGraph(A0 / PX0, PER/ETA);
	g1->SetTitle("todos os freaking pontos;px0 #times a0; custo");
	g1->SetMarkerStyle(20);
	g1->SetMarkerSize(2);
	g1->Draw("APL");
	*/

	/*
	auto fit_aux = [=](Double_t *x, Double_t *par) { // x[0] - v; x[1] - t
		//return par[0] + par[1] * x[0] + par[2] * x[1];
		return par[0] * atan2(x[1], x[0]) + par[1];
		//return par[0] * x[0] / x[1] + par[1];
	};
	
	TF2* fit = new TF2("fit2", fit_aux, A0.getMin().val()-0.5, A0.getMax().val()+0.5, PX0.getMin().val()-0.5, PX0.getMax().val()+0.5, 2);
	fit->SetParameters(1, 100);
	g->Fit(fit, "R");

	int pad_levels = 5;
	int n_levels = (PER.getMax() - PER.getMin()).val();// + 2 * pad_levels;
	double *levels = new double[n_levels];
	for(int i = 0; i < n_levels; i++) 
		levels[i] = PER.getMin().val() + i;
	fit->SetContour(n_levels, levels);
	fit->SetLineColor(kRed);

	fit->SetNpx(1000);
	fit->SetNpy(1000);
	fit->SetLineWidth(2);

	fit->Draw("SAME");

	char aux[128];
	sprintf(aux, "#T_{rot} = %3.2lf #times atan #frac{p_{x0}}{a_{0}} + %3.2lf (T_{las})", fit->GetParameter(0), fit->GetParameter(1));
	TLatex* latex1 = new TLatex(2, -10, aux);
	latex1->SetTextSize(0.03);
	latex1->Draw();

	sprintf(aux, "#chi^{2}/NDF = %4.3lf", fit->GetChisquare() / fit->GetNDF());
	TLatex* latex2 = new TLatex(2, -9.5, aux);
	latex2->SetTextSize(0.03);
	latex2->Draw();

	*/
	c1->SaveAs("MeanMomentum.png");

	//MyApp->Run();

	system("xdg-open MeanMomentum.png");

	return 0;
}


/*

MERGE THE COST FOR EACH ETA

*/
