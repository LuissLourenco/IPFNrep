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

int main(){
	
	int n_points, n_cols;
	double** values = ReadFile("../final_outputs/Data_Out_Finner.txt", &n_cols, &n_points, true);
	// file phi	r	period	raio_max	px_mean	x_final	eta
	DataSet PER(n_points, values[3]);
	DataSet PXM(n_points, values[5]);
	DataSet ETA(n_points, values[7]);

	values = ReadFile("../final_outputs/Logger_Finner", &n_cols, &n_points, true);
	// process	r	phi	px0	kdamp	T	N	pri	wave_type	tfwhm	stable	Eo	delta	w0	lambda	l	p
	DataSet PX0(n_points, values[3]);
	DataSet A0(n_points, values[11]);

	DataSet PER1, PXM1, ETA1, PX01, A01;
	for(int i = 0; i < PER.size(); i++){
		if(PER[i].val() < 0) continue;
		if(PX0[i] > A0[i]*Var(-0.2)+Var(1)) continue;
		if(PX0[i] > Var(-0.8)) continue;
		if(abs(ETA[i] - Var(30)) > 1) continue;
		PER1.append(PER[i]);
		PXM1.append(PXM[i]);
		ETA1.append(ETA[i]);
		PX01.append(PX0[i]);
		A01.append(A0[i]);
	}

	DataSet CUSTO = PER1 * PER1 * abs(PXM1);

	gStyle->SetPalette(kCMYK);
    gStyle->SetNumberContours(999);

    TApplication *MyApp = new TApplication("MyApp", NULL, NULL);

	TCanvas* c1 = new TCanvas("c", "", 4500, 3150);

	//TGraph2D *g = GetTGraph2D(A01, PX01, log(ETA1/CUSTO)+Var(13));
	TGraph2D *g = GetTGraph2D(A01, PX01, CUSTO);
	g->SetNpx(500);
	g->SetNpy(500);
	g->SetTitle(";a_{0}; p_{x0}");

	//g->Draw("TRI1");

	
	c1->cd();

	c1->SetRightMargin( 0.14);
	c1->SetLeftMargin(  0.06);
	c1->SetBottomMargin(0.08);
	c1->SetTopMargin(   0.08);

	TH1 *frame_graph = c1->cd(3)->DrawFrame(A0.getMin().val(), PX0.getMin().val(), A0.getMax().val(), PX0.getMax().val());
	frame_graph->GetXaxis()->SetLabelSize(0.03);
	frame_graph->GetXaxis()->SetTitleSize(0.06);
	frame_graph->GetXaxis()->SetTitleOffset(0.5);
	frame_graph->GetYaxis()->SetLabelSize(0.03);
	frame_graph->GetYaxis()->SetTitleSize(0.06);
	frame_graph->GetYaxis()->SetTitleOffset(0.3);
	frame_graph->GetZaxis()->SetMaxDigits(2);	
	frame_graph->SetTitle(";a_{0}; p_{x0}");

	g->Draw("same COLZ");

	c1->SetGrid(1, 1);
	c1->RedrawAxis("g");

	

	c1->SaveAs("plot.png");

	MyApp->Run();

	system("xdg-open plot.png");


	return 0;
}

