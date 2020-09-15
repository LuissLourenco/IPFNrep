#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <chrono>
#include <ctime> 

#include "theory4.c"
#include "DataAnalysis.cpp"

#include "TCanvas.h"
#include "TGraph2D.h"
#include "TApplication.h"
#include "TPolyLine3D.h"

using namespace std;



void input(double x01, double x02, double x03, double p01, double p02, double p03, double kdamp, double T, int N, int pri, double dx, int wave_type, double tfwhm, double stable, double Eo, double delta, double w0, double lambda, double n, double eta, int l, int p){
	FILE* itb = fopen("../../rk4_cooling/InputToBatch.txt","w");
	fprintf(itb, "trash|x01|x02|x03|p01|p02|p03|kdamp|T|N|pri|dx|wave_type|tfwhm|stable|E0|delta|w0|lambda|n|eta|l|p\n"); 
	fprintf(itb, "%.10e\n", x01);	//x01
	fprintf(itb, "%.10e\n", x02);	//x02
	fprintf(itb, "%.10e\n", x03);	//x03
	fprintf(itb, "%.10e\n", p01);	//p01
	fprintf(itb, "%.10e\n", p02);	//p02
	fprintf(itb, "%.10e\n", p03);	//p03
	fprintf(itb, "%.10e\n", kdamp);	//kdamp
	fprintf(itb, "%.10e\n", T);	//T
	fprintf(itb, "%i\n", N);	//N
	fprintf(itb, "%i\n", pri);	//pri
	fprintf(itb, "%.10e\n", dx);	//dx
	fprintf(itb, "%i\n", wave_type);	//wave_type
	fprintf(itb, "%.10e\n", tfwhm);	//tfwhm
	fprintf(itb, "%.10e\n", stable);	//stable
	fprintf(itb, "%.10e\n", Eo);	//Eo
	fprintf(itb, "%.10e\n", delta);	//delta
	fprintf(itb, "%.10e\n", w0);	//w0
	fprintf(itb, "%.10e\n", lambda);	//lambda
	fprintf(itb, "%.10e\n", n);	//n
	fprintf(itb, "%.10e\n", eta);	//eta
	fprintf(itb, "%i\n", l);	//l
	fprintf(itb, "%i", p);	//p
	fclose(itb);
}


int main(int argc, char **argv){

	/*
	int n_cols, n_points;
	double** values;
	double* bog;
	DataSet X(0,bog);
	DataSet Y(0,bog);
	DataSet Z(0,bog);
	*/

	int ns = 20;
	//double a0=M_PI/2;
	double r0=3;
	double a0=100;
	double p0=0;
	int p=2;
	int l=0;
	double phi0=0;


	double T=1000;
	int N=1000000;

	if(argc==7){//argv-> a0 r0 p0 p l ns
		sscanf(argv[1],"%lf",&a0);
		sscanf(argv[2],"%lf",&r0);
		//sscanf(argv[2],"%lf",&phi0);
		//phi0*= M_PI/180.;
		sscanf(argv[3],"%lf",&p0);
		sscanf(argv[4],"%i",&p);
		sscanf(argv[5],"%i",&l);
		sscanf(argv[6],"%i",&ns);
	}

	for(int i=1; i<=ns; i++){

		
		//input(0,0.5*i*cos(phi0),0.5*i*sin(phi0), p0,0,0,   0, 300, 100000, 10, 5e-3, 3, 50, 1000, a0, -1, 5, 1, 1, 1, l, p);
		input(0,r0*cos(i*2*M_PI/ns),r0*sin(i*2*M_PI/ns),  p0,0,0,   0, T, N, 10, 5e-3, 3, 50, 1000, a0, -1, 5, 1, 1, 1, l, p);
		run_theory3(i);
		/*
		values = ReadFile((string("Outputs/Out")+to_string(i)+string(".txt")).c_str(), &n_cols, &n_points, false, false);
		X = X.concat(DataSet(n_points,values[1]));
		Y = Y.concat(DataSet(n_points,values[2]));
		Z = Z.concat(DataSet(n_points,values[3]));
		*/
		cout<<i<<"/"<<ns<<"\t"<<flush;
		if(i%5==0) cout<<endl;
	}
	cout<<endl<<"gotovo"<<endl;





	/*
	TGraph2D* graph_y2 = GetTGraph2D(X,Y,Z);
	TGraph2D* graph_p2 = new TGraph2D(n_points,values[4],values[5],values[6]);

	graph_y2->SetTitle(";x;y;z");
	graph_p2->SetTitle(";px;py;pz");
	graph_y2->SetMarkerColor(kRed);
	graph_p2->SetMarkerColor(kRed);

	TH2F* h2 = new TH2F("h2", ";x;y;z", 1, 0, 5, 1, -15, 15);
	graph_y2->SetHistogram(h2);
	graph_y2->SetMinimum(-15);
	graph_y2->SetMaximum(15);

	graph_y2->GetXaxis()->SetLabelSize(0.04);
	graph_y2->GetYaxis()->SetLabelSize(0.04);
	graph_y2->GetZaxis()->SetLabelSize(0.04);
	graph_y2->GetXaxis()->SetTitleSize(0.05);
	graph_y2->GetYaxis()->SetTitleSize(0.05);
	graph_y2->GetZaxis()->SetTitleSize(0.05);
	graph_y2->GetXaxis()->SetTitleOffset(1.4);
	graph_y2->GetYaxis()->SetTitleOffset(1.4);
	graph_y2->GetZaxis()->SetTitleOffset(1.1);
	gStyle->SetOptStat(0);

	TPolyLine3D* pl3d1 = new TPolyLine3D(50);
	double w0=5;
	for(int i=0; i<50; i++) pl3d1->SetPoint(i, 0, w0*cos(i*2*M_PI/49), w0*sin(i*2*M_PI/49));
	TPolyLine3D* pl3d2 = new TPolyLine3D(2);
	TPolyLine3D* pl3d3 = new TPolyLine3D(2);
	pl3d2->SetPoint(0,0,-15,0); pl3d2->SetPoint(1,0,15,0);
	pl3d3->SetPoint(0,0,0,-15); pl3d3->SetPoint(1,0,0,15);

	pl3d1->SetLineColor(kGray+1);
	pl3d2->SetLineColor(kGray);
	pl3d3->SetLineColor(kGray);

	TApplication* MyRootApp;
	MyRootApp = new TApplication("MyRootApp", NULL, NULL);

	TCanvas* c1 = new TCanvas("c1", "", 1500, 1000);
	c1->Divide(2, 1);	
	c1->cd(1);
	graph_y2->Draw("P");
	pl3d1->Draw("SAME");
	pl3d2->Draw("SAME");
	pl3d3->Draw("SAME");
	c1->cd(2);
	graph_p2->Draw("P");
	c1->cd(1);

	gPad->WaitPrimitive();
	gPad->SaveAs("Plot.png");

	delete h2;
	delete graph_y2;
	delete graph_p2;
	delete c1;
	*/




	return 0;
}

