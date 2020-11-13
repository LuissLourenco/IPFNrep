#include "../src/DataAnalysis.cpp"
#include <dirent.h>
#include <string>
#include <vector>
#include "TLine.h"
#include "TApplication.h"

using namespace std;

int main(int argc, char** argv){

	int n_points, n_cols;
	double** values = ReadFile("Data_tfwhm.txt", &n_cols, &n_points, true);

	DataSet RISE = DataSet(n_points, values[0]);
	DataSet RI = DataSet(n_points, values[2]);
	DataSet TE = DataSet(n_points, values[3]);
	DataSet RM = DataSet(n_points, values[4]);
	DataSet PX = DataSet(n_points, values[5]);
	DataSet XF = DataSet(n_points, values[6]);
	DataSet ET = DataSet(n_points, values[7]);

	TCanvas* c1 = new TCanvas("c", "", 5000, 1000);
	//c1->SetRightMargin(0.0);
	//c1->SetLeftMargin(0.0);
	//c1->SetBottomMargin(0.0);
	//c1->SetTopMargin(0.0);

	c1->Divide(5, 1);

	c1->cd(1);
	TGraph* g1 = GetTGraph(RISE, TE);
	g1->SetMarkerStyle(20);
	g1->SetMarkerSize(0.2);
	g1->SetTitle("Period Elipse - tfwhm");
	g1->Draw("ALP");


	c1->cd(2);
	TGraph* g2 = GetTGraph(RISE, RM);
	g2->SetMarkerStyle(20);
	g2->SetMarkerSize(0.2);
	g2->SetTitle("Raio Max - tfwhm");
	g2->Draw("ALP");


	c1->cd(3);
	TGraph* g3 = GetTGraph(RISE, PX);
	g3->SetMarkerStyle(20);
	g3->SetMarkerSize(0.2);
	g3->SetTitle("PX_Mean - tfwhm");
	g3->Draw("ALP");


	c1->cd(4);
	TGraph* g4 = GetTGraph(RISE, XF);
	g4->SetMarkerStyle(20);
	g4->SetMarkerSize(0.2);
	g4->SetTitle("X final - tfwhm");
	g4->Draw("ALP");


	c1->cd(5);
	TGraph* g5 = GetTGraph(RISE, ET);
	g5->SetMarkerStyle(20);
	g5->SetMarkerSize(0.2);
	g5->SetTitle("Eta - tfwhm");
	g5->Draw("ALP");


	c1->SaveAs("Quick_Plot.png");

	return 0;
}