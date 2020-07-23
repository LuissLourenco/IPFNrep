#include <iostream>
#include <chrono> 
using namespace std::chrono; 

#include <ODEsolver.h>

#include <TFormula.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TMultiGraph.h>

using namespace std;


int main(int argc, char** argv){


	double* p0 = new double[6];
	for(int i=0; i<6; i++) 
		p0[i] = 0;
	ODEpoint point0(0, p0, 6);

	vector<TFormula> form;
	//a0 = 1
	//delta = 1 -> Pol Linear em Y
	form.push_back(TFormula("","x[3]"));
	form.push_back(TFormula("","x[4]"));
	form.push_back(TFormula("","x[5]"));
	//polarizacao circular
	form.push_back(TFormula("","-pow(1-x[3]*x[3]-x[4]*x[4]-x[5]*x[5], 3/2)*(-x[4]*sin(x[0]-x[6])/sqrt(2)+x[5]*cos(x[0]-x[6])/sqrt(2))"));
	form.push_back(TFormula("","-pow(1-x[3]*x[3]-x[4]*x[4]-x[5]*x[5], 3/2)*(-sin(x[0]-x[6])/sqrt(2)+x[3]*sin(x[0]-x[6])/sqrt(2))"));
	form.push_back(TFormula("","-pow(1-x[3]*x[3]-x[4]*x[4]-x[5]*x[5], 3/2)*(cos(x[0]-x[6])/sqrt(2)-x[3]*cos(x[0]-x[6])/sqrt(2))"));
	ODEsolver solver(form);


	auto start = high_resolution_clock::now(); 

	vector<ODEpoint> solucao = solver.RK4solver(point0, 0, 200, 0.001);

	auto stop = high_resolution_clock::now(); 
	auto duration = duration_cast<microseconds>(stop - start); 
	cout << " *** " << duration.count() << " ms ***" << endl; 


	int n = solucao.size();
	TGraph* g1 = new TGraph(n);
	TGraph* g2 = new TGraph(n);
	TGraph* g3 = new TGraph(n);
	TGraph* g4 = new TGraph(n);
	TGraph* g5 = new TGraph(n);
	for(int i=0; i<n; i++){
		g1->SetPoint(i, solucao[i].Get_VarTime()[0], solucao[i].Get_VarTime()[1]);
		g2->SetPoint(i, solucao[i].Get_VarTime()[0], solucao[i].Get_VarTime()[2]);
		g3->SetPoint(i, solucao[i].Get_VarTime()[6], solucao[i].Get_VarTime()[0]);
		g4->SetPoint(i, solucao[i].Get_VarTime()[6], solucao[i].Get_VarTime()[1]);
		g5->SetPoint(i, solucao[i].Get_VarTime()[6], solucao[i].Get_VarTime()[2]);
	}

	g1->SetMarkerColor(kBlue);
	g1->SetLineColor(kBlue);
	g2->SetMarkerColor(kBlue);
	g2->SetLineColor(kBlue);
	g3->SetMarkerColor(kBlue);
	g3->SetLineColor(kBlue);
	g4->SetMarkerColor(kBlue);
	g4->SetLineColor(kBlue);
	g5->SetMarkerColor(kViolet);
	g5->SetLineColor(kViolet);

	g1->SetTitle(";x;y");
	g2->SetTitle(";x;z");
	g3->SetTitle(";t;x");
	g4->SetTitle(";t;y");
	g5->SetTitle(";t;z");

	TCanvas* c1 = new TCanvas();	
	g1->Draw("APL");
	c1->SaveAs("plot1.png");

	TCanvas* c2 = new TCanvas();	
	g2->Draw("APL");
	c2->SaveAs("plot2.png");

	TCanvas* c3 = new TCanvas();	
	g3->Draw("APL");
	c3->SaveAs("plot3.png");

	TCanvas* c4 = new TCanvas();
	TMultiGraph* mg1 = new TMultiGraph();
	mg1->Add(g4);
	mg1->Add(g5);	
	mg1->SetTitle(";t;y,z");
	mg1->Draw("APL");
	c4->SaveAs("plot4.png");

	delete[] p0;
	delete g1;
	delete g2;
	delete g3;
	delete g4;
	delete g5;
	delete c1;
	delete c2;
	delete c3;
	delete c4;
	delete mg1;

  	return 0;
}
