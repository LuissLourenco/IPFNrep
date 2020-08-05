/*
#include "DataAnalysis.cpp"

void plotter(){

	int n_cols, n_points;
	double** values = ReadFile("Out0.txt", &n_cols, &n_points, false);

	DataSet X(n_points, values[1], 0.);
	DataSet Y(n_points, values[3], 0.);
	DataSet PX(n_points, values[4], 0.);
	DataSet PY(n_points, values[6], 0.);

	TCanvas* c1 = new TCanvas("c1", "", 1500, 750);
	
	TPad** pads = new TPad*[2]; 
	TGraphErrors* graph_y = GetTGraphErrors(X, Y);
	TGraphErrors* graph_p = GetTGraphErrors(PX, PY);
	//TLegend** leg = new TLegend*[6];

	pads[0] = new TPad("", "", 0.0, 0.0, 0.5, 1.0, 0, 0, 0);
	pads[0]->cd();
	pads[0]->SetGrid();
	pads[0]->SetRightMargin(0.01);
	pads[0]->SetLeftMargin(0.11);
	pads[0]->SetBottomMargin(0.11);
	pads[0]->SetTopMargin(0.01);

	graph_y->SetLineColor(2);
	graph_y->SetLineWidth(2);

	graph_y->SetTitle(";x;y");

	graph_y->GetXaxis()->SetLabelSize(0.05);
	graph_y->GetXaxis()->SetTitleSize(0.08);
	graph_y->GetXaxis()->SetTitleOffset(0.6);

	graph_y->GetYaxis()->SetLabelSize(0.05);
	graph_y->GetYaxis()->SetTitleSize(0.08);
	graph_y->GetYaxis()->SetTitleOffset(0.4);

	graph_y->Draw("AC");



	pads[1] = new TPad("", "", 0.5, 0.0, 1.0, 1.0, 0, 0, 0);
	pads[1]->cd();
	pads[1]->SetGrid();
	pads[1]->SetRightMargin(0.01);
	pads[1]->SetLeftMargin(0.11);
	pads[1]->SetBottomMargin(0.11);
	pads[1]->SetTopMargin(0.01);

	graph_p->SetLineColor(2);
	graph_p->SetLineWidth(2);

	graph_p->SetTitle(";p_{x};p_{y}");

	graph_p->GetXaxis()->SetLabelSize(0.05);
	graph_p->GetXaxis()->SetTitleSize(0.08);
	graph_p->GetXaxis()->SetTitleOffset(0.6);

	graph_p->GetYaxis()->SetLabelSize(0.05);
	graph_p->GetYaxis()->SetTitleSize(0.08);
	graph_p->GetYaxis()->SetTitleOffset(0.4);

	graph_p->Draw("AC");

	c1->cd();
	for(int i = 0; i < 2; i++) pads[i]->Draw();

	c1->SaveAs("Injected_Stopped.png");

}

*/
#include "DataAnalysis.cpp"

void plotter(){

	int n_cols, n_points;
	double** values;

	DataSet* X = new DataSet[3];
	DataSet* Y = new DataSet[3];
	DataSet* PX = new DataSet[3];
	DataSet* PY = new DataSet[3];

	for(int i = 0; i < 3; i++){
		values = ReadFile((string("Injected")+to_string(i)+string(".txt")).c_str(), &n_cols, &n_points, true);
		X[i] = DataSet(n_points, values[1], 0.);
		Y[i] = DataSet(n_points, values[3], 0.);
		PX[i] = DataSet(n_points, values[4], 0.);
		PY[i] = DataSet(n_points, values[6], 0.);
	}

	TCanvas* c1 = new TCanvas("c1", "", 1500, 1000);
	/*
	TMultiGraph* multi1 = new TMultiGraph();
	TGraphErrors* graph_y[6];
	for(int i = 0; i < 6; i++){
		graph_y[i] = GetTGraphErrors(R[i], Y[i]);
		graph_y[i]->SetLineColor(1+i);
		multi1->Add(graph_y[i]);
	}

	TMultiGraph* multi2 = new TMultiGraph();
	TGraphErrors* graph_p[6];
	for(int i = 0; i < 6; i++){
		graph_p[i] = GetTGraphErrors(R[i], P[i]);
		graph_p[i]->SetLineColor(1+i);
		multi2->Add(graph_p[i]);
	}

	c1->Divide(2, 1);
	c1->cd(1);
	multi1->Draw("APC");
	c1->cd(2);
	multi2->Draw("APC");
	*/

	TPad** pads = new TPad*[6]; 
	TMultiGraph** multi = new TMultiGraph*[6];
	TGraphErrors** graph_y = new TGraphErrors*[3];
	TGraphErrors** graph_p = new TGraphErrors*[3];
	TLegend** leg = new TLegend*[6];

	for(int i = 0; i < 3; i++){
		pads[2*i] = new TPad("", "", 0.0, 1-1./3*i, 0.5, 2./3-1./3*i, 0, 0, 0);
		pads[2*i]->cd();
		pads[2*i]->SetGrid();
		pads[2*i]->SetRightMargin(0.01);
		pads[2*i]->SetLeftMargin(0.11);
		pads[2*i]->SetBottomMargin(0.11);
		pads[2*i]->SetTopMargin(0.01);

		multi[2*i] = new TMultiGraph();
		graph_y[i] = GetTGraphErrors(X[i], Y[i]);

		graph_y[i]->SetLineColor(2+i);
		graph_y[i]->SetLineWidth(2);

		multi[2*i]->Add(graph_y[i]);

		multi[2*i]->SetTitle(";x;y");

		multi[2*i]->GetXaxis()->SetLabelSize(0.05);
		multi[2*i]->GetXaxis()->SetTitleSize(0.08);
		multi[2*i]->GetXaxis()->SetTitleOffset(0.6);
		multi[2*i]->GetYaxis()->SetLabelSize(0.05);
		multi[2*i]->GetYaxis()->SetTitleSize(0.08);
		multi[2*i]->GetYaxis()->SetTitleOffset(0.4);

		multi[2*i]->Draw("AC");

		//=============================================================

		pads[2*i+1] = new TPad("", "", 0.5, 1-1./3*i, 1.0, 2./3-1./3*i, 0, 0, 0);
		pads[2*i+1]->cd();
		pads[2*i+1]->SetGrid();
		pads[2*i+1]->SetRightMargin(0.01);
		pads[2*i+1]->SetLeftMargin(0.11);
		pads[2*i+1]->SetBottomMargin(0.11);
		pads[2*i+1]->SetTopMargin(0.01);

		multi[2*i+1] = new TMultiGraph();
		graph_p[i] = GetTGraphErrors(PX[i], PY[i]);

		graph_p[i]->SetLineColor(2+i);
		graph_p[i]->SetLineWidth(2);

		multi[2*i+1]->Add(graph_p[i]);

		multi[2*i+1]->SetTitle(";p_{x};p_{y}");
  
		multi[2*i+1]->GetXaxis()->SetLabelSize(0.05);
		multi[2*i+1]->GetXaxis()->SetTitleSize(0.08);
		multi[2*i+1]->GetXaxis()->SetTitleOffset(0.55);
		multi[2*i+1]->GetYaxis()->SetLabelSize(0.05);
		multi[2*i+1]->GetYaxis()->SetTitleSize(0.08);
		multi[2*i+1]->GetYaxis()->SetTitleOffset(0.4);

		multi[2*i+1]->GetXaxis()->SetLimits(-1.2, 3);
		multi[2*i+1]->GetYaxis()->SetRangeUser(-1.1, 1.1);

		multi[2*i+1]->Draw("AC");

		leg[2*i] = new TLegend(0.65, 0.45, 0.98, 0.65);
		leg[2*i]->SetTextSize(0.1);
		leg[2*i]->AddEntry(graph_p[i]  , (string("p_{x}(t=0) = ")+to_string(i-1)).c_str());
		leg[2*i]->DrawClone("SAME");


	}

	c1->cd();
	for(int i = 0; i < 6; i++) pads[i]->Draw();

	c1->SaveAs("Plot.png");

}