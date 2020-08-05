
#include "DataAnalysis.cpp"

void plotter(){

	int n_cols, n_points;
	double** values;

	DataSet* R = new DataSet[6];
	DataSet* Y = new DataSet[6];
	DataSet* P = new DataSet[6];

	double E0[6] = {1, -1, 1, -1, 1, -1};
	double w0[6] = {5, 5, 10, 10, 15, 15};

	for(int i = 1; i <= 6; i++){
		values = ReadFile((string("Data_Ponderomotive_v")+to_string(i)+string(".txt")).c_str(), &n_cols, &n_points, true);
		R[i-1] = DataSet(n_points, values[0], 0.);
		Y[i-1] = DataSet(n_points, values[1], 0.);
		P[i-1] = DataSet(n_points, values[2], 0.);
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
	TGraphErrors** graph_y = new TGraphErrors*[6];
	TGraphErrors** graph_p = new TGraphErrors*[6];
	TLegend** leg = new TLegend*[6];

	for(int i = 0; i < 3; i++){
		if(i != 2) pads[2*i] = new TPad("", "", 0.0, 1-0.3*i, 0.5, 0.7-0.3*i, 0, 0, 0);
		else pads[2*i] = new TPad("", "", 0.0, 1-0.3*i, 0.5, 0, 0, 0, 0);
		pads[2*i]->cd();
		pads[2*i]->SetGrid();
		pads[2*i]->SetRightMargin(0.01);
		pads[2*i]->SetLeftMargin(0.08);
		if(i != 2) pads[2*i]->SetBottomMargin(0.0);
		else pads[2*i]->SetBottomMargin(0.25);
		pads[2*i]->SetTopMargin(0.01);

		multi[2*i] = new TMultiGraph();
		graph_y[2*i] = GetTGraphErrors(R[2*i], Y[2*i]);
		graph_y[2*i+1] = GetTGraphErrors(R[2*i+1], Y[2*i+1]);

		graph_y[2*i]->SetLineColor(2);
		graph_y[2*i+1]->SetLineColor(4);

		multi[2*i]->Add(graph_y[2*i]);
		multi[2*i]->Add(graph_y[2*i+1]);

		multi[2*i]->SetTitle(";y_{0};y_{f}");

		multi[2*i]->GetXaxis()->SetLabelSize(0.06);
		multi[2*i]->GetXaxis()->SetTitleSize(0.08);
		//multi[2*i]->GetXaxis()->SetTitleOffset(-0.12);
		multi[2*i]->GetXaxis()->CenterTitle(true);

		multi[2*i]->GetYaxis()->SetLabelSize(0.06);
		multi[2*i]->GetYaxis()->SetTitleSize(0.08);
		multi[2*i]->GetYaxis()->SetTitleOffset(0.3);
		multi[2*i]->GetYaxis()->CenterTitle(true);

		multi[2*i]->Draw("ACP");

		if(i != 2) leg[2*i] = new TLegend(0.1, 0.75, 0.4, 0.92);
		else leg[2*i] = new TLegend(0.1, 0.81, 0.4, 0.95);
		leg[2*i]->AddEntry(graph_y[2*i]  , (string("E_{0} = 1  & w_{0} = ")+to_string(i*5+5)).c_str());
		leg[2*i]->AddEntry(graph_y[2*i+1], (string("E_{0} = -1 & w_{0} = ")+to_string(i*5+5)).c_str());
		leg[2*i]->DrawClone("SAME");

		//=============================================================

		if(i != 2) pads[2*i+1] = new TPad("", "", 0.5, 1-0.3*i, 1.0, 0.7-0.3*i, 0, 0, 0);
		else pads[2*i+1] = new TPad("", "", 0.5, 1-0.3*i, 1.0, 0, 0, 0, 0);		
		pads[2*i+1]->cd();
		pads[2*i+1]->SetGrid();
		pads[2*i+1]->SetRightMargin(0.01);
		pads[2*i+1]->SetLeftMargin(0.08);
		if(i != 2) pads[2*i+1]->SetBottomMargin(0.0);
		else pads[2*i+1]->SetBottomMargin(0.25);
		pads[2*i+1]->SetTopMargin(0.01);

		multi[2*i+1] = new TMultiGraph();
		graph_p[2*i] = GetTGraphErrors(R[2*i], P[2*i]);
		graph_p[2*i+1] = GetTGraphErrors(R[2*i+1], P[2*i+1]);

		graph_p[2*i]->SetLineColor(2);
		graph_p[2*i+1]->SetLineColor(4);

		multi[2*i+1]->Add(graph_p[2*i]);
		multi[2*i+1]->Add(graph_p[2*i+1]);

		multi[2*i+1]->SetTitle(";y_{0};            p_{y}^{max}");
  
		multi[2*i+1]->GetXaxis()->SetLabelSize(0.06);
		multi[2*i+1]->GetXaxis()->SetTitleSize(0.08);
		//multi[2*i+1]->GetXaxis()->SetTitleOffset(-0.12);
		multi[2*i+1]->GetXaxis()->CenterTitle(true);

		multi[2*i+1]->GetYaxis()->SetLabelSize(0.06);
		multi[2*i+1]->GetYaxis()->SetTitleSize(0.08);
		multi[2*i+1]->GetYaxis()->SetTitleOffset(0.35);
		multi[2*i+1]->GetYaxis()->CenterTitle(true);

		multi[2*i+1]->Draw("ACP");

	}

	c1->cd();
	for(int i = 0; i < 6; i++) pads[i]->Draw();

	c1->SaveAs("Plot.png");

}