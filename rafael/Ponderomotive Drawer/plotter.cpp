
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



	DataSet* RM = new DataSet[3];
	DataSet* YM = new DataSet[3];
	DataSet* PM = new DataSet[3];

	TGraph** graph_ym = new TGraph*[3];
	TGraph** graph_pm = new TGraph*[3];

	for(int i = 0; i < 3; i++){
		RM[i] = (R[2*i] + R[2*i+1]) / Var(2);
		YM[i] = (Y[2*i] + Y[2*i+1]) / Var(2);
		PM[i] = (P[2*i] + P[2*i+1]) / Var(2);
	}







	TCanvas* c1 = new TCanvas("c1", "", 1500, 1000);
	
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

		graph_ym[i] = GetTGraphErrors(RM[i], YM[i]);
		graph_ym[i]->SetLineColor(3);
		graph_ym[i]->SetLineWidth(3);
		multi[2*i]->Add(graph_ym[i]);

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

		multi[2*i]->Draw("AC");

		if(i != 2) leg[2*i] = new TLegend(0.1, 0.7, 0.4, 0.92);
		else leg[2*i] = new TLegend(0.1, 0.76, 0.4, 0.95);
		leg[2*i]->AddEntry(graph_y[2*i]  , (string("E_{0} = 10  & w_{0} = ")+to_string(i*5+5)).c_str());
		leg[2*i]->AddEntry(graph_y[2*i+1], (string("E_{0} = -10 & w_{0} = ")+to_string(i*5+5)).c_str());
		leg[2*i]->AddEntry(graph_ym[i], "Average");
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

		graph_pm[i] = GetTGraphErrors(RM[i], PM[i]);		
		graph_pm[i]->SetLineWidth(3);
		graph_pm[i]->SetLineColor(3);
		multi[2*i+1]->Add(graph_pm[i]);

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

		multi[2*i+1]->Draw("AC");

	}

	c1->cd();
	for(int i = 0; i < 6; i++) pads[i]->Draw();

	c1->SaveAs("Plot.png");
	
	/*

	int r_min = -20;
	int n_data = 41;

	DataSet* X = new DataSet[n_data];
	DataSet* Y = new DataSet[n_data];
	DataSet* PX = new DataSet[n_data];
	DataSet* PY = new DataSet[n_data];

	int n_cols, n_points;
	double** values;

	for(int i = 0; i < n_data; i++){
		values = ReadFile(("Movement_r"+to_string(r_min+i)+".000000.txt").c_str(), &n_cols, &n_points, true);
		X[i] = DataSet(n_points, values[1]);
		Y[i] = DataSet(n_points, values[2]);
		PX[i] = DataSet(n_points, values[3]);
		PY[i] = DataSet(n_points, values[4]);
	}

	TCanvas* c1 = new TCanvas("c1", "", 1500, 500);
	TMultiGraph* multi = new TMultiGraph();
	TGraph** graph = new TGraph*[n_data];
	for(int i = 0; i < n_data; i++){
		graph[i] = GetTGraph(X[i], Y[i]);
		graph[i]->SetLineColor(i%3+2);
		graph[i]->SetLineWidth(2);
		multi->Add(graph[i]);
	}

	multi->SetTitle(";x;y");

  	multi->GetXaxis()->SetLabelSize(0.06);
	multi->GetXaxis()->SetTitleSize(0.08);
	multi->GetXaxis()->SetTitleOffset(0.2);

	multi->GetYaxis()->SetLabelSize(0.06);
	multi->GetYaxis()->SetTitleSize(0.08);
	multi->GetYaxis()->SetTitleOffset(0.2);
	multi->GetYaxis()->CenterTitle(true);

	multi->GetXaxis()->SetLimits(0, 75);
	multi->GetYaxis()->SetRangeUser(-85, 85);

	c1->cd();
	c1->SetGrid();
	c1->SetRightMargin(0.01);
	c1->SetLeftMargin(0.04);
	c1->SetBottomMargin(0.08);
	c1->SetTopMargin(0.01);
	multi->Draw("AC");

	c1->SaveAs("PlotTrajectories.png");
*/
	// E0 = 10 w0 = 10 lambda = 1 t = 150

}