
#include "DataAnalysis.cpp"

void Data_Ponderomotive_kdamp(){
	
	int n_cols, n_points;
	double** values;

	DataSet* R = new DataSet[6];
	DataSet* Y = new DataSet[6];
	DataSet* P = new DataSet[6];
	DataSet* E = new DataSet[6];

	double E0[6] = {1, -1, 1, -1, 1, -1};
	double w0[6] = {5, 5, 10, 10, 15, 15};

	for(int i = 1; i <= 6; i++){
		values = ReadFile((string("Data_Ponderomotive_kdamp_v")+to_string(i)+string(".txt")).c_str(), &n_cols, &n_points, true);
		R[i-1] = DataSet(n_points, values[0], 0.);
		Y[i-1] = DataSet(n_points, values[1], 0.) - R[i-1] + Var(25.);
		P[i-1] = DataSet(n_points, values[2], 0.);
		E[i-1] = DataSet(n_points, values[5], 0.);
	}

	TCanvas* c1 = new TCanvas("c1", "", 1500, 1000);
	
	TPad** pads = new TPad*[9]; 
	TMultiGraph** multi = new TMultiGraph*[9];
	TGraphErrors** graph_y = new TGraphErrors*[9];
	TGraphErrors** graph_p = new TGraphErrors*[9];
	TGraphErrors** graph_e = new TGraphErrors*[9];
	TLegend** leg = new TLegend*[3];

	for(int i = 0; i < 3; i++){
		if(i != 2) pads[3*i] = new TPad("", "", 0.0, 1-0.3*i, 1./3, 0.7-0.3*i, 0, 0, 0);
		else pads[3*i] = new TPad("", "", 0.0, 1-0.3*i, 1./3, 0, 0, 0, 0);
		pads[3*i]->cd();
		pads[3*i]->SetGrid();
		pads[3*i]->SetRightMargin(0.01);
		pads[3*i]->SetLeftMargin(0.1);
		if(i != 2) pads[3*i]->SetBottomMargin(0.0);
		else pads[3*i]->SetBottomMargin(0.25);
		pads[3*i]->SetTopMargin(0.01);

		multi[3*i] = new TMultiGraph();
		graph_y[3*i] = GetTGraphErrors(R[2*i], Y[2*i]);
		graph_y[3*i+1] = GetTGraphErrors(R[2*i+1], Y[2*i+1]);

		graph_y[3*i]->SetLineColor(2);
		graph_y[3*i+1]->SetLineColor(4);
		graph_y[3*i]->SetLineWidth(2);
		graph_y[3*i+1]->SetLineWidth(2);


		multi[3*i]->Add(graph_y[3*i]);
		multi[3*i]->Add(graph_y[3*i+1]);

		multi[3*i]->SetTitle(";y_{0};y_{f}");

		multi[3*i]->GetXaxis()->SetLabelSize(0.06);
		multi[3*i]->GetXaxis()->SetTitleSize(0.08);
		//multi[2*i]->GetXaxis()->SetTitleOffset(-0.12);
		multi[3*i]->GetXaxis()->CenterTitle(true);

		multi[3*i]->GetYaxis()->SetLabelSize(0.06);
		multi[3*i]->GetYaxis()->SetTitleSize(0.08);
		multi[3*i]->GetYaxis()->SetTitleOffset(0.5);
		multi[3*i]->GetYaxis()->CenterTitle(true);

		if(i == 2) multi[3*i]->GetYaxis()->SetLabelSize(0.044);
		if(i == 2) multi[3*i]->GetYaxis()->SetTitleSize(0.06);
		if(i == 2) multi[3*i]->GetYaxis()->SetTitleOffset(0.7);

		multi[3*i]->GetXaxis()->SetLimits(-24.9, 24.9);
		multi[3*i]->GetYaxis()->SetRangeUser(-0.0022, 0.0022);

		multi[3*i]->Draw("AC");

		if(i != 2) leg[i] = new TLegend(0.12, 0.7, 0.55, 0.92);
		else leg[i] = new TLegend(0.12, 0.76, 0.55, 0.95);
		if(i != 2) leg[i]->SetTextSize(0.051);
		else leg[i]->SetTextSize(0.04);
		leg[i]->AddEntry(graph_y[3*i]  , (string("w_{0} = ")+to_string(i*5+5) + " & k_{damp} = 0").c_str());
		leg[i]->AddEntry(graph_y[3*i+1], (string("w_{0} = ")+to_string(i*5+5) + " & k_{damp} = 1.18E-8").c_str());
		leg[i]->DrawClone("SAME");

		//=============================================================

		if(i != 2) pads[3*i+1] = new TPad("", "", 1./3, 1-0.3*i, 2./3, 0.7-0.3*i, 0, 0, 0);
		else pads[3*i+1] = new TPad("", "", 1./3, 1-0.3*i, 2./3, 0, 0, 0, 0);		
		pads[3*i+1]->cd();
		pads[3*i+1]->SetGrid();
		pads[3*i+1]->SetRightMargin(0.01);
		pads[3*i+1]->SetLeftMargin(0.1);
		if(i != 2) pads[3*i+1]->SetBottomMargin(0.0);
		else pads[3*i+1]->SetBottomMargin(0.25);
		pads[3*i+1]->SetTopMargin(0.01);

		multi[3*i+1] = new TMultiGraph();
		graph_p[3*i] = GetTGraphErrors(R[2*i], P[2*i]);
		graph_p[3*i+1] = GetTGraphErrors(R[2*i+1], P[2*i+1]);

		graph_p[3*i]->SetLineColor(2);
		graph_p[3*i+1]->SetLineColor(4);
		graph_p[3*i]->SetLineWidth(2);
		graph_p[3*i+1]->SetLineWidth(2);
		graph_p[3*i+1]->SetLineStyle(7);

		multi[3*i+1]->Add(graph_p[3*i]);
		multi[3*i+1]->Add(graph_p[3*i+1]);

		multi[3*i+1]->SetTitle(";y_{0};p_{y}^{max}");
  
		multi[3*i+1]->GetXaxis()->SetLabelSize(0.06);
		multi[3*i+1]->GetXaxis()->SetTitleSize(0.08);
		//multi[3*i+1]->GetXaxis()->SetTitleOffset(-0.12);
		multi[3*i+1]->GetXaxis()->CenterTitle(true);

		multi[3*i+1]->GetYaxis()->SetLabelSize(0.06);
		multi[3*i+1]->GetYaxis()->SetTitleSize(0.08);
		multi[3*i+1]->GetYaxis()->SetTitleOffset(0.5);
		multi[3*i+1]->GetYaxis()->CenterTitle(true);

		if(i == 2) multi[3*i+1]->GetYaxis()->SetLabelSize(0.044);
		if(i == 2) multi[3*i+1]->GetYaxis()->SetTitleSize(0.06);
		if(i == 2) multi[3*i+1]->GetYaxis()->SetTitleOffset(0.7);

		multi[3*i+1]->GetXaxis()->SetLimits(-24.9, 24.9);
		multi[3*i+1]->GetYaxis()->SetRangeUser(-0.2, 8.4);

		multi[3*i+1]->Draw("AC");

		//=============================================================

		if(i != 2) pads[3*i+2] = new TPad("", "", 2./3, 1-0.3*i, 1.0, 0.7-0.3*i, 0, 0, 0);
		else pads[3*i+2] = new TPad("", "", 2./3, 1-0.3*i, 1.0, 0, 0, 0, 0);		
		pads[3*i+2]->cd();
		pads[3*i+2]->SetGrid();
		pads[3*i+2]->SetRightMargin(0.01);
		pads[3*i+2]->SetLeftMargin(0.11);
		if(i != 2) pads[3*i+2]->SetBottomMargin(0.0);
		else pads[3*i+2]->SetBottomMargin(0.25);
		pads[3*i+2]->SetTopMargin(0.01);

		multi[3*i+2] = new TMultiGraph();
		graph_e[3*i] = GetTGraphErrors(R[2*i], E[2*i]);
		graph_e[3*i+1] = GetTGraphErrors(R[2*i+1], E[2*i+1]);

		graph_e[3*i]->SetLineColor(2);
		graph_e[3*i+1]->SetLineColor(4);
		graph_e[3*i]->SetLineWidth(2);
		graph_e[3*i+1]->SetLineWidth(2);

		multi[3*i+2]->Add(graph_e[3*i]);
		multi[3*i+2]->Add(graph_e[3*i+1]);

		multi[3*i+2]->SetTitle(";y_{0};#frac{E_{f}}{E_{i}}   ");
  
		multi[3*i+2]->GetXaxis()->SetLabelSize(0.06);
		multi[3*i+2]->GetXaxis()->SetTitleSize(0.08);
		//multi[3*i+2]->GetXaxis()->SetTitleOffset(-0.12);
		multi[3*i+2]->GetXaxis()->CenterTitle(true);

		multi[3*i+2]->GetYaxis()->SetLabelSize(0.06);
		multi[3*i+2]->GetYaxis()->SetTitleSize(0.08);
		multi[3*i+2]->GetYaxis()->SetTitleOffset(0.43);
		multi[3*i+2]->GetYaxis()->CenterTitle(true);

		if(i == 2) multi[3*i+2]->GetYaxis()->SetLabelSize(0.044);
		if(i == 2) multi[3*i+2]->GetYaxis()->SetTitleSize(0.05);
		if(i == 2) multi[3*i+2]->GetYaxis()->SetTitleOffset(0.7);

		multi[3*i+2]->GetXaxis()->SetLimits(-24.9, 24.9);
		multi[3*i+2]->GetYaxis()->SetRangeUser(1050, 2050);

		multi[3*i+2]->GetYaxis()->ChangeLabel(1,-1,-1,-1,-1,-1," ");
		multi[3*i+2]->GetYaxis()->ChangeLabel(2,-1,-1,-1,-1,-1,"60%");
		multi[3*i+2]->GetYaxis()->ChangeLabel(3,-1,-1,-1,-1,-1," ");
		multi[3*i+2]->GetYaxis()->ChangeLabel(4,-1,-1,-1,-1,-1,"70%");
		multi[3*i+2]->GetYaxis()->ChangeLabel(5,-1,-1,-1,-1,-1," ");
		multi[3*i+2]->GetYaxis()->ChangeLabel(6,-1,-1,-1,-1,-1,"80%");
		multi[3*i+2]->GetYaxis()->ChangeLabel(7,-1,-1,-1,-1,-1," ");
		multi[3*i+2]->GetYaxis()->ChangeLabel(8,-1,-1,-1,-1,-1,"90%");
		multi[3*i+2]->GetYaxis()->ChangeLabel(9,-1,-1,-1,-1,-1," ");
		multi[3*i+2]->GetYaxis()->ChangeLabel(10,-1,-1,-1,-1,-1,"100%");

		multi[3*i+2]->Draw("AC");

	}

	c1->cd();
	for(int i = 0; i < 9; i++) pads[i]->Draw();

	c1->SaveAs("Data_Ponderomotive_kdamp.png");
	
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