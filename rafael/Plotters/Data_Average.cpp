#include "DataAnalysis.cpp"

void Data_Average(){

	int n_cols, n_points;
	double** values = ReadFile("Data_Average.txt", &n_cols, &n_points, true);
	DataSet X(n_points, values[0]);
	DataSet Y(n_points, values[1]);

	TCanvas* c1 = new TCanvas("c1", "", 1200, 800);
	TGraph* graph = GetTGraph(X, Y);
	graph->SetMarkerColor(1);
	graph->SetMarkerStyle(24);
	graph->SetMarkerSize(1.5);

	graph->SetTitle(";a_{0};p_{x}^{ave}     ");

  	graph->GetXaxis()->SetLabelSize(0.06);
	graph->GetXaxis()->SetTitleSize(0.08);
	graph->GetXaxis()->SetTitleOffset(0.35);

	graph->GetYaxis()->SetLabelSize(0.06);
	graph->GetYaxis()->SetTitleSize(0.08);
	graph->GetYaxis()->SetTitleOffset(0.5);

	graph->GetXaxis()->SetLimits(0, 10.9);
	graph->GetYaxis()->SetRangeUser(-3.6, 1.1);

	TF1* fit_linear = new TF1("fit_linear", "pol1", 1, 10);
	graph->Fit(fit_linear, "0QR");
	fit_linear->SetLineColor(4);
	fit_linear->SetLineWidth(3);

	TF1* fit_quadratic = new TF1("fit_quadratic", "pol2", 0, 1);
	graph->Fit(fit_quadratic, "0QR");
	fit_quadratic->SetLineColor(2);
	fit_quadratic->SetLineWidth(3);

	c1->cd();
	c1->SetGrid();
	c1->SetRightMargin(0.01);
	c1->SetLeftMargin(0.1);
	c1->SetBottomMargin(0.08);
	c1->SetTopMargin(0.01);
	graph->Draw("AP");
	fit_linear->Draw("SAME");
	fit_quadratic->Draw("SAME");

	TLegend* leg = new TLegend(0.11, 0.8, 0.98, 0.98);
	leg->AddEntry(graph, "Data Simulated");
	leg->AddEntry(fit_quadratic, ("y = "+ to_string(fit_quadratic->GetParameter(0)) + " + " + to_string(fit_quadratic->GetParameter(1)) + " #times x + " + to_string(fit_quadratic->GetParameter(2)) + " #times x^{2}").c_str());
	leg->AddEntry(fit_linear, ("y = "+ to_string(fit_linear->GetParameter(0)) + " + " + to_string(fit_linear->GetParameter(1)) + " #times x").c_str());
	leg->DrawClone("SAME");

	c1->SaveAs("Data_Average.png");

	//AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA3333333333333333333333333333333333333
	//AAAAAAAAAAAAAAAA

}