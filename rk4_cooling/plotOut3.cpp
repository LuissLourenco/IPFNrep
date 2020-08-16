#include "DataAnalysis.cpp"

void plotOut3(){

	int n_cols, n_points;
	double** values = ReadFile("Out3.txt", &n_cols, &n_points, false);

	TGraph* graph_y = GetTGraph(DataSet(n_points, values[1]), DataSet(n_points, values[2]));
	TGraph* graph_p = GetTGraph(DataSet(n_points, values[4]), DataSet(n_points, values[5]));
	
	graph_y->SetTitle(";x;y");
	graph_p->SetTitle(";px;py");
	//graph_y->SetMarkerColor(kAzure+9);
	//graph_p->SetMarkerColor(kAzure+9);
	graph_y->SetMarkerColor(kRed);
	graph_p->SetMarkerColor(kRed);

	graph_y->GetYaxis()->SetMaxDigits(3);
	graph_y->GetYaxis()->SetTitleOffset(1);

	TCanvas* c1 = new TCanvas("c1", "", 1500, 1000);
	c1->Divide(2, 1);	
	c1->cd(1);
	graph_y->Draw("AP");
	c1->cd(2);
	graph_p->Draw("AP");

	c1->SaveAs("Plot.png");

}