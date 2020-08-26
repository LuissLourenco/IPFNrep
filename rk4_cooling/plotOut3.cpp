#include "DataAnalysis.cpp"

void plotOut3(){


	FILE *foo;
	char trash[128];
	double x01,x02,x03,p01,p02,p03, kdamp, T;
	long long int N; int pri;
	double dx, dy; int wave_type;
	double tfwhm, stable, Eo, delta, w0, lambda, n, eta;
	int l,p;
	foo=fopen("InputToBatch.txt","r");
	fscanf(foo,"%s %lf %lf %lf %lf %lf %lf %lf %lf %lli %i %lf %lf %i %lf %lf %lf %lf %lf %lf %lf %lf %i %i", 
				trash, &x01, &x02, &x03, &p01, &p02, &p03, &kdamp, &T, &N, &pri, &dx, &dy, 
				&wave_type, &tfwhm, &stable, &Eo, &delta, &w0, &lambda, &n, &eta, &l, &p);
	fclose(foo);





	int n_cols, n_points;
	double** values = ReadFile("Out3.txt", &n_cols, &n_points, false);

	TGraph* graph_y = GetTGraph(DataSet(n_points, values[1]), DataSet(n_points, values[2]));
	TGraph* graph_p = GetTGraph(DataSet(n_points, values[4]), DataSet(n_points, values[5]));
		
	TGraph2D* graph_y2 = new TGraph2D(n_points,values[1],values[2],values[3]);
	TGraph2D* graph_p2 = new TGraph2D(n_points,values[4],values[5],values[6]);

	graph_y->SetTitle(";x;y");
	graph_p->SetTitle(";px;py");
	//graph_y->SetMarkerColor(kAzure+9);
	//graph_p->SetMarkerColor(kAzure+9);
	graph_y->SetMarkerColor(kRed);
	graph_p->SetMarkerColor(kRed);

	graph_y->GetYaxis()->SetMaxDigits(3);
	graph_y->GetYaxis()->SetTitleOffset(1);
	graph_p->GetYaxis()->SetMaxDigits(3);
	graph_p->GetYaxis()->SetTitleOffset(1);

	graph_y2->SetTitle(";x;y;z");
	graph_p2->SetTitle(";px;py;pz");
	graph_y2->SetMarkerColor(kRed);
	graph_p2->SetMarkerColor(kRed);

	graph_y2->GetZaxis()->SetRange(-1,1);
	graph_p2->GetZaxis()->SetRangeUser(-1,1);

	TCanvas* c1 = new TCanvas("c1", (string("wave_type=")+to_string(wave_type)).c_str(), 1500, 1000);
	c1->Divide(2, 1);	
	c1->cd(1);
	graph_y->Draw("AP");
	//graph_y2->Draw("P");
	c1->cd(2);
	graph_p->Draw("AP");
	//graph_p2->Draw("P");

	c1->SaveAs("Plot.png");
}