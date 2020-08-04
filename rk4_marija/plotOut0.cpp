
#include "Read.cpp"

void plotOut0(){

	double p01, p02, p03, Eo, delta, kdamp, k, T, xgrid, ygrid;
	long long int N;
	int pri;
	FILE *foo;
	foo=fopen("InputTotBatch.txt","r");
	fscanf(foo,"%lf %lf %lf %lf %lf %lf %lf %lf %lld %i %lf %lf ", &p01, &p02, &p03, &Eo, &delta, &kdamp,&k, &T, &N, &pri, &xgrid, &ygrid);

	// t x y z px py pz gama
	vector<vector<double>> sol = ReadFile2Vec("Out0.txt");
	int n = sol.size();

	TGraph2D* trajetoria = new TGraph2D(n);
	for (int i=0; i<n; i++){
		trajetoria->SetPoint(i, sol[i][1], sol[i][2], sol[i][3]);
	}

	TGraph2D* trajetoria_teorica = new TGraph2D(n);
	for (int i=0; i<n; i++){
		trajetoria_teorica->SetPoint(i, Eo*Eo/4.*(sol[i][0]-sol[i][1]+(delta*delta-1./2.)*sin(2.*(sol[i][0]-sol[i][1]))),
										delta*Eo*sin(sol[i][0]-sol[i][1]),
										-sqrt(1-delta*delta)*Eo*cos(sol[i][0]-sol[i][1]));
	}

	TGraph2D* momentos = new TGraph2D(n);
	for (int i=0; i<n; i++){
		momentos->SetPoint(i, sol[i][4], sol[i][5], sol[i][6]);
	}	

	TGraph2D* momentos_teoricos = new TGraph2D(n);
	for (int i=0; i<n; i++){
		momentos_teoricos->SetPoint(i, Eo*Eo/4.*(1.+(2.*delta*delta-1.)*cos(2.*(sol[i][0]-sol[i][1]))),
										delta*Eo*cos(sol[i][0]-sol[i][1]),
										sqrt(1-delta*delta)*Eo*sin(sol[i][0]-sol[i][1]));
	}

	trajetoria->SetTitle(";x;y;z");
	trajetoria->SetMarkerColor(kRed);
	trajetoria_teorica->SetMarkerColor(kBlue);

	momentos->SetTitle(";px;py;pz");
	momentos->SetMarkerColor(kRed);
	momentos_teoricos->SetMarkerColor(kBlue);



	trajetoria->SetMinimum(-1.);
	trajetoria->SetMaximum(1.);
	trajetoria_teorica->SetMinimum(-1.);
	trajetoria_teorica->SetMaximum(1.);

	momentos->SetMinimum(-1.);
	momentos->SetMaximum(1.);
	momentos_teoricos->SetMinimum(-1.);
	momentos_teoricos->SetMaximum(1.);



	TCanvas* c1 = new TCanvas("c1", "", 2400, 1200);
	c1->Divide(2,1);

	c1->cd(1);
	trajetoria->Draw("P");
	trajetoria_teorica->Draw("SAME");

	c1->cd(2);
	momentos->Draw("P");
	momentos_teoricos->Draw("SAME");

	momentos->GetXaxis()->SetTitleOffset(1.75);
	momentos->GetYaxis()->SetTitleOffset(1.75);
	momentos->GetZaxis()->SetTitleOffset(1.55);


	/*
	c1->Clear();
	c1->Divide(2,1);
	TGraph* graph_pos = new TGraph(n);
	for (int i=0; i<n; i++){
		graph_pos->SetPoint(i, sol[i][1], sol[i][3]);
	}
	TGraph* graph_vel = new TGraph(n);
	for (int i=0; i<n; i++){
		graph_vel->SetPoint(i, sol[i][4], sol[i][6]);
	}
	graph_pos->SetTitle("Trajectory;x;y");	
	graph_pos->GetYaxis()->SetTitleOffset(0.55);
	graph_vel->SetTitle("Momentum;p_{x};p_{y}");	
	graph_vel->GetYaxis()->SetTitleOffset(0.55);
	c1->cd(1);
	graph_pos->Draw("ACP");
	c1->cd(2);
	graph_vel->Draw("ACP");
	*/
	c1->SaveAs("Plot.png");

}