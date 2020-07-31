
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
		trajetoria_teorica->SetPoint(i, Eo*Eo/4*(sol[i][0]-sol[i][1]+(delta*delta-1/2)*cos(2*(sol[i][0]-sol[i][1]))),
										delta*Eo*sin(sol[i][0]-sol[i][1]),
										-sqrt(1-delta*delta)*Eo*cos(sol[i][0]-sol[i][1]));
	}

	TGraph2D* momentos = new TGraph2D(n);
	for (int i=0; i<n; i++){
		momentos->SetPoint(i, sol[i][4], sol[i][5], sol[i][6]);
	}	

	trajetoria->SetTitle(";x;y;z");
	trajetoria->SetMarkerColor(kRed);
	trajetoria_teorica->SetMarkerColor(kBlue);

	momentos->SetTitle(";px;py;pz");
	momentos->SetMarkerColor(kRed);

	TCanvas* c1 = new TCanvas("c1", "", 2400, 1200);
	c1->Divide(2,1);

	c1->cd(1);
	trajetoria->Draw("P");
	trajetoria_teorica->Draw("SAME");

	c1->cd(2);
	momentos->Draw("P");


	c1->SaveAs("Plot.png");

}