
#include "Read.cpp"

void plotOut0(){

	// t x y px py pz gama
	vector<vector<double>> sol = ReadFile2Vec("Out0.txt");
	int n = sol.size();

	TGraph* trajetoria = new TGraph(n);
	for (int i=0; i<n; i++){
		trajetoria->SetPoint(i, sol[i][1], sol[i][2]);
	}

	TGraph* momentos = new TGraph(n);
	for (int i=0; i<n; i++){
		momentos->SetPoint(i, sol[i][3], sol[i][4]);
	}	

	trajetoria->SetTitle("trajetoria;x;y");
	trajetoria->SetMarkerStyle(8);
	trajetoria->SetMarkerSize(1);
	trajetoria->SetMarkerColor(kRed);

	momentos->SetTitle("momentos;px;py");
	momentos->SetMarkerStyle(8);
	momentos->SetMarkerSize(1);
	momentos->SetMarkerColor(kRed);


	TCanvas* c1 = new TCanvas("c1", "", 1200, 600);
	trajetoria->Draw("AP");

	TCanvas* c2 = new TCanvas("c2", "", 1200, 600);
	momentos->Draw("AP");

}