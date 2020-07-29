
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

	trajetoria->SetTitle(";x;y");
	trajetoria->SetMarkerColor(kRed);

	momentos->SetTitle(";px;py");
	momentos->SetMarkerColor(kRed);


	TCanvas* c1 = new TCanvas("c1", "", 1200, 600);
	c1->Divide(2,1);

	c1->cd(1);
	trajetoria->Draw("AP");

	c1->cd(2);
	momentos->Draw("AP");

}