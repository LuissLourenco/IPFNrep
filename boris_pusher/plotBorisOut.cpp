
#include "Read.cpp"

void plotBorisOut(){

	// t x y px py pz gama
	vector<vector<double>> sol = ReadFile2Vec("BorisOut.txt");
	int n = sol.size();

	TGraph2D* trajetoria = new TGraph2D(n);
	for (int i=0; i<n; i++){
		trajetoria->SetPoint(i, sol[i][1], sol[i][2], sol[i][3]);
	}

	TGraph2D* momentos = new TGraph2D(n);
	for (int i=0; i<n; i++){
		momentos->SetPoint(i, sol[i][4], sol[i][5], sol[i][6]);
	}	

	trajetoria->SetTitle(";x;y;z");
	trajetoria->SetMarkerColor(kRed);

	momentos->SetTitle(";px;py;pz");
	momentos->SetMarkerColor(kRed);


	TCanvas* c1 = new TCanvas("c1", "", 1200, 600);
	c1->Divide(2,1);

	c1->cd(1);
	trajetoria->Draw("P");

	c1->cd(2);
	momentos->Draw("P");

}