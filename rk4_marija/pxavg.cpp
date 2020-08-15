
#include "trial.c"
#include "Read.cpp"
#include <thread>

double simulate_px(double a0){

	double px;

	ofstream itb;
	itb.open("InputToBatch.txt");
	itb<<"trash|x01|x02|x03|p01|p02|p03|T|N|pri|E0|w0"<<endl;
	itb<<"0."<<endl;
	itb<<"0."<<endl;
	itb<<"0."<<endl;
	itb<<a0*a0*0.5<<endl;
	itb<<-a0<<endl;
	itb<<"0."<<endl;
	itb<<"1000."<<endl;
	itb<<"1000000"<<endl;
	itb<<"10"<<endl;
	itb<<a0<<endl;
	itb<<"10."<<endl;
	itb.close();

	main();

	auto movimento = ReadFile2Vec("Out0.txt");
	int n = movimento.size();
	px = 0;
	for(int i=0 ; i<n; i++)
		px += movimento[i][4]/(double)n;

	return px;
}


void pxavg(){


	vector<double> a0_vec;
	vector<double> px_vec;
	
	double da0 = .1;
	double a0 = .1;

	bool simulate = false;

	if (simulate){

		while(a0<=10.){
			
			a0_vec.push_back(a0);
			px_vec.push_back(simulate_px(a0));

			a0+=da0;
		}


		ofstream out;
		out.open("pxavg_results.txt");
		int n = a0_vec.size();
		for(int i=0; i<n; i++)
			out<<a0_vec[i]<<"\t"<<px_vec[i]<<endl;
		out.close();
	}


	if (!simulate){
		a0_vec.clear();
		px_vec.clear();
		auto resultados = ReadFile2Vec("pxavg_results.txt");
		int nr = resultados.size();
		for(int i=0; i<nr; i++){
			a0_vec.push_back(resultados[i][0]);
			px_vec.push_back(resultados[i][1]);
		}
	}


	auto g1 = new TGraph(n);
	for(int i=0; i<a0_vec.size(); i++)
		g1->SetPoint(i, a0_vec[i], px_vec[i]);
	g1->SetTitle(";a_{0};#bar{p_{x}}");
	g1->SetMarkerStyle(8);
	g1->SetMarkerSize(1);
	g1->SetMarkerColor(kViolet);
	g1->GetXaxis()->SetLabelSize(0.05);
	g1->GetYaxis()->SetLabelSize(0.05);
	g1->GetXaxis()->SetTitleSize(0.05);
	g1->GetYaxis()->SetTitleSize(0.05);
	g1->GetXaxis()->SetTitleOffset(0.7);
	g1->GetYaxis()->SetTitleOffset(0.7);

	auto f1 = new TF1("f1", "[0]*x*x", 0, 11);
	f1->SetParameter(0,0.25);

	auto f2 = new TF1("f2", "[0]*x*x", 0, 11);
	g1->Fit(f2,"r");
	f2->SetLineColor(kBlue);

	auto c1 = new TCanvas("c1","",1500,1000);
	g1->Draw("AP");
	f1->Draw("SAME");
	f2->Draw("SAME");
	auto l1 = new TLegend(.75, .15, .87, .27);
	l1->AddEntry(g1, "simulated");
	l1->AddEntry(f1, "a_{0}^{2}/4");
	l1->AddEntry(f2, "Fit: p_{0}*x*x");
	l1->Draw();

	gPad->WaitPrimitive();
	c1->SaveAs("pxavg.png");

	delete g1;
	delete f1;
	delete c1;
	delete l1;

}