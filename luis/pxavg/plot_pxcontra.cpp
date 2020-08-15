
#include "Read.cpp"


void plot_pxcontra(){


	auto pontos = ReadFile2Vec("Data_Average.txt");
	int n = pontos.size();

	auto f1 = new TF1("f1", "-x*x/4/sqrt(1+x*x/2)", 0, 11);
	f1->SetLineColor(kBlue);

	auto g1 = new TGraph(n);
	for(int i=0; i<n; i++)
		g1->SetPoint(i, pontos[i][0], pontos[i][1]);
	g1->SetTitle(";a_{0};#bar{p_{x}}");
	g1->SetMarkerStyle(8);
	g1->SetMarkerSize(1);
	g1->SetMarkerColor(kViolet+6);
	g1->GetXaxis()->SetLabelSize(0.05);
	g1->GetYaxis()->SetLabelSize(0.05);
	g1->GetXaxis()->SetTitleSize(0.07);
	g1->GetYaxis()->SetTitleSize(0.07);
	g1->GetXaxis()->SetTitleOffset(0.5);
	g1->GetYaxis()->SetTitleOffset(0.5);

	auto pt1 = new TPaveText(7.9,-1.8,10.1,-0.7);
	pt1->AddText("p_{x} = #frac{a_{0}^{2}}{4#gamma_{0}}");
	((TText*)pt1->GetListOfLines()->Last())->SetTextSize(0.07);
	//auto t1 = new TLatex(8,-3,"p_{x} = #frac{a_{0}^{2}}{4#gamma_{0}}");
	//t1->SetTextSize(0.05);

	auto c1 = new TCanvas("c1", "", 1500, 1000);
	c1->SetGrid();

	g1->Draw("AP");
	f1->Draw("SAME");
	pt1->Draw();



	gPad->WaitPrimitive();
	c1->SaveAs("pxavg_contra.png");

	delete c1;
	delete f1;
	delete g1;

}


