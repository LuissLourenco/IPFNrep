
#include "Read.cpp"

void plotOut0(){

	double x01, x02, x03, p01, p02, p03, Eo, delta, T;
	long long int N;
	int pri;
	char trash[64];
	FILE *foo;
	foo=fopen("InputToBatch.txt","r");
	fscanf(foo,"%s %lf %lf %lf %lf %lf %lf %lf %lli %i ", 
				trash, &x01, &x02, &x03, &p01, &p02, &p03, &T, &N, &pri);
	fclose(foo);


	// t x y z px py pz gama
	vector<vector<double>> sol = ReadFile2Vec("Out0.txt");
	int n = sol.size();

	delta = 0.7071;
	Eo = 10;


	double avgx=0;
	TGraph2D* trajetoria = new TGraph2D(n);
	for (int i=0; i<n; i++){
		avgx +=sol[i][1]/(double)n;
		trajetoria->SetPoint(i, sol[i][1], sol[i][2], sol[i][3]);
	}
	cout<<avgx<<endl;
	TGraph2D* momentos = new TGraph2D(n);
	for (int i=0; i<n; i++){
		momentos->SetPoint(i, sol[i][4], sol[i][5], sol[i][6]);
	}		

	//LAB FRAME
	/*TGraph2D* trajetoria_teorica = new TGraph2D(n);
	for (int i=0; i<n; i++){
		trajetoria_teorica->SetPoint(i, Eo*Eo/4.*(sol[i][0]-sol[i][1]+(delta*delta-1./2.)*sin(2.*(sol[i][0]-sol[i][1]))),
										delta*Eo*sin(sol[i][0]-sol[i][1]),
										-sqrt(1-delta*delta)*Eo*cos(sol[i][0]-sol[i][1]));
	}

	TGraph2D* momentos_teoricos = new TGraph2D(n);
	for (int i=0; i<n; i++){
		momentos_teoricos->SetPoint(i, Eo*Eo/4.*(1.+(2.*delta*delta-1.)*cos(2.*(sol[i][0]-sol[i][1]))),
										delta*Eo*cos(sol[i][0]-sol[i][1]),
										sqrt(1-delta*delta)*Eo*sin(sol[i][0]-sol[i][1]));
	}*/
	//REST FRAME
	TPolyLine3D* trajetoria_teorica = new TPolyLine3D(n);
	double gama0 = sqrt(1.+Eo*Eo/2.);
	double q = Eo/2./gama0;
	for (int i=0; i<n; i++){
		trajetoria_teorica->SetPoint(i, avgx+(delta*delta-0.5)*q*q*sin(2*(sol[i][0]-sol[i][1])),
										2.*delta*q*sin(sol[i][0]-sol[i][1]),
										-2*sqrt(1-delta*delta)*q*cos(sol[i][0]-sol[i][1]));
	}
	TGraph2D* momentos_teoricos = new TGraph2D(n);
	for (int i=0; i<n; i++){
		momentos_teoricos->SetPoint(i, (2*delta*delta-1)*Eo*Eo/4./gama0*cos(2*(sol[i][0]-sol[i][1])),
										delta*Eo*cos(sol[i][0]-sol[i][1]),
										sqrt(1-delta*delta)*Eo*sin(sol[i][0]-sol[i][1]));
	}

	/*
	double avgx=0;

	TGraph* trajetoria = new TGraph(n);
	for (int i=0; i<n; i++){
		avgx += sol[i][1]/(double)n;
		trajetoria->SetPoint(i, sol[i][1], sol[i][2]);
	}


	TGraph* momentos = new TGraph(n);
	for (int i=0; i<n; i++){
		momentos->SetPoint(i, sol[i][4], sol[i][5]);
	}	

	//LAB FRAME 
	TGraph* trajetoria_teorica = new TGraph(n);
	for (int i=0; i<n; i++){
		trajetoria_teorica->SetPoint(i, Eo*Eo/4.*(sol[i][0]-sol[i][1]+(delta*delta-1./2.)*sin(2.*(sol[i][0]-sol[i][1]))),
										delta*Eo*sin(sol[i][0]-sol[i][1])
										);
	}
	TGraph* momentos_teoricos = new TGraph(n);
	for (int i=0; i<n; i++){
		momentos_teoricos->SetPoint(i, Eo*Eo/4.*(1.+(2.*delta*delta-1.)*cos(2.*(sol[i][0]-sol[i][1]))),
										delta*Eo*cos(sol[i][0]-sol[i][1])
										);
	}

	//REST FRAME
	TGraph* trajetoria_teorica = new TGraph(n);
	double gama0 = sqrt(1.+Eo*Eo/2.);
	double q = Eo/2./gama0;
	for (int i=0; i<n; i++){
		trajetoria_teorica->SetPoint(i, (delta*delta-0.5)*q*q*sin(2*(sol[i][0]-sol[i][1])),
										2*delta*q*sin(sol[i][0]-sol[i][1])
										);
	}
	TGraph* momentos_teoricos = new TGraph(n);
	for (int i=0; i<n; i++){
		momentos_teoricos->SetPoint(i, (2*delta*delta-1)*Eo*Eo/4./gama0*cos(2*(sol[i][0]-sol[i][1])),
										delta*Eo*cos(sol[i][0]-sol[i][1])
										);
	}
	*/

	trajetoria->SetTitle(";x;y;z");
	trajetoria->SetMarkerColor(kPink);
	//trajetoria_teorica->SetMarkerColor(kAzure+9);
	//trajetoria_teorica->SetLineColor(kAzure+9);
	trajetoria_teorica->SetLineColorAlpha(kAzure+9,1);
	trajetoria_teorica->SetLineWidth(4);

	momentos->SetTitle(";px;py;pz");
	momentos->SetMarkerColor(kPink);
	momentos_teoricos->SetMarkerColor(kAzure+9);
	momentos_teoricos->SetLineColor(kAzure+9);

	trajetoria->GetXaxis()->SetMaxDigits(3);
	//trajetoria->GetXaxis()->SetRangeUser(avgx-3, avgx+3);
	//trajetoria_teorica->SetMarkerStyle(8);
	//trajetoria_teorica->SetMarkerSize(0.5);
	momentos_teoricos->SetMarkerStyle(8);
	momentos_teoricos->SetMarkerSize(0.5);
	trajetoria->GetXaxis()->SetTitleOffset(1.5);
	trajetoria->GetYaxis()->SetTitleOffset(0.5);

	/*
	trajetoria->SetMinimum(-1.);
	trajetoria->SetMaximum(1.);
	trajetoria_teorica->SetMinimum(-1.);
	trajetoria_teorica->SetMaximum(1.);

	momentos->SetMinimum(-1.);
	momentos->SetMaximum(1.);
	momentos_teoricos->SetMinimum(-1.);
	momentos_teoricos->SetMaximum(1.);
	*/


	TCanvas* c1 = new TCanvas("c1", "", 1800, 600);
	c1->Divide(3,1);

	c1->cd(1)->SetPad(0, 0, 0.45, 1);
	c1->cd(2)->SetPad(0.44, 0, 0.9, 1);
	c1->cd(3)->SetPad(0.9, 0, 1, 1);

	c1->cd(1);
	c1->cd(1)->SetGrid();
	trajetoria->Draw("P");
	trajetoria_teorica->Draw("SAME4");

	c1->cd(2);
	c1->cd(2)->SetGrid();
	momentos->Draw("P");
	momentos_teoricos->Draw("SAME P");



	TLatex* t[5];
	t[0] = new TLatex(.0, .8, "a_{0} = 10");
	t[1] = new TLatex(.0, .77, "delta = 0.7071");
	t[2] = new TLatex(.0, .74, "#vec{x}_{0}=0, #vec{p}_{0}=-3.5#hat{x}");
	t[3] = new TLatex(.0, .71, "t_{subida}=20pi");
	t[4] = new TLatex(.0, .68, "t_{stable}=300");

	for(int i=0; i<5; i++) t[i]->SetTextSize(0.1);

	c1->cd(3);
	for(int i=0; i<5; i++)
		t[i]->Draw();

	/*
	momentos->GetXaxis()->SetTitleOffset(1.75);
	momentos->GetYaxis()->SetTitleOffset(1.75);
	momentos->GetZaxis()->SetTitleOffset(1.55);

	trajetoria->GetXaxis()->SetTitleOffset(1.75);
	trajetoria->GetYaxis()->SetTitleOffset(1.85);
	trajetoria->GetZaxis()->SetTitleOffset(1.55);
	*/


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

	/*
	gPad->WaitPrimitive();

	delete trajetoria;
	delete trajetoria_teorica;
	delete momentos_teoricos;
	delete momentos;
	delete c1;
	*/
}