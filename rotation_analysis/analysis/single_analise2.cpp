#include "../src/DataAnalysis.cpp"
#include <string>
#include <vector>

using namespace std;

void mood(string file_in, string plot_out, string file_out){

	int n_points, n_cols;
	double** values = ReadFile(file_in.c_str(), &n_cols, &n_points, true);

	//cout << "Saving to " << plot_out << "\n          " << file_out << endl;

	// prints t, x, y, z, px, py, pz, gamma, theta, p_theta
	DataSet T(n_points, values[0]);
	DataSet X(n_points, values[1]);
	DataSet Y(n_points, values[2]);
	DataSet Z(n_points, values[3]);
	DataSet PX(n_points, values[4]);
	DataSet PY(n_points, values[5]);
	DataSet PZ(n_points, values[6]);
	DataSet GAMMA(n_points, values[7]);

	if(sqrt(Y*Y+Z*Z).getMax().val() > 10){
		cout << "Radius diverges! Exiting..." << endl;
		return;
	}

	TCanvas* canvas = new TCanvas("", "", 2000, 1000);
	canvas->Divide(2, 2);
	TGraph* trajectory; double r_max = sqrt(Y[0]*Y[0]+Z[0]*Z[0]).val()/2;
	TGraph* sigma; TGraph* sigma2;
	TGraph* elipse; TF1* fit_elipse; TEllipse* elipse_draw;

	int N = 100;
	int start = 0;
	char name[64];

	double** dft_out = computeDft((T[1]-T[0]).val(), T.size(), Y.array(), 1000);
	int index = DataSet(1000, dft_out[2]).getMaxI();
	double osc_per = DataSet(1000, dft_out[0])[index].val();

	double raio_max;
	raio_max = sqrt(Y*Y+Z*Z).subDataSet((int)(n_points/2), n_points-1).getMax().val();

	N = osc_per/(T[1]-T[0]).val()*0.95;

	DataSet TSIGMA(1, 0.), SIGMA(1, 0.);
	DataSet PHI, R;
	while(start+N < n_points){

		trajectory = GetTGraph(Y.subDataSet(start, start+N), Z.subDataSet(start, start+N));

		PHI = atan2(Z.subDataSet(start, start+N), Y.subDataSet(start, start+N));
		R = sqrt(Y.subDataSet(start, start+N)*Y.subDataSet(start, start+N) + Z.subDataSet(start, start+N)*Z.subDataSet(start, start+N));

		elipse = GetTGraph(PHI, R);
		fit_elipse = new TF1("", "[0]*[1]/sqrt([0]*[0]*sin(x+[2])*sin(x+[2])+[1]*[1]*cos(x+[2])*cos(x+[2]))", -M_PI, M_PI);
		fit_elipse->SetParameter(0, r_max);
		fit_elipse->SetParameter(1, r_max);
		elipse->Fit(fit_elipse, "Q");
		if(fit_elipse->GetParameter(0) < fit_elipse->GetParameter(1)){
			double aux = fit_elipse->GetParameter(0);
			fit_elipse->SetParameter(0, fit_elipse->GetParameter(1));
			fit_elipse->SetParameter(1, aux);
			fit_elipse->SetParameter(2, fit_elipse->GetParameter(2)+M_PI/2);
		}
		while(fit_elipse->GetParameter(2) > M_PI) 
			fit_elipse->SetParameter(2, fit_elipse->GetParameter(2)-M_PI);
		while(fit_elipse->GetParameter(2) < 0) 
			fit_elipse->SetParameter(2, fit_elipse->GetParameter(2)+M_PI);


		TSIGMA.append(Var(T[start]));
		SIGMA.append(Var(fit_elipse->GetParameter(2)));
		
		start+=N;
	}	
	TSIGMA = TSIGMA.subDataSet(1, TSIGMA.size());
	SIGMA = SIGMA.subDataSet(1, SIGMA.size());

	// CHECK THE ELLIPSE FITTING ======================================
	canvas->cd(1);
	gPad->SetGrid();
	elipse->GetXaxis()->SetLimits(-M_PI, M_PI);
	elipse->GetYaxis()->SetRangeUser(1E-3, r_max);
	elipse->SetTitle("1) Fitando a Elipse;#theta;r");
	elipse->Draw("AP");
	fit_elipse->Draw("SAME");

	// CHECK THE ELLIPSE FITTING ======================================
	canvas->cd(2);
	gPad->SetGrid();
	trajectory->GetXaxis()->SetLimits(-r_max, r_max);
	trajectory->GetYaxis()->SetRangeUser(-r_max, r_max);
	trajectory->SetTitle("2) Checkando a Elipse;y;z");
	trajectory->Draw("ACP");
	elipse_draw = new TEllipse(0, 0, fit_elipse->GetParameter(0), fit_elipse->GetParameter(1), 0, 360, -fit_elipse->GetParameter(2)*180/M_PI);
	elipse_draw->SetFillStyle(0);
	elipse_draw->SetLineColor(3);
	elipse_draw->Draw("SAME");
	

	// SIGMA WITH PERIOD NOT NORMALIZED ======================================
	canvas->cd(3);
	gPad->SetGrid();
	sigma = GetTGraph(TSIGMA, SIGMA);
	sigma->SetTitle("3) Desenhando movimento;t;#sigma");
	sigma->Draw("ALP");	

	
	// OBTER SIGMA_dot ======================================
	

	//contagem aproximada de picos	

	//SIGMA.cut(SIGMA.size()/2);
	//TSIGMA.cut(TSIGMA.size()/2);

	/*
	int n_picos=0;
	bool aux=false;
	for(int i=0; i<SIGMA.size()-1; i++){
		if(!aux)
			if(SIGMA[i+1].val() > 0.8*M_PI){ aux = true; n_picos++; continue; }
		if(aux)
			if(SIGMA[i+1].val() < 0.8*M_PI){ aux = false; continue;}
	}
	cout << "nr picos: " << n_picos << endl;
	double Tguess = TSIGMA[-1].val()/(n_picos);
	cout << "Period guess: " << Tguess << endl;
	double p1guess = TSIGMA[SIGMA.getMaxI()].val() * M_PI / Tguess;
	cout << "T at max: " << TSIGMA[SIGMA.getMaxI()].val() << endl;
	cout << "P1 guess: " << p1guess << endl;


	canvas->cd(4);
	gPad->SetGrid();
	sigma2 = GetTGraph(TSIGMA, SIGMA);
	sigma2->SetTitle("4) Fit;t;#sigma");
	sigma2->GetYaxis()->SetRangeUser(-0.2,M_PI+0.2);
	sigma2->SetMarkerStyle(8);
	sigma2->SetMarkerSize(0.5);

	TF1* picardlindelof = new TF1("pp","pi/2+1*atan(1/tan(x*pi/[0] - [1]))",0,TSIGMA[-1].val());
	picardlindelof->SetParLimits(0, Tguess*0.95,Tguess*1.05);
	//picardlindelof->FixParameter(0, Tguess);
	picardlindelof->FixParameter(1, p1guess);

	picardlindelof->SetNpx(2000);
	
	sigma2->Draw("APL");
	sigma2->Fit(picardlindelof,"R");
	osc_per = 2 * picardlindelof->GetParameter(0);
	*/

	for(int i=0; i<SIGMA.size()-1; i++){
		if(abs(SIGMA[i+1].val() - SIGMA[i].val()) > 0.7*M_PI){
			//cout << i << endl;
			int aux = round(abs(SIGMA[i+1].val() - SIGMA[i].val()) / M_PI);
			for(int j=i+1; j<SIGMA.size(); j++){
				SIGMA[j] = SIGMA[j] - Var((double)aux*M_PI);
			}
		}
	}

	TF1* reta = new TF1("reta", "[0]*x+[1]",0,TSIGMA[-1].val());
	sigma2 = GetTGraph(TSIGMA, SIGMA);
	sigma2->SetTitle("4) Fit;t;#sigma");
	//sigma2->GetYaxis()->SetRangeUser(-0.2,M_PI+0.2);
	sigma2->SetMarkerStyle(8);
	sigma2->SetMarkerSize(0.5);

	canvas->cd(4);
	gPad->SetGrid();
	sigma2->Draw("APL");
	sigma2->Fit("reta","QR");

	osc_per = -2*M_PI / reta->GetParameter(0);


	// DATA FOR MAINTENANCE ======================================
	canvas->cd();
	char text_aux[516];
	sprintf(text_aux, "FILE = <%s>, r = %lf, phi = %lf", file_in.c_str(), r_max*2, atan2(Z[0].val(), Y[0].val()));
	TLatex* text = new TLatex(0.03, 0.5, text_aux);
	text->SetTextSize(0.015);
	text->Draw();

	// SAVE PLOT ======================================
	canvas->SaveAs(plot_out.c_str());

	// SAVE LOG ======================================
	FILE* fout = fopen(file_out.c_str(), "w");
	fprintf(fout, "%.14e", atan2(Z[0].val(), Y[0].val()));
	fprintf(fout, "\t%.14e", r_max*2);
	fprintf(fout, "\t%.14e", osc_per);
	fprintf(fout, "\t%.14e", raio_max);

	double px_mean = PX.subDataSet(n_points/4, n_points-1).getMean().val();
	fprintf(fout, "\t%.14e", px_mean);

	fclose(fout);

	// DELETE STUUUUUUFF ======================================
	for(int i = 0; i < n_cols; i++) 
		delete[] values[i];
	delete[] values;
	delete canvas;
	delete trajectory; 
	delete sigma;
	delete sigma2;
	delete elipse;
	for(int i = 0; i < 3; i++)
		delete[] dft_out[i];
	delete[] dft_out;
		
	//delete picardlindelof;
	delete reta;

	delete text;

}

int main(int argc, char** argv){

	if(argc == 4){
		mood(argv[1], argv[2], argv[3]);
		return 1;
	}

	return 0;

}