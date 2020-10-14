#include "../src/DataAnalysis.cpp"
#include <string>
#include <vector>

using namespace std;

void mood(string file_in, string plot_out, string file_out){

	int n_points, n_cols;
	double** values = ReadFile(file_in.c_str(), &n_cols, &n_points, true);

	cout << "Saving to " << plot_out << "\n          " << file_out << endl;

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
	canvas->Divide(3, 2);
	TGraph* trajectory; double r_max = sqrt(Y[0]*Y[0]+Z[0]*Z[0]).val()/2;
	TGraph* sigma;
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
	elipse->GetXaxis()->SetLimits(-M_PI, M_PI);
	elipse->GetYaxis()->SetRangeUser(1E-3, r_max);
	elipse->SetTitle("1) Fitando a Elipse;#theta;r");
	elipse->Draw("AP");
	fit_elipse->Draw("SAME");

	// CHECK THE ELLIPSE FITTING ======================================
	canvas->cd(2);
	trajectory->GetXaxis()->SetLimits(-r_max, r_max);
	trajectory->GetYaxis()->SetRangeUser(-r_max, r_max);
	trajectory->SetTitle("2) Checkando a Elipse;y;z");
	trajectory->Draw("ACP");
	elipse_draw = new TEllipse(0, 0, fit_elipse->GetParameter(0), fit_elipse->GetParameter(1), 0, 360, -fit_elipse->GetParameter(2)*180/M_PI);
	elipse_draw->SetFillStyle(0);
	elipse_draw->SetLineColor(3);
	elipse_draw->Draw("SAME");

	// FOURIER THE PERIOD ======================================
	canvas->cd(3);
	int dft_to_compute = SIGMA.size()/5;
	dft_out = computeDft((TSIGMA[1]-TSIGMA[0]).val(), SIGMA.size(), SIGMA.array(), dft_to_compute);
	index = DataSet(dft_to_compute, dft_out[2]).getMaxI();

	dft_out = computeDft2((TSIGMA[1]-TSIGMA[0]).val(), SIGMA.size(), SIGMA.array(), DataSet(dft_to_compute, dft_out[0])[index+1].val(), DataSet(dft_to_compute, dft_out[0])[index-1].val(), 0.1, &dft_to_compute);
	index = DataSet(dft_to_compute, dft_out[2]).getMaxI(); 
	osc_per = DataSet(dft_to_compute, dft_out[0])[index].val();

	TGraph* dft_g = GetTGraph(DataSet(dft_to_compute, dft_out[0]), DataSet(dft_to_compute, dft_out[2]));
	dft_g->SetTitle("3) Fouriando o periodo;T;AMP");
	dft_g->Draw("ALP");

	// SIGMA WITH PERIOD NOT NORMALIZED ======================================
	canvas->cd(4);
	sigma = GetTGraph(TSIGMA, SIGMA);
	sigma->SetTitle("4) Desenhando cenas nao normalizadas;t;#sigma");
	sigma->Draw("ALP");	canvas->cd(2);
	TF1* per_fit = new TF1("", "1.5+cos([0]*x)", 0, TSIGMA[-1].val());
	per_fit->SetParameter(0, 2*M_PI/osc_per);
	per_fit->Draw("SAME");

	// GET THE AMPLITUDE AND MEAN OF OSCILATION ======================================
	canvas->cd(5);
	bool centered_on_zero;
	TH1D* hist = new TH1D("", "", 400, -M_PI, M_PI);
	for(int i = 0; i < SIGMA.size(); i++){
		hist->Fill(SIGMA[i].val());
		hist->Fill(SIGMA[i].val()-M_PI);
	}
	if(hist->Integral(290, 310) < 2){
		centered_on_zero = true;
		hist = new TH1D("", "5) Distribuicao dos sigmas para ver a amplitude", 400, -M_PI/2, M_PI/2);
		for(int i = 0; i < SIGMA.size(); i++){
			hist->Fill(SIGMA[i].val());
			hist->Fill(SIGMA[i].val()-M_PI);
		}
	}else{
		centered_on_zero = false;
		hist = new TH1D("", "5) Distribuicao dos sigmas para ver a amplitude", 400, 0, M_PI);
		for(int i = 0; i < SIGMA.size(); i++){
			hist->Fill(SIGMA[i].val());
			hist->Fill(SIGMA[i].val()-M_PI);
		}
	}
	hist->Draw();
	int nq = 100;
	double* xq = new double[nq];
	double* piq = new double[nq];
	double* yq = new double[nq];
	for(int i = 0; i < nq; i++){
		xq[i] = (double)(i+1)/nq;
		piq[i] = (double)(i+1)/nq*M_PI;
	}
	hist->GetQuantiles(nq,yq,xq);
	for(int i = 0; i < nq; i++)
		piq[i] = hist->GetMaximum()*(i+1)/nq;
	TGraph* auxg = new TGraph(nq, yq, piq);
	auxg->Draw("SAME PC");
	double amp_min = yq[5];
	double amp_max = yq[95];

	// FINAL FUCKING PRODUCT ======================================
	canvas->cd(6);
	TGraph* final_graph = new TGraph(TSIGMA.size());
	final_graph->SetTitle("6) Finalmente...;t;#sigma");
	if(centered_on_zero){
		for(int i = 0; i < TSIGMA.size(); i++){
			if(SIGMA[i].val() > M_PI/2)
				final_graph->SetPoint(i, TSIGMA[i].val(), SIGMA[i].val()-M_PI);
			else
				final_graph->SetPoint(i, TSIGMA[i].val(), SIGMA[i].val());
		}
	}else{
		for(int i = 0; i < TSIGMA.size(); i++){
			final_graph->SetPoint(i, TSIGMA[i].val(), SIGMA[i].val());
		}
	}
	final_graph->Draw("ACP");
	TF1* final_fit = new TF1("", "[0] + [1]*sin([2]*x+[3])", 0, TSIGMA[-1].val());
	final_fit->FixParameter(0, (amp_max+amp_min)/2);
	final_fit->FixParameter(1, (amp_max-amp_min)/2);
	final_fit->FixParameter(2, 2*M_PI/osc_per);
	final_graph->Fit(final_fit, "Q");
	final_fit->Draw("SAME");

	// DATA FOR MAINTENANCE ======================================
	canvas->cd();
	char text_aux[516];
	sprintf(text_aux, "FILE = <%s>, r = %lf, phi = %lf", file_in.c_str(), r_max*2, atan2(Z[0].val(), Y[0].val()));
	TLatex* text = new TLatex(0.1, 0.49, text_aux);
	text->SetTextSize(0.02);
	text->Draw();

	// SAVE PLOT ======================================
	canvas->SaveAs(plot_out.c_str());

	// SAVE LOG ======================================
	FILE* fout = fopen(file_out.c_str(), "w");
	fprintf(fout, "%.14e", Y[0].val());
	fprintf(fout, "\t%.14e", Z[0].val());
	fprintf(fout, "\t%.14e", atan2(Z[0].val(), Y[0].val()));
	fprintf(fout, "\t%.14e", r_max*2);
	fprintf(fout, "\t%.14e", osc_per);
	fprintf(fout, "\t%.14e", amp_min);
	fprintf(fout, "\t%.14e", amp_max);
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
	delete elipse;
	for(int i = 0; i < 3; i++)
		delete[] dft_out[i];
	delete[] dft_out;
	delete dft_g;
	delete per_fit;
	delete hist;
	delete[] xq;
	delete[] piq;
	delete[] yq;
	delete auxg;
	delete final_graph;
	delete final_fit;
	delete text;

}

int main(int argc, char** argv){

	if(argc == 4){
		mood(argv[1], argv[2], argv[3]);
		return 1;
	}

	return 0;

}