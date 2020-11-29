#include "../src/DataAnalysis.cpp"
#include <string>
#include <vector>

using namespace std;

void mood(string file_in, string plot_out, string file_out){

	int n_points, n_cols;
	double** values = ReadFile(file_in.c_str(), &n_cols, &n_points, true);

	// READ DATA
	DataSet T(n_points, values[0]);
	DataSet X(n_points, values[1]);
	DataSet Y(n_points, values[2]);
	DataSet Z(n_points, values[3]);
	DataSet PX(n_points, values[4]);
	DataSet PY(n_points, values[5]);
	DataSet PZ(n_points, values[6]);
	DataSet GAMMA(n_points, values[7]);

	double starting_radius = sqrt(Z[0]^2+Y[0]^2).val();

	// DISCARD FIRST POINTS
	double per_cut = 0.1;
	T = T.subDataSet(n_points*per_cut, n_points);
	X = X.subDataSet(n_points*per_cut, n_points);
	Y = Y.subDataSet(n_points*per_cut, n_points);
	Z = Z.subDataSet(n_points*per_cut, n_points);
	PX = PX.subDataSet(n_points*per_cut, n_points);
	PY = PY.subDataSet(n_points*per_cut, n_points);
	PZ = PZ.subDataSet(n_points*per_cut, n_points);
	GAMMA = GAMMA.subDataSet(n_points*per_cut, n_points);
	n_points = T.size();

	// CHECK IF IT DOES NOT EXPLODE
	if(sqrt(Y*Y+Z*Z).getMax().val() > 10){
		cout << "Radius diverges! Exiting..." << endl;
		return;
	}

	TCanvas* canvas = new TCanvas("", "", 2000, 1000);
	canvas->Divide(2, 2);

	// DRAW Y(T)
	canvas->cd(1);
	TGraph* g = GetTGraph(T, Y);
	g->Draw("AP");

	// DRAW DFT[Y(T)]
	canvas->cd(2);
	int n_samples;
	double** dft_data = computeDft2((T[1]-T[0]).val(), T.size(), Y.array(), 0.1, T[-1].val(), (T[1]-T[0]).val(), &n_samples);
	TGraph* dft = GetTGraph(DataSet(n_samples, dft_data[0]), DataSet(n_samples, dft_data[2]));
	dft->Draw("ALP");

	// CHECK PERIOD OS ELLIPSIS
	canvas->cd(1);
	int index = DataSet(n_samples, dft_data[2]).getMaxI();
	double osc_per = DataSet(n_samples, dft_data[0])[index].val();
	TF1* f = new TF1("", "0.2*cos(x*[0])", 0, 500);
	f->SetParameter(0, 2*M_PI/osc_per);
	f->SetNpx(1000);
	f->SetLineColor(kGreen);
	f->Draw("SAME");

	// GET ENVELOPE OF ABS(Y(T))
	canvas->cd(3);
	int wid = osc_per / (T[1]-T[0]).val() / 2;
	DataSet T2, Y2;
	int aux = 0;
	while(aux < n_points-wid){
		T2.append(T.subDataSet(aux, aux+wid).getMean());
		Y2.append(abs(Y).subDataSet(aux, aux+wid).getMax());
		aux += wid;
	}
	TGraph* g2 = GetTGraph(T2, Y2);
	g2->SetLineColor(kRed);
	TGraph* gc = GetTGraph(T, abs(Y));

	// GET ARCCOS(ENV(Y(T)))
	DataSet ARCCOS;
	double max_r = Y2.getMax().val();
	for(int i = 0; i < T2.size(); i++){
		ARCCOS.append(Var(acos(Y2[i].val()/max_r)));
	}

	TMultiGraph *multi = new TMultiGraph();
	multi->Add(gc);
	multi->Add(g2);
	multi->Draw("ALP");

	// COMPUTE DFT AND FINNER DFT TO GET PERIOD OF ROTATION
	canvas->cd(4);
	dft_data = computeDft2((T2[1]-T2[0]).val(), T2.size(), ARCCOS.array(), osc_per, T[-1].val(), (T[1]-T[0]).val(), &n_samples);
	index = DataSet(n_samples, dft_data[2]).getMaxI();
	osc_per = DataSet(n_samples, dft_data[0])[index].val();
	dft_data = computeDft2((T2[1]-T2[0]).val(), T2.size(), ARCCOS.array(), osc_per-(T[1]-T[0]).val(), osc_per+(T[1]-T[0]).val(), (T[1]-T[0]).val()/100, &n_samples);

	// DRAW DFT
	TGraph* dft2 = GetTGraph(DataSet(n_samples, dft_data[0]), DataSet(n_samples, dft_data[2]));
	dft2->Draw("ALP");

	// CHECK PERIOD OF ROTATION
	canvas->cd(3);
	index = DataSet(n_samples, dft_data[2]).getMaxI();
	osc_per = DataSet(n_samples, dft_data[0])[index].val();
	TF1* f2 = new TF1("", "0.6+0.3*cos(x*[0])", 0, 500);
	f2->SetParameter(0, 2*M_PI/osc_per);
	f2->SetNpx(1000);
	f2->SetLineColor(kGreen);
	f2->Draw("SAME");

	// SAVE PLOT
	canvas->SaveAs(plot_out.c_str());

	// SAVE LOG 
	double t_rot = 2*M_PI/f2->GetParameter(1)*2;
	double t_ell = 2*M_PI/f->GetParameter(1);

	FILE* fout = fopen(file_out.c_str(), "w");
	fprintf(fout, "%.14e", atan2(Z[0].val(), Y[0].val())); //phi0
	fprintf(fout, "\t%.14e", starting_radius); //r0  <- e mesmo, acredita bro
	fprintf(fout, "\t%.14e", t_rot);
	fprintf(fout, "\t%.14e", qrt(Z.getMax()^2+Y.getMax()^2).val());

	double px_mean = PX.getMean().val();
	fprintf(fout, "\t%.14e", px_mean); //px_mean
	fprintf(fout, "\t%.14e", X[-1].val()); //X_final
	fprintf(fout, "\t%.14e", osc_per/T_elipse); //eta

	fclose(fout);

}

int main(int argc, char** argv){

	if(argc == 4){
		mood(argv[1], argv[2], argv[3]);
		return 1;
	}

	return 0;

}