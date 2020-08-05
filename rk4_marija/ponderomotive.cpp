#include "trial.c"
#include "Read.cpp"

void ponderomotive(){
	
	double r = -20;
	double r_max = 20;
	double dr = 0.1;
	FILE *file;

	vector<vector<double>> sol;
	int n;
	int n_data = (int)((r_max - r)/dr);

	double** data_x  = new double*[n_data];
	double** data_y  = new double*[n_data];
	double** data_px = new double*[n_data];
	double** data_py = new double*[n_data];

	double* r_arr = new double[n_data];
	double* y_f = new double[n_data];
	double* p_y_max = new double[n_data];

	int i = 0;
	while(r <= r_max){
		file = fopen("InputToBatch.txt","w");
		fprintf(file, "trash|x01|x02|x03|p01|p02|p03|T|N|pri\n"); 
		fprintf(file, "0.\n%lf\n0.\n0.\n0.\n0.\n150\n20000\n100", r);
		fclose(file);
	
		main();

		sol = ReadFile2Vec("Out0.txt");
		n = sol.size();

		data_x [i] = new double[n];
		data_y [i] = new double[n];
		data_px[i] = new double[n];
		data_py[i] = new double[n];	
		y_f[i] = sol[n-1][2];
		p_y_max[i] = 0;
		for(int j = 0; j < n; j++){
			data_x [i][j] = sol[j][1];
			data_y [i][j] = sol[j][2];
			data_px[i][j] = sol[j][4];
			data_py[i][j] = sol[j][5];
			if(abs(data_py[i][j]) > p_y_max[i]) p_y_max[i] = abs(data_py[i][j]);
		}

		cout << "Step " << i << " of " << n_data << "; r = " << r << endl;

		r_arr[i] = r;
		r += dr;
		i++;
	}

	file = fopen("Data_Ponderomotive.txt","w");
	fprintf(file, "r\ty_f\tp_y_max"); 
	for(int i = 0; i < n_data; i++){
		fprintf(file, "\n%lf\t%lf\t%lf", 
				r_arr[i], y_f[i], p_y_max[i]);
	}
	fclose(file);

	TCanvas* c1 = new TCanvas("c1", "", 2400, 1200);
	TGraph** graph = new TGraph*[n_data];
	TMultiGraph* multi = new TMultiGraph();
	for(int i = 0; i < n_data; i++){
		graph[i] = new TGraph(n, data_x[i], data_y[i]);
		graph[i]->SetLineColor((1+i)%10);
		graph[i]->SetMarkerColor((1+i)%10);
		multi->Add(graph[i]);
	}	
	multi->Draw("ACP");
	c1->SaveAs("Ponderomotive.png");
	

	
}