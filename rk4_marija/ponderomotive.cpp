#include "trial.c"
#include "Read.cpp"

void ponderomotive(){
	
	double r_min = -20;
	double r_max = 20;
	double dr = 0.01;
	double r = r_min;
	FILE *file;

	vector<vector<double>> sol;
	int n;
	int n_data = (int)((r_max - r)/dr);

	double* r_arr = new double[n_data];
	double* y_f = new double[n_data];
	double* p_y_max = new double[n_data];
	double* p_y_f = new double[n_data];
	double* p_x_f = new double[n_data];

	int counter = 0;
	double E0_list[6] = {10, -10, 10, -10, 10, -10};
	double w0_list[6] = {5, 5, 10, 10, 15, 15};

	for(counter = 2; counter < 6; counter++){
		int i = 0;
		r = r_min;
		while(r <= r_max){
			file = fopen("InputToBatch.txt","w");
			fprintf(file, "trash|x01|x02|x03|p01|p02|p03|T|N|pri|E0|w0\n"); 
			fprintf(file, "0.\n%lf\n0.\n0.\n0.\n0.\n150\n100000\n100\n%lf\n%lf", r, E0_list[counter], w0_list[counter]);
			fclose(file);
		
			main();

			sol = ReadFile2Vec("Out0.txt");
			n = sol.size();

			y_f[i] = sol[n-1][2];
			p_x_f[i] = sol[n-1][4];
			p_y_f[i] = sol[n-1][5];
			p_y_max[i] = 0;
			for(int j = 0; j < n; j++){
				if(abs(sol[j][5]) > p_y_max[i]) p_y_max[i] = abs(sol[j][5]);
			}

			cout << "Counter " << counter+1 << " of 6\t";
			cout << "Step " << i << " of " << n_data << "; r = " << r << endl;

			r_arr[i] = r;
			r += dr;
			i++;
		}

		file = fopen(("Data_Ponderomotive_v"+to_string(counter+1)+".txt").c_str(),"w");
		fprintf(file, "r\ty_f\tp_y_max\tp_x_f\tp_y_f\tw0=%lf|lambda=1|E0=%lf|time=150"
						, w0_list[counter], E0_list[counter]); 
		for(int i = 0; i < n_data; i++){
			fprintf(file, "\n%lf\t%lf\t%lf\t%lf\t%lf", 
					r_arr[i], y_f[i], p_y_max[i], p_x_f[i], p_y_f[i]);
		}
		fclose(file);

	}

/*

	double r = -20;
	double r_max = 20;
	double dr = 1;
	FILE *file;

	vector<vector<double>> sol;
	int n;
	int n_data = (int)((r_max - r)/dr);

	while(r <= r_max){
		file = fopen("InputToBatch.txt","w");
		fprintf(file, "trash|x01|x02|x03|p01|p02|p03|T|N|pri\n"); 
		fprintf(file, "0.\n%lf\n0.\n0.\n0.\n0.\n150\n20000\n100", r);
		fclose(file);
	
		main();

		sol = ReadFile2Vec("Out0.txt");
		n = sol.size();

		file = fopen(("Movement_r"+to_string(r)+".txt").c_str(),"w");
		fprintf(file, "t\tx\ty\tpx\tpy"); 
		for(int i = 0; i < n; i++){
			fprintf(file, "\n%lf\t%lf\t%lf\t%lf\t%lf", 
					sol[i][0], sol[i][1], sol[i][2], sol[i][4], sol[i][5]);
		}
		fclose(file);

		r += dr;
	}

	*/
	

	
}