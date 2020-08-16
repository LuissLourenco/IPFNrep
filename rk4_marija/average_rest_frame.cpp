#include "trial.c"
#include "DataAnalysis.cpp"

void average_rest_frame(){
	
	double a0_min = 0.1;
	double a0_max = 10;
	double da0 = 0.1;
	double a0 = a0_min;
	FILE *file;

	double px, pxs;
	DataSet X;
	double x_min, x_max;
	int n_cols, n_points;
	double** values;

	double dx_max = 0.0001;

	FILE *filewrite = fopen("Data_Average.txt","w");
	fprintf(filewrite, "a0\tpx"); 

	while(a0 <= a0_max){

		//px = 0;
		pxs = 0.01;

		while(true){

			file = fopen("InputToBatch.txt","w");
			fprintf(file, "trash|x01|x02|x03|p01|p02|p03|T|N|pri|E0|w0\n"); 
			fprintf(file, "0.\n0.\n0.\n%lf\n0.\n0.\n150\n100000\n100\n%lf\n0.", px, a0);
			fclose(file);

			main();

			values = ReadFile("Out0.txt", &n_cols, &n_points, false);

			X = DataSet(n_points, values[1]);
			x_max = abs(X.getMax().val());
			x_min = abs(X.getMin().val());

			if(abs(x_max - x_min) > dx_max){
				px -= (x_max - x_min) * pxs;
			}else{
				break;
			}

			cout << x_min << " - " << x_max << " => " << px << endl;

		}

		fprintf(filewrite, "\n%lf\t%lf", a0, px);
		cout << a0 << " => " << px << endl;
		a0 += da0;

	}
	
	fclose(filewrite);
	
}