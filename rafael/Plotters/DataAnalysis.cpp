#include "DataAnalysis.h"

void printProgress(double percentage){
	int PBWIDTH = 60;
	char* PBSTR = new char[PBWIDTH];
	for(int i = 0; i < PBWIDTH; i++) PBSTR[i] = '|';
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush(stdout);
    if (val >= 100) cout << endl;
}

//////////////////////////////////////////////////////////////

ostream& operator<<(ostream& os, const Var& v){
    os << v.value << " \\pm " << v.error;
    return os;
}


//////////////////////////////////////////////////////////////

void print_latex_line(string name, Var *v, int n, int p){
		cout << "$" << name << "$";
		for(int i = 0; i < n; i++) 
				cout << " & " << v[i].print_latex(p);
		cout << " \\\\" << endl;

}


double **ReadFile(string src_file, int *n_cols, int *n_points, bool header=false){

		cout << "Reading values from: <" << src_file;
		if (header) cout << "> with header" << endl;
		else cout << "> without header" << endl;

		vector<string> line_vec;
		ifstream file(src_file);

		if(!file.is_open()){
				cout << "Ficheiro para leitura inexistente!" << endl;
				return 0;
		}

		string line;
		while(getline(file, line)) line_vec.push_back(line);
		file.close();

		if (header) line_vec.erase(line_vec.begin());

		for(int i = 0; i < line_vec.size(); i++){
				for(int j = 0; j < line_vec[i].size(); j++)
						if (line_vec[i][j] == 44) line_vec[i][j] = 46;
		}

		double aux;
		int n_col = 0;
		size_t aux_i = 0, aux_p = 0;
		while(aux_p < line_vec[0].size()){
				aux = stod(line_vec[0].substr(aux_p), &aux_i);
				n_col++;
				aux_p += aux_i;
		}

		double **res = new double*[n_col];
		for(int i = 0; i < n_col; i++) res[i] = new double[line_vec.size()];

		int j;
		for(int i = 0; i < line_vec.size(); i++){
				j = 0;
				aux_i = 0;
				aux_p = 0;
				while(aux_p < line_vec[i].size()){
						res[j][i] = stod(line_vec[i].substr(aux_p), &aux_i);
						aux_p += aux_i;
						j++;
				}
		}

		*n_cols = n_col;
		*n_points = line_vec.size();
		return res;

}


TGraphErrors *GetTGraphErrors(int n_points, Var* X, Var* Y){

	double *x = new double[n_points];
	double *y = new double[n_points];
	double *ex = new double[n_points];
	double *ey = new double[n_points];
	for(int i = 0; i < n_points; i++){
		x[i] = X[i].val();
		y[i] = Y[i].val();
		ex[i] = X[i].err();
		ey[i] = Y[i].err();
	}

	return new TGraphErrors(n_points, x, y, ex, ey);
	
}


TGraphErrors *GetTGraphErrors(DataSet X, DataSet Y){

	if(X.size() != Y.size()){
		cout << "Tamanho incompatível para TGraphErrors" << endl;
		return new TGraphErrors();
	}

	int n_points = X.size();

	double *x = new double[n_points];
	double *y = new double[n_points];
	double *ex = new double[n_points];
	double *ey = new double[n_points];
	for(int i = 0; i < n_points; i++){
		x[i] = X[i].val();
		y[i] = Y[i].val();
		ex[i] = X[i].err();
		ey[i] = Y[i].err();
	}

	return new TGraphErrors(n_points, x, y, ex, ey);
	
}

TGraph *GetTGraph(DataSet X, DataSet Y){

	if(X.size() != Y.size()){
		cout << "Tamanho incompatível para TGraphErrors" << endl;
		return new TGraph();
	}

	int n_points = X.size();

	double *x = new double[n_points];
	double *y = new double[n_points];
	for(int i = 0; i < n_points; i++){
		x[i] = X[i].val();
		y[i] = Y[i].val();
	}

	return new TGraph(n_points, x, y);
	
}

TGraph2DErrors *GetTGraph2DErrors(DataSet X, DataSet Y, DataSet Z){

	if(X.size() != Y.size() || X.size() != Z.size() || Y.size() != Z.size()){
		cout << "Tamanho incompatível para TGraph2DErrors" << endl;
		return new TGraph2DErrors();
	}

	int n_points = X.size();

	double *x = new double[n_points];
	double *y = new double[n_points];
	double *z = new double[n_points];
	double *ex = new double[n_points];
	double *ey = new double[n_points];
	double *ez = new double[n_points];
	for(int i = 0; i < n_points; i++){
		x[i] = X[i].val();
		y[i] = Y[i].val();
		z[i] = Z[i].val();
		ex[i] = X[i].err();
		ey[i] = Y[i].err();
		ez[i] = Z[i].err();
	}

	return new TGraph2DErrors(n_points, x, y, z, ex, ey, z);
	
}

TGraph2D *GetTGraph2D(DataSet X, DataSet Y, DataSet Z){

	if(X.size() != Y.size() || X.size() != Z.size() || Y.size() != Z.size()){
		cout << "Tamanho incompatível para TGraph2D" << endl;
		return new TGraph2D();
	}

	int n_points = X.size();

	double *x = new double[n_points];
	double *y = new double[n_points];
	double *z = new double[n_points];
	for(int i = 0; i < n_points; i++){
		x[i] = X[i].val();
		y[i] = Y[i].val();
		z[i] = Z[i].val();
	}

	return new TGraph2D(n_points, x, y, z);
	
}

Var GetVar(TF1 *f, int par){return Var(f->GetParameter(par), f->GetParError(par));};
double GetChi2NDF(TF1* f){return f->GetChisquare() / f->GetNDF();}
