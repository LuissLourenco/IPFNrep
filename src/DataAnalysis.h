#ifndef _DataAnalysis_
#define _DataAnalysis_

#include <iostream>
#include <cstdio>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <ctime>
#include <dirent.h>

#include "TApplication.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TAxis.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPad.h"
#include "TGraph2DErrors.h"
#include "TF2.h"
#include "TApplication.h"
#include "TSystem.h"
#include "TMarker.h"
#include "TArrow.h"
#include "TGaxis.h"
#include "TColorGradient.h"
#include "TStyle.h"
#include "TText.h"
#include "TH2.h"
#include "TPolyLine3D.h"
#include "TEllipse.h"
#include "TBox.h"

using namespace std;

void printProgress(double percentage);

class Var{
	public:
		Var(double x=0, double y=0):value(x), error(y){
			//if(value == 0) value = 1E-60;
			//if(error == 0) error = 1E-60;			
		};
		~Var(){};

		string get(){return (to_string(value)+string(" \\pm ")+to_string(error)).c_str();}
		double val(){return value;}
		double err(){return error;}
		string print_latex(int p){cout << "$" << setprecision(p) << value; cout << "\\pm" << setprecision(p) << error << "$"; return "";}
		void setError(double e){error = e;}

		Var operator=(const Var& m1){value=m1.value; error=m1.error; return *this;}
		Var operator+(const Var& m1){return Var(value+m1.value, sqrt(error*error+m1.error*m1.error));}
		Var operator-(const Var& m1){return Var(value-m1.value, sqrt(error*error+m1.error*m1.error));}
		Var operator*(const Var& m1){return Var(value*m1.value, abs(value*m1.value) * sqrt(pow(error/value,2)+pow(m1.error/m1.value,2)));}
		Var operator/(const Var& m1){return Var(value/m1.value, abs(value/m1.value) * sqrt(pow(error/value,2)+pow(m1.error/m1.value,2)));}
		Var operator^(const Var& m1){return Var(pow(value, m1.value), abs(pow(value, m1.value)) * sqrt(pow(m1.value/value*error, 2) + pow(log(value)*m1.error, 2)));}

		bool operator<(const Var& m1){return value < m1.value;}
		bool operator>(const Var& m1){return value > m1.value;}
		bool operator<=(const Var& m1){return value <= m1.value;}
		bool operator>=(const Var& m1){return value >= m1.value;}
		bool operator==(const Var& m1){return value == m1.value;}
		bool operator!=(const Var& m1){return value != m1.value;}
		friend Var abs(const Var& v){return Var(abs(v.value), v.error);};
		friend Var sin(const Var& v){return Var(sin(v.value), abs(cos(v.value)*v.error));};
		friend Var cos(const Var& v){return Var(cos(v.value), abs(sin(v.value)*v.error));};
		friend Var sqrt(const Var& v){return Var(v.value, v.error)^Var(1./2, 0);};
		

		friend ostream& operator<<(ostream& os, const Var& v);

	private:
		double value;
		double error;

};

class DataSet{
	public:
		DataSet(int n=0, double v=0, double e=0):n(n){
			if(n == 0) {
				empty = true;
				n = 1;
			}
			else empty = false;
			var = new Var[n];
			for(int i = 0; i < n; i++) var[i] = Var(v, e);
		};
		DataSet(int n, double *values, double error=0):n(n){
			var = new Var[n];
			for(int i = 0; i < n; i++) var[i] = Var(values[i], error);
		};
		DataSet(int n, Var *values):n(n){
			var = new Var[n];
			for(int i = 0; i < n; i++) var[i] = values[i];
		};
		~DataSet(){/*delete[] var;*/};

		int size(){ return n;} 
		int size() const { return n;}
		Var getMax(){
			Var res = var[0];
			for(int i = 1; i < n; i++)
				if(var[i] > res)
					res = var[i];
			return res;
		}
		int getMaxI(){
			Var res = var[0];
			int resi = 0;
			for(int i = 1; i < n; i++)
				if(var[i] > res){
					res = var[i];
					resi = i;
				}
			return resi;
		}
		Var getMin(){
			Var res = var[0];
			for(int i = 1; i < n; i++)
				if(var[i] < res)
					res = var[i];
			return res;
		}
		int getMinI(){
			Var res = var[0];
			int resi = 0;
			for(int i = 1; i < n; i++)
				if(var[i] < res){
					res = var[i];
					resi = i;
				}
			return resi;
		}
		Var getMean(){
			double sum = 0;
			double max_error = 0;
			for(int i = 0; i < n; i++) 
				sum+=var[i].val();
			sum /= n;
			for(int i = 0; i < n; i++){
				if(var[i].err() > max_error) 
					max_error = var[i].err();
				if(abs(sum - var[i].val()) > max_error)
					max_error = abs(sum - var[i].val());
			}
			return Var(sum, max_error);
		}
		DataSet compress(int n_per_n){
			DataSet *aux;
			Var* v_aux = new Var[n_per_n];
			Var* res = new Var[n/n_per_n];
			for(int i = 0; i < n/n_per_n; i++){
				for(int j = 0; j < n_per_n; j++){
					v_aux[j] = var[i*n_per_n+j];
				}
				aux = new DataSet(n_per_n, v_aux);
				res[i] = aux->getMean();
			}
			return DataSet(n/n_per_n, res);
		}
		DataSet concat(DataSet s){
			DataSet res(n+s.n);
			for(int i = 0; i < n; i++) res[i] = var[i];
			for(int i = 0; i < s.n; i++) res[n+i] = s[i];
			return res;
		}
		double* array(){
			double* res = new double[n];
			for(int i = 0; i < n; i++) res[i] = var[i].val();
			return res;
		}
		DataSet setError(double e){
			for(int i = 0; i < n; i++) var[i].setError(e);
			return *this;
		}
		DataSet append(Var v){
			if(empty) n = 0;
			empty = false;
			Var *aux = new Var[n+1];
			for(int i = 0; i < n; i++){
				aux[i] = var[i];
			}
			aux[n] = Var(v.val(), v.err());
			var = aux;
			n+=1;
			return *this;
		}
		void print(){
			for(int i = 0; i < n; i++) cout << i+1 << ": " << var[i] << endl;
		}
		DataSet cut(int c){
			Var* aux = new Var[c];
			for(int i = 0; i < c; i++) aux[i] = var[i];
			n = c;
			delete[] var;
			var = new Var[n];
			for(int i = 0; i < n; i++) var[i] = aux[i];
			return *this;
		}
		DataSet subDataSet(int a, int b){
			Var* aux = new Var[b-a];
			for(int i = a; i < b; i++) aux[i-a] = var[i];
			return DataSet(b-a, aux);
		}



    	Var& operator[](int i){return var[(n+i)%n];}
		Var operator[](int i) const{return var[(n+i)%n];}
		DataSet& operator=(const DataSet &v){
		  	n = v.n;
		  	delete[] var;
		  	var = new Var[n];
		  	for(int i = 0; i < n; i++) var[i] = v.var[i];
		  	return *this;
		}	
		DataSet operator+(const DataSet &v){
			if(n != v.n){
				cout << "Tamanhos de vetores incompatíveis!" << endl;
				return DataSet();
			}
			Var *aux = new Var[n];
			for(int i = 0; i < n; i++) aux[i] = var[i] + v.var[i];
			DataSet res(n, aux);
			delete[] aux;
			return res;
		}
		DataSet operator+(Var v){ 
			Var *aux = new Var[n];
			for(int i = 0; i < n; i++) aux[i] = var[i] + v;
			DataSet res(n, aux);
			delete[] aux;
			return res;
		}
		DataSet operator-(const DataSet &v){ 
			if(n != v.n){
				cout << "Tamanhos de vetores incompatíveis!" << endl;
				return DataSet();
			}
			Var *aux = new Var[n];
			for(int i = 0; i < n; i++) aux[i] = var[i] - v.var[i];
			DataSet res(n, aux);
			delete[] aux;
			return res;
		}
		DataSet operator-(Var v){ 
			Var *aux = new Var[n];
			for(int i = 0; i < n; i++) aux[i] = var[i] - v;
			DataSet res(n, aux);
			delete[] aux;
			return res;
		}
		DataSet operator*(const DataSet &v){ 
			if(n != v.n){
				cout << "Tamanhos de vetores incompatíveis!" << endl;
				return DataSet();
			}
			Var *aux = new Var[n];
			for(int i = 0; i < n; i++) aux[i] = var[i] * v.var[i];
			DataSet res(n, aux);
			delete[] aux;
			return res;
		}
		DataSet operator*(Var v){ 
			Var *aux = new Var[n];
			for(int i = 0; i < n; i++) aux[i] = var[i] * v;
			DataSet res(n, aux);
			delete[] aux;
			return res;
		}
		DataSet operator/(const DataSet &v){ 
			if(n != v.n){
				cout << "Tamanhos de vetores incompatíveis!" << endl;
				return DataSet();
			}
			Var *aux = new Var[n];
			for(int i = 0; i < n; i++) aux[i] = var[i] / v.var[i];
			DataSet res(n, aux);
			delete[] aux;
			return res;
		}
		DataSet operator/(Var v){ 
			Var *aux = new Var[n];
			for(int i = 0; i < n; i++) aux[i] = var[i] / v;
			DataSet res(n, aux);
			delete[] aux;
			return res;
		}

		friend DataSet sqrt(const DataSet& v){
			Var *aux = new Var[v.size()];
			for(int i = 0; i < v.size(); i++) aux[i] = v[i]^Var(1./2, 0);
			DataSet res(v.size(), aux);
			delete[] aux;
			return res;
		};

		friend DataSet abs(const DataSet& v){
			Var *aux = new Var[v.size()];
			for(int i = 0; i < v.size(); i++) aux[i] = abs(v[i]);
			DataSet res(v.size(), aux);
			delete[] aux;
			return res;
		};

		friend DataSet atan2(const DataSet& v1, const DataSet& v2){
			if(v1.size() != v2.size()){
				cout << "Tamanhos incompatíveis de DataSet! (atan2)" << endl;
				return DataSet();
			}
			Var *aux = new Var[v1.size()];
			for(int i = 0; i < v1.size(); i++) aux[i] = atan2(v1[i].val(), v2[i].val());
			DataSet res(v1.size(), aux);
			delete[] aux;
			return res;
		};

	private:
		int n;
		Var *var;
		bool empty;

};

void FLAG(string s="AM HERE"){cout << "==========" << s << "==========" << endl;}
void FLAG(double s){cout << "==========" << s << "==========" << endl;}
void FLAG(int s){cout << "==========" << s << "==========" << endl;}
void print_latex_line(string name, Var *v, int n, int p);

double **ReadFile(string src_file, int *n_cols, int *n_points, bool header, bool print);

TGraphErrors *GetTGraphErrors(int n_points, Var* X, Var* Y);
TGraphErrors *GetTGraphErrors(DataSet X, DataSet Y);
TGraph *GetTGraph(DataSet X, DataSet Y);

TGraph2DErrors *GetTGraph2DErrors(DataSet X, DataSet Y, DataSet Z);
TGraph2D *GetTGraph2D(DataSet X, DataSet Y, DataSet Z);
Var GetVar(TF1 *f, int par);
double GetChi2NDF(TF1* f);

string* list_dir(const char *path, int *n_files);


double** computeDft(double dt, int N, double* in, int stop=0);
double** computeDft2(double dt, int N, double* in, double T_min, double T_max, double rate);
//returns period, frequency, value



#endif
