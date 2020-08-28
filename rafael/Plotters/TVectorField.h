#ifndef _TVectorField_
#define _TVectorField_

#include <iostream>
#include <cstdio>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <ctime>

#include "DataAnalysis.h"

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
#include "TH1D.h"
#include "TColor.h"

using namespace std;

class TVectorField{
	public:
		TVectorField(int nr, double* x1r, double* y1r, double* x2r, double* y2r);
		TVectorField(DataSet x1r, DataSet y1r, DataSet x2r, DataSet y2r);
		~TVectorField();

		void Draw(string options);
		TH1D* GetTH1(){return f;}

		void SetLimits(double v1, double v2, double v3, double v4);
		
	private:
		int n;
		double* x1;
		double* y1;
		double* x2;
		double* y2;
		TH1D* f;
		TGraph2D* graph;
		TArrow** arr;
		double x_min, x_max, y_min, y_max;
		double scale = 0.5;

		void Draw_A();
		void Draw_F();

};


#endif
