#ifndef _ODEsolver_H
#define _ODEsolver_H

#include <iostream>
#include <vector>
#include <TFormula.h>
#include "ODEpoint.h"

using namespace std;

class ODEsolver{
  public:
    ODEsolver(vector<TFormula> Form);
    ~ODEsolver();
    vector<ODEpoint> Eulersolver(const ODEpoint& P0, double xmin, double xmax, double h_step); //DONE
    vector<ODEpoint> Heun(const ODEpoint& P0, double xmin, double xmax, double h_step); //DONE
    vector<ODEpoint> RK2solver(const ODEpoint& P0, double xmin, double xmax, double h_step); //DONE
    //vector<ODEpoint> RK4_AdapStep(const ODEpoint& P0, double xmin, double xmax, double h_step);
    vector<ODEpoint> RK4solver(const ODEpoint& P0, double xmin, double xmax, double h_step);

    void SetODEfunc(vector<TFormula> Form);

  private:
    ODEpoint EULER_iterator (const ODEpoint&, double step); //DONE
    ODEpoint Heun_iterator (const ODEpoint&, double step); //DONE
    ODEpoint RK2_iterator (const ODEpoint&, double step); //DONE
    //ODEpoint RK4_AS_iterator(const ODEpoint&, double step, vector<vector<double> >& K);
    ODEpoint RK4_iterator(const ODEpoint&, double step);

    vector<TFormula> F;
};

#endif
