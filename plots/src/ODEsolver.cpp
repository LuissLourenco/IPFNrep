#include "ODEsolver.h"

//CONSTRUCTOR
ODEsolver::ODEsolver(vector<TFormula> Form){
  F = Form;
}


//DESTRUCTOR
ODEsolver::~ODEsolver(){

}


//METHODS
vector<ODEpoint> ODEsolver::Eulersolver(const ODEpoint& P0, double xmin, double xmax, double h_step){ //DONE

  vector<ODEpoint> points;
  points.push_back(P0);

  double x = xmin;
  while(x < xmax){
    points.push_back(EULER_iterator(points[points.size()-1], h_step));
    x += h_step;
  }

  return points;
}

vector<ODEpoint> ODEsolver::Heun(const ODEpoint& P0, double xmin, double xmax, double h_step){ //DONE

  vector<ODEpoint> points;
  points.push_back(P0);

  double x = xmin;
  while(x < xmax){
    points.push_back(Heun_iterator(points[points.size()-1], h_step));
    x += h_step;
  }

  return points;
}


vector<ODEpoint> ODEsolver::RK2solver(const ODEpoint& P0, double xmin, double xmax, double h_step){ //DONE

  vector<ODEpoint> points;
  points.push_back(P0);

  double x = xmin;
  while(x < xmax){
    points.push_back(RK2_iterator(points[points.size()-1], h_step));
    x += h_step;
  }

  return points;

}

vector<ODEpoint> ODEsolver::RK4solver(const ODEpoint& P0, double xmin, double xmax, double h_step){ //DONE

  vector<ODEpoint> points;
  double x = xmin;
  //points.push_back(P0);
  points.push_back(RK4_iterator(P0, h_step));
    x += h_step;
  while(x<xmax){
    points.push_back(RK4_iterator(points[points.size()-1], h_step));
    x += h_step;
  }

  return points;
}

/*vector<ODEpoint> ODEsolver::RK4_AdapStep(const ODEpoint& P0, double xmin, double xmax, double h_step){

  }*/


void ODEsolver::SetODEfunc(vector<TFormula> Form){
  F = Form;
}



//PRIVATE METHODS
ODEpoint ODEsolver::EULER_iterator (const ODEpoint& pi, double step){ //DONE
  double *y_next = new double[pi.GetNdim()];
  double *var = pi.Get_Var_ptr();
  for(int i = 0; i < pi.GetNdim(); i++){
    y_next[i] = var[i] + step * F[i].EvalPar(pi.Get_VarTime());
  }
  return ODEpoint(pi.Get_Time()+step, y_next, pi.GetNdim());
}

ODEpoint ODEsolver::Heun_iterator (const ODEpoint& pi, double step){ //DONE

  vector<double> y_next;
  vector<double> var = pi.Get_Var_vec();

  for(int i = 0; i < pi.GetNdim(); i++){
    y_next.push_back(var[i] + step/2 * (F[i].EvalPar(pi.Get_VarTime()) + F[i].EvalPar(EULER_iterator(pi, step).Get_VarTime())));
  }
  return ODEpoint(pi.Get_Time()+step, y_next);

}

ODEpoint ODEsolver::RK2_iterator (const ODEpoint& pi, double step){ //DONE

  vector<double> y_next;
  vector<double> var = pi.Get_Var_vec();

  for(int i = 0; i < pi.GetNdim(); i++){
    y_next.push_back(var[i] + step * F[i].EvalPar(EULER_iterator(pi, step/2).Get_VarTime()));
  }
  return ODEpoint(pi.Get_Time()+step, y_next);

}

/*ODEpoint ODEsolver::RK4_AS_iterator(const ODEpoint& pi, double step, vector<vector<double> >& K){

  }*/

ODEpoint ODEsolver::RK4_iterator(const ODEpoint& pi, double step){

  int N = pi.GetNdim();
  double* y_next = new double[N]; 
  double* k1 = new double[N];
  double* k2 = new double[N];
  double* k3 = new double[N];
  double* k4 = new double[N];
  double* var1 = pi.Get_VarTime();
  double* var2 = pi.Get_VarTime();
  double* var3 = pi.Get_VarTime();
  double* var4 = pi.Get_VarTime();

  for(int i=0; i<N; i++){
    k1[i] = step * F[i].EvalPar(var1);
    var2[i] += 0.5*k1[i];
  }
  var2[N] += step/2;

  for(int i=0; i<N; i++){
    k2[i] = step * F[i].EvalPar(var2);
    var3[i] += 0.5*k2[i];
  }
  var3[N] += step/2;
  
  for(int i=0; i<N; i++){
    k3[i] = step * F[i].EvalPar(var3);
    var4[i] += k3[i];
  }
  var4[N] += step;

  for(int i=0; i<N; i++){
    k4[i] = step * F[i].EvalPar(var4);
    y_next[i] = var1[i] + (k1[i]+2*k2[i]+2*k3[i]+k4[i])/6;
  }

  ODEpoint res(pi.Get_Time()+step, y_next, N);

  delete[] y_next;
  delete[] var1;
  delete[] var2;
  delete[] var3;
  delete[] var4;
  delete[] k1;
  delete[] k2;
  delete[] k3;
  delete[] k4;


  return res;
}

TFormula& ODEsolver::GetFormula(int i){ 
  if(i < 0) i = 0;
  else if(i >= F.size()) i = F.size()-1;
  return F[i];
}