#include "ODEpoint.h"

//CONSTRUCTORS
ODEpoint::ODEpoint(){
  t = 0;
  var.push_back(0);
  Ndim = 1;
}

ODEpoint::ODEpoint(double tval, double* funct, int Ndimf){
  t = tval;
  for(int i = 0; i < Ndimf; i++) var.push_back(funct[i]);
  Ndim = Ndimf;
}

ODEpoint::ODEpoint(double tval, vector<double> funct){
  t = tval;
  Ndim = funct.size();
  var = funct;
}

ODEpoint::ODEpoint(const ODEpoint& P){
  t = P.t;
  Ndim = P.Ndim;
  var = P.var;
  //for(int i = 0; i < Ndim; i++) var.push_back(P.var[i]);
}

//DESTRUCTOR
ODEpoint::~ODEpoint(){

}


//GETS AND SETS
vector<double> ODEpoint::Get_Var_vec() const{
  return var;
}

double* ODEpoint::Get_Var_ptr() const{
  double *res = new double[Ndim];
  for(int i = 0; i < Ndim; i++) res[i] = var[i];
  return res;
}

double* ODEpoint::Get_VarTime() const{
  double *res = new double[Ndim+1];
  for(int i = 0; i < Ndim; i++) res[i] = var[i];
  res[Ndim] = t;
  return res;
}

int ODEpoint::GetNdim() const{
  return Ndim;
}

double ODEpoint::Get_Time() const{
  return t;
}	

void ODEpoint::Set_Time(double tval){
  t = tval;
}

void ODEpoint::Set_Var(vector<double> funct){
  var = funct; 
  Ndim = var.size();
}



//OPERATORS 
ODEpoint& ODEpoint::operator=(const ODEpoint& P){
  Ndim = P.Ndim;
  t = P.t;
  var.clear();
  var = P.var;
  return *this;
}

const double& ODEpoint::operator[](int i) const{
  return var[i];
}

double& ODEpoint::operator[] (int i){
  return var[i];
}

void ODEpoint::Print() const{
  cout << "t = " << t << "\t y = (";
  for(int i = 0; i < Ndim; i++){
    if(i != 0) cout << ", ";
    cout << var[i]; 
  } 
  cout << ")" << endl;
}
