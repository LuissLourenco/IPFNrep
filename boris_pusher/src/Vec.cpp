#include <iostream>
#include "Vec.h"

#include <iomanip>      // std::setprecision


using namespace std;

//----------------------------------------------------------------------
//Constuctors & Destructor

Vec::Vec(int n, double a){ //DONE
  N = n;
  entries = new double[n];
  for(int i = 0; i < N; i++) entries[i] = a;
}

Vec::Vec(int n, double *a){ //DONE
  N = n;
  entries = new double[n];
  for(int i = 0; i < N; i++) entries[i] = a[i];
}

Vec::Vec(const Vec &v0){ //DONE
  N = v0.N;
  entries = new double[N];
  for(int i = 0; i < N; i++) entries[i] = v0.entries[i];
}

Vec::~Vec(){ //DONE
  delete[] entries;
}

//Functions

void Vec::SetEntries (int n, double* a){ //DONE
  N = n;
  delete[] entries;
  entries = new double[n];
  for(int i = 0; i < N; i++) entries[i] = a[i];
}

void Vec::Print(){ //DONE
  cout << "(";
  for(int i = 0; i < N; i++){
    cout << setprecision(20) << entries[i];
    if(i != N-1) cout << ", ";
  }
  cout << ")" << endl;
}

int Vec::size(){ return N;} //DONE

int Vec::size() const { return N;} //DONE

double Vec::dot(const Vec &v){ //DONE
  if(N != v.N){
    cout << "Tamanhos de vetores incompatíveis!" << endl;
    return 0;
  }
  double res = 0;
  for(int i = 0; i < N; i++) res += entries[i] * v.entries[i];
  return res;
}

void Vec::swap(int i1, int i2){ //DONE
  if(i1 < 0 || i1 >= N || i2 < 0 || i2 >= N){
    cout << "Índices fora do intervalo." << endl;
    return;
  }
  double aux = entries[i1];
  entries[i1] = entries[i2];
  entries[i2] = aux;
}


//Overloaded Functions

Vec& Vec::operator=(const Vec &v){ //DONE
  N = v.N;
  delete[] entries;
  entries = new double[N];
  for(int i = 0; i < N; i++) entries[i] = v.entries[i];
  return *this;
}

Vec& Vec::operator+=(const Vec &v){ //DONE
  if(N != v.N){
    cout << "Tamanhos de vetores incompatíveis!" << endl;
    return *this;
  }
  for(int i = 0; i < N; i++) entries[i] += v.entries[i];
  return *this;
}

Vec Vec::operator+(const Vec &v){ //DONE
  if(N != v.N){
    cout << "Tamanhos de vetores incompatíveis!" << endl;
    return Vec();
  }
  double *aux = new double[N];
  for(int i = 0; i < N; i++) aux[i] = entries[i] + v.entries[i];
  Vec res(N, aux);
  delete[] aux;
  return res;
}

Vec& Vec::operator-=(const Vec &v){ //DONE
  if(N != v.N){
    cout << "Tamanhos de vetores incompatíveis!" << endl;
    return *this;
  }
  for(int i = 0; i < N; i++) entries[i] -= v.entries[i];
  return *this;
}

Vec Vec::operator-(const Vec &v){ //DONE
  if(N != v.N){
    cout << "Tamanhos de vetores incompatíveis!" << endl;
    return Vec(); 
  }
  double *aux = new double[N];
  for(int i = 0; i < N; i++) aux[i] = entries[i] - v.entries[i];
  Vec res(N, aux);
  delete[] aux;
  return res;
}

double& Vec::operator[](int i){ //DONE
  if(i < 0) i = 0;
  else if(i >= N) i = N-1;
  return entries[i];
}

double Vec::operator[](int i) const{
  if(i < 0) i = 0;
  else if(i >= N) i = N-1;
  return entries[i];
}

Vec Vec::operator-(){ //DONEInvalid read of size 8

  double *aux = new double[N];
  for(int i = 0; i < N; i++) aux[i] = - entries[i];
  Vec res(N, aux);
  delete[] aux;
  return res;
}

Vec Vec::operator+(){ //EPAH IDK
  return *this;
}

Vec Vec::operator*(const Vec &v){ //DONE
  if(N != v.N){
    cout << "Tamanhos de vetores incompatíveis!" << endl;
    return Vec();
  }
  double *aux = new double[N];
  for(int i = 0; i < N; i++) aux[i] = entries[i] * v.entries[i];
  Vec res(N, aux);
  delete[] aux;
  return res;
}

Vec Vec::operator*(double c){ //DONE
  double *aux = new double[N];
  for(int i = 0; i < N; i++) aux[i] = entries[i] * c;
  Vec res(N, aux);
  delete[] aux;
  return res;
}

Vec& Vec::operator*=(const Vec &v){ //DONE
  if(N != v.N){
    cout << "Tamanhos de vetores incompatíveis!" << endl;
    return *this;
  }
  for(int i = 0; i < N; i++) entries[i] *= v.entries[i];
  return *this;
}

Vec& Vec::operator*=(double c){ //DONE
  for(int i = 0; i < N; i++) entries[i] *= c;
  return *this;
}



//Getters & Setters
double Vec::At(int i){ //DONE
  if(i < 0) i = 0;
  else if(i >= N) i = N-1;
  return entries[i];
}

void Vec::SetAt(int i, double a){ //DONE
  if(i < 0) i = 0;
  else if(i >= N) i = N-1;
  entries[i] = a;
}

