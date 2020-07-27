#ifndef _ODEpoint_H
#define _ODEpoint_H

#include <iostream>
#include <vector>

using namespace std;

class ODEpoint{   

public: 
    ODEpoint(); //default constructor 
    ~ODEpoint(); //destructor    
    ODEpoint(double tval, double* funct, int Ndimf); //using double  *   
    ODEpoint(double tval, vector<double> funct); //using vector<double>  
    ODEpoint(const ODEpoint&); //copy  constructor 
    
    //member    access  funcEons    
    vector<double> Get_Var_vec() const;  //return    the y1,…,yNdim  dependent   variables   
    double *Get_Var_ptr() const;  //same but  as  double  *
    double *Get_VarTime() const;  //first the y1,…,yNdim  then    t    
    double* Get(){return Get_VarTime();};   
    int GetNdim() const;  //return    the number  of  dependent   variables   
    double Get_Time() const;  //return    only    the independent variable    (t) 
    void Set_Time(double tval);  //Set   dependent   variable    
    void Set_Var(vector<double> funct);   

    //operators 
    ODEpoint &operator=(const ODEpoint& P); 
    const double &operator[](int i) const;  
    double& operator[] (int i); 
    void Print() const;  

private:    
    double t;  
    vector<double> var;    
    int Ndim; 

};  

#endif
