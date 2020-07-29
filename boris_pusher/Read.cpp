#ifndef _READTOOL_H_
#define _READTOOL_H_

#include <iostream>
#include <string>
#include <cstddef>
#include <vector>
#include <fstream>
#include <sstream>

using namespace std;


vector<string> ReadFile2String(string src){

  vector<string> res;
  ifstream file(src);

  if(!file.is_open()){
    cout<<"File not found"<<endl;
    return res;
  }

  string line;

  while(getline(file, line)){
    res.push_back(line);
    if(line.empty())
      break;
  }

  file.close();

  return res;
}

vector<vector<double>> ReadFile2Vec(string src){

  vector<string> vstr = ReadFile2String(src);
  vector<vector<double>> res;

  string aux1;
  vector<double> aux;

  for(int i=0; i<(int)vstr.size(); i++){

    stringstream ss(vstr[i]);

    aux.clear();

    while(ss>>aux1)
      aux.push_back(stod(aux1));

    res.push_back(aux);
  }

  return res;
}	

double** ReadColumns2Array(string src){

  vector<vector<double>> vec = ReadFile2Vec(src);

  int ncol = vec[0].size();
  int nelm = vec.size();

  double** res = new double*[ncol];
  for(int i=0; i<ncol; i++)
    res[i] = new double[nelm];

  for(int i=0; i<ncol; i++)
    for(int j=0; j<nelm; j++)
      res[i][j] = vec[j][i];

  return res;
}


#endif

