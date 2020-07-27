#ifndef _Vec_H_
#define _Vec_H_

#include <utility>
#include <vector>
#include <cmath>
#include <cstdarg>


using namespace std;

class Vec{

  public:

    //Constructors
    Vec(int n=1, double a=0); //DONE
    Vec(int n, double *a); //DONE
    Vec(const Vec &v0); //DONE 
    //Destructor
    ~Vec(); //DONE

    Vec(int n, ...);

    //Functions
    void SetEntries (int n, double*); //DONE
    void Print(); //DONE
    int size(); //DONE
    int size() const; //DONE
    double dot(const Vec &v); //DONE
    void swap(int i1, int i2); //DONE

    //Overloaded Functions
    Vec& operator=(const Vec &v); //DONE
    Vec& operator+=(const Vec &v); //DONE
    Vec operator+(const Vec &v); //DONE
    Vec& operator-=(const Vec &v); //DONE
    Vec operator-(const Vec &v); //DONE
    double& operator[](int i); //DONE
    double operator[](int i) const; 
    Vec operator-(); //DONE
    Vec operator+(); //EPAH IDK
    Vec operator*(const Vec &v); //DONE
    Vec operator*(double c); //DONE
    Vec& operator*=(const Vec &v); //DONE
    Vec& operator*=(double c); //DONE

    Vec operator%(const Vec &v); //DONE <= CROSS PRODUCT

    //Getters and Setters
    double At(int i); //DONE
    void SetAt(int i, double a); //DONE


  private:
    int N;
    double *entries;

};

#endif
