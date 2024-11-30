#ifndef POLYNOM_H
#define POLYNOM_H
#include <vector>
#include <iostream>
#include <algorithm>
using namespace std;

class Polynom {
private:
	vector<double> coeff;
public:
	Polynom (const vector<double>& vector_coeff);
	Polynom operator + (const Polynom& Q) const;
	Polynom operator * (double scalar) const ;
	Polynom operator * (const Polynom& Q) const;
	double operator() (double x) const;
};

#endif