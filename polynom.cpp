#include <vector>
#include <iostream>
#include <algorithm>
#include "polynom.h"
using namespace std;

class Polynom{
private :
    vector<double> coeff;

public:
    /* constructor */
    Polynom(const vector<double>& vector_coeff) : coeff(vector_coeff) { };

    /* add to a polynom a polynom Q */
    Polynom operator + (const Polynom & Q) const{
        int degree = max(this->coeff.size(), Q.coeff.size());
        vector<double> result(degree,0.0);

        for (int i=0; i<this->coeff.size(); i++){
            result[i] = this->coeff[i];
        }
        for (int j=0; j<Q.coeff.size(); j++){
            result[j] += Q.coeff[j];
        }
        return Polynom(result);
    }

    /* multiply a polynom by a scalar */
   Polynom operator * (double scalar) const{
       vector<double> result = coeff;
        for (int i = 0; i < result.size(); i++){
            result[i] *= scalar;
        }
        return Polynom(result);
    }

    /* multiply a polynom by a polynom Q */
    Polynom operator * (const Polynom & Q) const {
        int degree = this->coeff.size() + Q.coeff.size() -1 ;
        vector<double> result(degree,0.0);
        for (int i = 0; i < this->coeff.size(); i++){
            for (int j = 0; j < Q.coeff.size(); j++){
                result[i + j] += this->coeff[i] * Q.coeff[j];
            }
        }
        return Polynom(result);
    }
    
    /* evaluate a polynom at a point x using the Horner's method */
    double operator() (double x) const{
    double result = 0.0;
    for (int i = coeff.size() - 1; i >= 0; i--) {
    result = result * x + coeff[i];
    }
    return result;
}   



