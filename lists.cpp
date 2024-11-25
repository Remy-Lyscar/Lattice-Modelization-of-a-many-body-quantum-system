#include <vector>
#include <iostream>
#include "lists.h"
using namespace std;

class Lists {
private:
	vector<double> list;
	
public:
	/* constructor */
	Lists(const vector<double>& values) : value(values) { };

    /* add the value x in ascending order in the list */
    void add(double x) {
        int i = 0;
        while (i < list.size() && list[i] < x) {
            ++i;
        }
        list.insert(list.begin() + i, x);
    }

    /*apply a function f to each value of the list*/
    void apply(double (*f)(double)) {
        for (int i = 0; i < list.size(); i++) {
            list[i] = f(list[i]);
        }
    }
