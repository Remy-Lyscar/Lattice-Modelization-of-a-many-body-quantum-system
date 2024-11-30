#ifndef LISTS_H
#define LISTS_H

class Lists {
private:
	vector<double> list;

public :
	Lists(const vector<double>& values);
	void add(double x);
	void apply(double (*f)(double));
};
