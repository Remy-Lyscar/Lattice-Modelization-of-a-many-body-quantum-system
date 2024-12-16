#ifndef RANDOMVECTOR_H
#define RANDOMVECTOR_H

#include <Eigen/Dense>

enum class DistributionType {
    Uniform,
    Normal,
    Binomial,
    Exponential,
    Poisson,
    Gamma,
    Cauchy,
    Bernoulli
};

class RandomVector {
public:
    RandomVector(int size, int n, unsigned int seed = 42);

    Eigen::VectorXd uniform() const;
    Eigen::VectorXd normal() const;
    Eigen::VectorXd binomial(int t, double p) const;
    Eigen::VectorXd exponential(double lambda) const;
    Eigen::VectorXd poisson(double mean) const;
    Eigen::VectorXd gamma(double alpha, double beta) const;
    Eigen::VectorXd cauchy(double a, double b) const;
    Eigen::VectorXd bernoulli(double p) const;

private:
    int size;
    int n;
    unsigned int seed;
};

#endif // RANDOMVECTOR_H
