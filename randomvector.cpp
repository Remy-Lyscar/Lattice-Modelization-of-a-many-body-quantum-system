#include <Eigen/Dense>
#include <random>

#include "randomvector.h"

RandomVector::RandomVector(int size, int n, unsigned int seed) : size(size), n(n), seed(seed) {}

Eigen::VectorXd RandomVector::uniform() const {
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    Eigen::VectorXd v(size);
    for (int i = 0; i < size; i++) {
        v(i) = dist(gen);
    }
    double sum = v.sum();
    v = (v / sum) * n;
    return v;
}

Eigen::VectorXd RandomVector::normal() const {
    std::mt19937 gen(seed);
    std::normal_distribution<double> dist(0.0, 1.0);
    Eigen::VectorXd v(size);
    for (int i = 0; i < size; i++) {
        v(i) = dist(gen);
    }
    double sum = v.sum();
    v = (v / sum) * n;
    return v;
}

Eigen::VectorXd RandomVector::binomial(int t, double p) const {
    std::mt19937 gen(seed);
    std::binomial_distribution<int> dist(t, p);
    Eigen::VectorXd v(size);
    for (int i = 0; i < size; i++) {
        v(i) = dist(gen);
    }
    double sum = v.sum();
    v = (v / sum) * n;
    return v;
}

Eigen::VectorXd RandomVector::exponential(double lambda) const {
    std::mt19937 gen(seed);
    std::exponential_distribution<double> dist(lambda);
    Eigen::VectorXd v(size);
    for (int i = 0; i < size; i++) {
        v(i) = dist(gen);
    }
    double sum = v.sum();
    v = (v / sum) * n;
    return v;
}

Eigen::VectorXd RandomVector::poisson(double mean) const {
    std::mt19937 gen(seed);
    std::poisson_distribution<int> dist(mean);
    Eigen::VectorXd v(size);
    for (int i = 0; i < size; i++) {
        v(i) = dist(gen);
    }
    double sum = v.sum();
    v = (v / sum) * n;
    return v;
}

Eigen::VectorXd RandomVector::gamma(double alpha, double beta) const {
    std::mt19937 gen(seed);
    std::gamma_distribution<double> dist(alpha, beta);
    Eigen::VectorXd v(size);
    for (int i = 0; i < size; i++) {
        v(i) = dist(gen);
    }
    double sum = v.sum();
    v = (v / sum) * n;
    return v;
}

Eigen::VectorXd RandomVector::cauchy(double a, double b) const {
    std::mt19937 gen(seed);
    std::cauchy_distribution<double> dist(a, b);
    Eigen::VectorXd v(size);
    for (int i = 0; i < size; i++) {
        v(i) = dist(gen);
    }
    double sum = v.sum();
    v = (v / sum) * n;
    return v;
}

Eigen::VectorXd RandomVector::bernoulli(double p) const {
    std::mt19937 gen(seed);
    std::bernoulli_distribution dist(p);
    Eigen::VectorXd v(size);
    for (int i = 0; i < size; i++) {
        v(i) = dist(gen);
    }
    double sum = v.sum();
    v = (v / sum) * n;
    return v;
}
