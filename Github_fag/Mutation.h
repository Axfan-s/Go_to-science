#pragma once
#include <optional>
#include <variant>
#include <random>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/special_functions/detail/erf_inv.hpp>
#include "Structures.h"



long double Normal(long double gene, double a, double b, double Var) {
    long double mut_gene;
    double q = my_random_d(0, 1);
    double eta = my_random_d(0, 1);
    long double s = boost::math::erf(b - gene);
    //std::cout << s<<" ";
    long double Expr = eta * (boost::math::erf((b - gene) / sqrt(2 * Var))) + (1 - eta) * boost::math::erf((a - gene) / sqrt(2 * Var));
    long double R;
    if (Expr >= 1) {
        R = b;
    }
    else if (Expr <= -1) {
        R = a;
    }
    else {
        R = sqrt(2 * Var) * (boost::math::erf_inv(Expr)) + gene;
    }

    mut_gene = R;
    return R;
}

long double Cauchy(long double gene, double a, double b, double Var) {
    long double mut_gene;
    double q = my_random_d(0, 1);
    double eta = my_random_d(0, 1);
    if (gene < a) {
        gene = a;
    }
    else if (gene > b) {
        gene = b;
    }
    double eL = atan((a - gene) / Var);
    double eR = atan((b - gene) / Var);
    double A = eR - eL;
    long double B = 1.0 / (eR - eL);
    long double R = Var * tan((eta * A + atan((a - gene) / Var))) + gene;

    if (R < a) {
        R = a;
    }
    if (R > b) {
        R = b;
    }
    mut_gene = R;
    return mut_gene;
}