#pragma once
#include "Structures.h"
#include <algorithm>
#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;
/*Point Gradient(Point z, long double (*J)(Point)) {
    double dy = 0.000001;
    Point Grad(z.getDimension());
    Point newPoint;
    double curFitness = J(z);
    double newFitness;
    for (int i = 0; i < z.getDimension(); ++i) {
        newPoint = z;
        newPoint[i] = z[i]+dy;
        newFitness = J(newPoint);
        Grad[i] = (newFitness-curFitness) / dy;
        curFitness = newFitness;
    }
    return Grad;
}*/
Point Gradient(Point z, long double (*J)(Point), long double h = 1e-6) {
    int dim = z.getDimension();
    Point Grad(dim);
    long double fz = J(z);  // сохраняем значение в исходной точке

    for (int i = 0; i < dim; ++i) {
        Point zp = z;
        zp[i] += h;
        long double fzp = J(zp);
        Grad[i] = (fzp - fz) / h;
    }
    return Grad;
}

std::vector<double> Gradient_Search(double start_learning_rate, long double (*J)(Point),
    int dim, double a, double b, int max_iterations = 100) {
    Point z(dim);
    z.MakePoint(a, b);
    std::vector<double> Fitness;
    double prev_J = J(z);
    double curr_J = prev_J;
    int iter = 0;
    while (iter < max_iterations) {

        Point grad = Gradient(z, J);
        Point newPoint = z - grad * start_learning_rate;
        for (int i = 0; i < dim; ++i) {
            if (newPoint[i] < a) {
                newPoint[i] = a;
            }
            else if (newPoint[i] > b) {
                newPoint[i] = b;
            }
        }
        curr_J = J(newPoint);
        if (abs(prev_J - curr_J) < 0.00000000001) {
            while (iter < max_iterations) {
                Fitness.push_back(curr_J);
                iter += 1;

            }
            break;
        }
        if (curr_J < prev_J) {
            z = newPoint;
            prev_J = curr_J;
            Fitness.push_back(curr_J);
            iter += 1;

        }
        else {
            start_learning_rate *= 0.5;
        }

    }
    if (curr_J < 0) {
        for (int i = 0; i < dim; ++i) {
            std::cout << z.getPoint()[i] << ' ';
        }
        std::cout << std::endl;
    }
    return Fitness;
}

std::vector<double> Global_Gradient_Search(double start_learning_rate,
    long double (*J)(Point), int dim, double a, double b,
    int max_iterations_per_point,
    int Point_count) {
    std::vector<double> Result(max_iterations_per_point, INT_MAX);
    std::vector<double> Result1;
    for (int i = 0; i < Point_count; ++i) {
        Result1 = Gradient_Search(0.5, J, dim, a, b, max_iterations_per_point);
        for (int j = 0; j < max_iterations_per_point; ++j) {
            Result[j] = (Result[j] > Result1[j]) ? Result1[j] : Result[j];
        }
    }
    return Result;
}

void Middle_Gradient(
    int Point_count,

    int max_iterations,
    long double (*J)(Point),
    int dim,
    double a,
    double b,
    const std::string& plot_name,
    int num_runs = 10) {
    std::vector<double> Total(max_iterations, 0.0);
    std::vector<double> curResult;
    for (int i = 0; i < num_runs; ++i) {
        curResult = Global_Gradient_Search(0.5,
            J, dim, a, b,
            max_iterations,
            Point_count);
        for (int j = 0; j < max_iterations; ++j) {
            Total[j] += curResult[j];
        }
    }
    for (int j = 0; j < max_iterations; ++j) {
        Total[j] *= 0.1;
    }
    std::vector<int> Fitness_call;
    int call = max_iterations * (dim + 1) * Point_count;
    for (int i = 0; i < max_iterations; ++i) {

        Fitness_call.push_back(call);
        call += max_iterations * (dim + 1) * Point_count;
    }
    // Рисуем все запуски полупрозрачными
    /*for (int run = 0; run < num_runs; ++run) {
        plt::plot(point_numbers, all_histories[run],
            { {"alpha", "0.2"}, {"color", "gray"} });
    }*/

    // Рисуем среднюю линию поверх
    plt::named_plot(plot_name, Fitness_call, Total);

    plt::xlabel("Average call of Fitness Function");
    plt::ylabel("Average Fitness value");
    plt::legend();
    plt::grid(true);
    //plt::savefig(plot_name + ".png");
    //plt::show();
}
std::vector<long double> Monte_Karlo(int dim, int T, long double (*J)(Point), double a, double b) {

    long double Fitness_Best = INT_MAX;
    std::vector<long double> Fitness_list;
    for (int i = 0; i < T; ++i) {
        Point cur(dim);
        for (int j = 0; j < dim; ++j) {
            cur[j] = my_random_d(a, b);

        }
        if (J(cur) < Fitness_Best) {
            Fitness_Best = J(cur);

        }
        Fitness_list.push_back(Fitness_Best);
    }
    return Fitness_list;
}
void MiddleMK(int dim, int T, long double (*J)(Point), double a, double b, const std::string& plot_name) {
    std::vector<long double> M;
    Point middle(T);
    std::vector<int>Iteration;
    for (int i = 0; i < 10; ++i) {
        M = Monte_Karlo(dim, T, J, a, b);
        Iteration = {};
        middle = middle + M;
    }
    middle = middle * 0.1;
    for (int j = 0; j < T; ++j) {
        Iteration.push_back(j);
    }

    plt::named_plot(plot_name, Iteration, middle.getPoint());

}