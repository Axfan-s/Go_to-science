#pragma once
#include "Structures.h"
#include <math.h>
#define pi 3.1415926535897932384626433832795028841971693
long double Rasstrigin(Point z) {
    double A = 10;
    double sum = 0;
    for (size_t i = 0; i < z.getDimension(); ++i) {
        sum += A - A * cos(2 * pi * z[i]) + z[i] * z[i];
    }
    return sum;
}

long double Griewank(Point z) {
    double sum = 0;
    double multiply = 1.0;

    for (size_t i = 0; i < z.getDimension() - 1; ++i) {
        sum += z[i + 1] * z[i + 1];
        double m = cos(z[i + 1] * 1.0 / sqrt(i + 1));
        multiply *= m;
    }
    return   sum * 1.0 / 4000 - 10 * multiply + 10;
}

long double Styblinsky_Tang(Point  z) {
    double A = 39.16617;
    double sum = 0;
    for (size_t i = 0; i < z.getDimension(); ++i) {
        sum += 0.5 * (z[i] * z[i] * z[i] * z[i] - 16 * z[i] * z[i] + 5 * z[i]);
    }
    sum += A * z.getDimension();
    return sum;
}

long double Levy(Point z) {
    double sum = 0;
    for (size_t i = 0; i < z.getDimension() - 1; ++i) {
        double w = 1 + (z[i] - 1) * 1.0 / 4;
        sum += (w - 1) * (w - 1) * (1 + 10 * sin(pi * w + 1) * sin(pi * w + 1));

    }
    double w = 1 + (z[z.getDimension() - 1] - 1) * 1.0 / 4;
    double sum_t = sin(pi * (1 + (z[0] - 1) * 1.0 / 4)) * sin(pi * (1 + (z[0] - 1) * 1.0 / 4)) + sum + (w - 1) * (w - 1) * (1 + sin(2 * pi * w) * sin(2 * pi * w));
    return sum_t;
}

long double Schwefel(Point x) {
    int d = x.getDimension();
    double sum = 0;
    for (int i = 0; i < d; ++i) {
        sum += x[i] * sin(sqrt(fabs(x[i])));
    }
    return 418.9829 * d - sum;
}

long double Ackley(Point x) {
    int d = x.getDimension();

    // Стандартные параметры
    const long double a = 20.0L;
    const long double b = 0.2L;
    const long double c = 2.0L * pi;  // M_PI из <cmath>

    long double sum_sq = 0.0L;   // сумма квадратов (для первого члена)
    long double sum_cos = 0.0L;  // сумма косинусов (для второго члена)

    for (int i = 0; i < d; ++i) {
        sum_sq += x[i] * x[i];
        sum_cos += cos(c * x[i]);
    }

    // Первый член: -a * exp(-b * sqrt(mean_sq))
    long double mean_sq = sum_sq / d;
    long double term1 = -a * exp(-b * sqrt(mean_sq));

    // Второй член: -exp(mean_cos)
    long double mean_cos = sum_cos / d;
    long double term2 = -exp(mean_cos);

    return term1 + term2 + a + exp(1.0L);
}
double CEC2017_1(Point z) {
    int D = z.getDimension();
    double f = 0.0;
    double penalty_coefficient = 1e4;
    // Целевая функция
    for (int i = 0; i < D; ++i) {
        double inner_sum = 0.0;
        for (int j = 0; j <= i; ++j) {
            inner_sum += z[j];
        }
        f += inner_sum * inner_sum;
    }

    // Ограничение
    double g = 0.0;
    for (int i = 0; i < D; ++i) {
        g += z[i] * z[i] - 5000.0 * std::cos(0.1 * pi * z[i]) - 4000.0;
    }

    // Штрафной подход
    double penalty = 0.0;
    if (g > 0) {
        penalty = penalty_coefficient * g;
    }

    return f + penalty;
}

double CEC2017_3(Point  z) {
    double sum;
    double sum2 = 0;
    double total;
    for (int i = 0; i < z.getDimension(); ++i) {
        sum = 0;
        for (int j = 0; j < i; ++j) {
            sum += z[j] * z[j];
        }
        sum2 += sum + abs(z[i] * std::sin(0.1 * pi * z[i]));
    }
    double sum3 = 0;
    for (int i = 0; i < z.getDimension(); ++i) {
        sum3 += z[i] * z[i] - 5000 * std::cos(0.1 * pi * z[i]) - 4000;
    }
    total = sum2;
    if (sum3 >= 0) {
        total += sum3;
    }
    return total;
}

double CEC2017_6(Point z) {
    //double sum;
    double sum2 = 0;
    double total;
    for (int i = 0; i < z.getDimension(); ++i) {
        sum2 += z[i] * z[i] - 10 * std::cos(2 * pi * z[i]) + 10 + 10 * abs(z[i] * (std::sin(z[i]) + std::cos(z[i])) + 10 * abs(z[i] * (std::sin(pi * z[i]) + std::cos(pi * z[i])))) + 10 * abs(std::sin(2 * sqrt(z[i])));
    }
    total = sum2;
    return total;
}

double CEC2017_11(Point z) {
    double sum2 = 0;
    double total;
    double mult = 1;
    for (int i = 0; i < z.getDimension() - 1; i++) {
        sum2 += 10 * (z[i] - z[i + 1]) * (z[i] - z[i + 1]);
    }
    for (int i = 0; i < z.getDimension(); i++) {
        mult *= z[i];
    }
    if (mult > 0) {
        sum2 += 10 * mult;
    }
    total = sum2;
    return total;
}

double CEC2017_9(Point z) {
    double sum2 = *std::max_element(begin(z.getPoint()), end(z.getPoint()));
    double total;
    double mult = 1;
    for (int i = 0; i < z.getDimension() / 2 - 1; i++) {
        sum2 += 10 * (z[2 * i] * z[2 * i] - z[2 * i + 1]) * (z[2 * i] * z[2 * i] - z[2 * i + 1]);
    }
    for (int i = 0; i < z.getDimension() / 2; i++) {
        mult *= z[2 * i];
    }
    if (mult > 0) {
        sum2 += 10 * mult;
    }
    total = sum2;
    return total;
}

double Phi(size_t k) {
    return sqrt(1 + k * 1.0 / 10000);
}
double Phi1(size_t k) {
    return 1 + log(1 + k * 1.0 / 1000) / 10000;
}
double Phi2(size_t k) {
    return 1 + log(1 + log(1 + k * 1.0 / 10000));
}
double Phi3(size_t k) {
    return sqrt(k * 1.0 / 10000);
}
double Sigma1(size_t k) {
    return 500000.0 / (k + 1);
}
double Sigma11(size_t k) {
    return 50000.0 / (k + 1);
}
double Sigma2(size_t k) {
    return 0.000000001 / sqrt(k + 1);
}
double Sigma3(size_t k) {
    return 100.0 / (log2(k + 1));
}
double Sigma4(size_t k) {
    return 500.0 / (log2(k + 1) * log2(k + 1));
}
double Sigma5(size_t k) {
    return 10000.0 * (log2(k + 1) * log2(k + 1)) / (k + 1);
}
double Sigma_SM(size_t k) {
    return pow(1.0 / k, 0.05);
}
double H(size_t k) {
    return (1 - tanh((k + 1) * 0.01));
}
double H1(size_t k) {
    return std::exp(-pow(k, 0.01) / 10000);
}
double H2(size_t k) {
    return 0.005 * (sin(0.07 * k) + 1) / sqrt(sqrt(0.07 * k+1));
}
double H3(size_t k) {
    return 0.0000005 * (sin(0.07 * k) + 1) / sqrt(sqrt(0.07 * k + 1));
}
long double g(Point z, Point shift, double global_min) {
    long double sum = 0;
    for (int i = 0; i < z.getDimension(); ++i) {
        sum += (z[i] - shift[i]) * (z[i] - shift[i]);
    }
    return sum + global_min;
}

/*double expLevy(Point z) {
    return 1.0 / (1 + exp(Levy(z) / z.getDimension()) * 1.0);
}

double expStyblinsky_Tang(Point z) {
    return 1.0 / (1 + exp(Styblinsky_Tang(z) / z.getDimension()) * 1.0);
}
double expGriewank(Point z) {
    return 1.0 / (1 + exp(Griewank(z) / z.getDimension()) * 1.0);
}
double expRasstrigin(Point z) {
    return 1.0 / (1 + exp(Rasstrigin(z) / z.getDimension()) * 1.0);
}

double expCEC2017_1(Point z) {
    return 1.0 / (1 + exp(CEC2017_1(z) / z.getDimension()) * 1.0);
}

double expCEC2017_3(Point z) {
    return 1.0 / (1 + exp(CEC2017_3(z) / z.getDimension()) * 1.0);
}

double expCEC2017_6(Point z) {
    return 1.0 / (1 + exp(CEC2017_6(z) / z.getDimension()) * 1.0);
}
double expCEC2017_9(Point z) {
    return 1.0 / (1 + exp(CEC2017_9(z) / z.getDimension()) * 1.0);
}
double expCEC2017_11(Point z) {
    return 1.0 / (1 + exp(CEC2017_11(z) / z.getDimension()) * 1.0);
}*/
long double exponential_Fitness(double F, int dim) {
    return 1.0 / (1 + exp(F / dim) * 1.0);
}

    class GKLS {
    private:
        int n;                          // размерность
        int m;                          // количество минимумов (включая M0)
        long double a, b;               // границы области
        long double t = 0;                  // значение в M0
        std::vector<Point> M;           // центры минимумов
        std::vector<long double> Rho;   // радиусы шаров (для M0 — 0)
        std::vector<long double> F;
    public:
        GKLS(int dim, int num_minima, long double left, long double right, long double t_value = 0.0)
            : n(dim), m(num_minima), a(left), b(right), t(t_value) {

            if (m < 2) {
                throw std::invalid_argument("GKLS: num_minima must be at least 2");
            }

            // M0 — центр параболоида
            Point T(n);
            for (int i = 0; i < n; ++i) {
                T[i] = (a + b) * 0.5;
            }
            M.push_back(T);
            Rho.push_back(0.0);  // у M0 нет радиуса
            F.push_back(t);
            generateDomains();
        }

        /*void generateDomains() {
            const long double EPS = 1e-9;
            long double max_radius = (std::abs(b) + std::abs(a)) * 0.1;
            const int MAX_ATTEMPTS = 10000;

            for (int i = 1; i < m; ++i) {
                bool success = false;
                F.push_back(my_random_d(t, t + 1));
                for (int attempt = 0; attempt < MAX_ATTEMPTS; ++attempt) {
                    // Генерируем центр и радиус
                    Point Q(n);
                    long double new_rho = my_random_d(0.01, max_radius);
                    for (int j = 0; j < n; ++j) {
                        Q[j] = my_random_d(a, b);
                    }

                    // Проверяем пересечение с существующими шарами
                    bool has_conflict = false;
                    for (size_t k = 0; k < M.size(); ++k) {
                        long double dist = distance(M[k], Q);
                        if (dist < Rho[k] + new_rho + EPS) {
                            has_conflict = true;
                            break;
                        }
                    }

                    // Если конфликта нет — сохраняем шар и выходим
                    if (!has_conflict) {
                        M.push_back(Q);
                        Rho.push_back(new_rho);
                        success = true;
                        break;
                    }
                }

                if (!success) {
                    throw std::runtime_error("GKLS: cannot generate non-overlapping balls");
                }
            }
        }*/
        void generateDomains() {
            const long double EPS = 1e-9;
            long double max_radius = (std::abs(b) + std::abs(a)) * 0.5;
            const int MAX_ATTEMPTS = 10000;

            for (int i = 1; i < m; ++i) {
                bool success = false;
                F.push_back(my_random_d(t + 3, t + 7));
                // Генерируем центр и радиус
                Point Q(n);
                long double new_rho = my_random_d(0.01, max_radius);
                for (int j = 0; j < n; ++j) {
                    Q[j] = my_random_d(a, b);
                }
                M.push_back(Q);
            }

            for (int i = 1; i < m; ++i) {
                long double Rho_min = 1000000000;
                for (int j = 0; j < m; ++j) {
                    if (i != j) {
                        long double D = distance(M[i], M[j]);
                        if (D < Rho_min) {
                            Rho_min = D;
                        }
                    }
                    else {
                        continue;
                    }
                }
                Rho.push_back(0.5 * Rho_min);
            }
        }

        long double GKLS_example(Point x) {
            long double y;
            bool gotcha = false;
            for (int i = 1; i < m; ++i) {

                if (distance(x, M[i]) < Rho[i]) {
                    Point kusok = x - M[i];
                    Point kusok2 = M[0] - M[i];
                    long double norma = kusok.norm();
                    long double A = kusok2.norm() * kusok2.norm() + t - F[i];
                    long double Sc = ScalarMultiply(x - M[i], M[0] - M[i]);
                    long double Scal_Mult = (1.0 - Sc * (2.0 / Rho[i])) * norma;
                    y = Scal_Mult + (1.0 / (Rho[i] * Rho[i])) * A * norma * norma + F[i];
                    gotcha = true;
                }
            }
            if (!gotcha) {
                Point kusok = x - M[0];
                long double norma = kusok.norm();
                y = norma * norma + t;
            }
            return y;
        }
};