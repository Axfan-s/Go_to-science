#pragma once
#include "Structures.h"
#include "matplotlibcpp.h"
#include <cmath>
#include <limits>
#include <algorithm>
#include <random>
namespace plt = matplotlibcpp;
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<double> dist(0.0, 1.0);

double my_rand2() {
    static std::uniform_int_distribution<int> int_dist(0, 1000000);
    return int_dist(gen);
}

// Âñïîìîãàòåëüíûå ôóíêöèè äëÿ ðàáîòû ñ ãðàíèöàìè
Point createBoundMask(const Point& point, double bound, bool greater) {
    Point mask(point.getDimension());
    for (size_t i = 0; i < point.getDimension(); ++i) {
        if (greater) {
            mask[i] = (point[i] > bound) ? 1.0 : 0.0;
        }
        else {
            mask[i] = (point[i] < bound) ? 1.0 : 0.0;
        }
    }
    return mask;
}

std::vector<long double> SnakeOptimization(
    size_t N,
    size_t T,
    long double (*objective_func)(Point),
    size_t dim,
    double lb,
    double ub,
    std::vector<int>& iteration_counts,
    bool time_flag = false)
{
    const double eps = std::numeric_limits<double>::epsilon();
    const std::vector<int> vec_flag = { 1, -1 };
    const double threshold = 0.25;
    const double threshold2 = 0.6;
    const double C1 = 0.5;
    const double C2 = 0.05;
    const double C3 = 2.0;

    int total_evaluations = 0;
    std::vector<Point> Record;
    // 1. Èíèöèàëèçàöèÿ ïîïóëÿöèè
    Population X(N, dim);
    X.GeneratePopulation(lb, ub);

    std::vector<long double> fitness(N);
    for (size_t i = 0; i < N; ++i) {
        fitness[i] = objective_func(X[i]);
    }
    total_evaluations += N;

    // Íàõîæäåíèå ãëîáàëüíîãî ëó÷øåãî
    size_t global_best_idx = 0;
    long double global_best_fit = fitness[0];
    for (size_t i = 1; i < N; ++i) {
        if (fitness[i] < global_best_fit) {
            global_best_fit = fitness[i];
            global_best_idx = i;
        }
    }

    Point X_food = X[global_best_idx];

    // 2. Ðàçäåëåíèå íà ìóæñêèå è æåíñêèå îñîáè
    size_t Nm = N / 2;        // Ìóæñêèå îñîáè
    size_t Nf = N - Nm;       // Æåíñêèå îñîáè

    Population Xm(Nm, dim);
    std::vector<double> fitness_m(Nm);
    for (size_t i = 0; i < Nm; ++i) {
        Xm[i] = X[i];
        fitness_m[i] = fitness[i];
    }

    Population Xf(Nf, dim);
    std::vector<long double> fitness_f(Nf);
    for (size_t i = 0; i < Nf; ++i) {
        Xf[i] = X[Nm + i];
        fitness_f[i] = fitness[Nm + i];
    }

    // Íàõîæäåíèå ëó÷øèõ â êàæäîé ãðóïïå
    size_t best_m_idx = 0;
    long double best_m_fit = fitness_m[0];
    size_t worst_m_idx = 0;
    long double worst_m_fit = fitness_m[0];

    for (size_t i = 1; i < Nm; ++i) {
        if (fitness_m[i] < best_m_fit) {
            best_m_fit = fitness_m[i];
            best_m_idx = i;
        }
        if (fitness_m[i] > worst_m_fit) {
            worst_m_fit = fitness_m[i];
            worst_m_idx = i;
        }
    }

    size_t best_f_idx = 0;
    long double best_f_fit = fitness_f[0];
    size_t worst_f_idx = 0;
    long double worst_f_fit = fitness_f[0];

    for (size_t i = 1; i < Nf; ++i) {
        if (fitness_f[i] < best_f_fit) {
            best_f_fit = fitness_f[i];
            best_f_idx = i;
        }
        if (fitness_f[i] > worst_f_fit) {
            worst_f_fit = fitness_f[i];
            worst_f_idx = i;
        }
    }

    Point X_best_m = Xm[best_m_idx];
    Point X_best_f = Xf[best_f_idx];

    std::vector<long double> best_history(T);

    // 3. Îñíîâíîé öèêë îïòèìèçàöèè
    int start = clock();
    for (size_t t = 0; t < T; ++t) {
        // Ðàñ÷åò ïàðàìåòðîâ òåìïåðàòóðû è Q
        long double temp = std::exp(-static_cast<long double>(t) / T);
        long double Q = C1 * std::exp((static_cast<long double>(t) - T) / T);
        if (Q > 1.0) Q = 1.0;

        // Ñîçäàíèå íîâûõ ïîïóëÿöèé
        Population X_new_m(Nm, dim);
        Population X_new_f(Nf, dim);

        if (Q < threshold) {
            // ========== ÔÀÇÀ ÈÑÑËÅÄÎÂÀÍÈß (áåç ïèùè) ==========

            // Ìóæñêèå îñîáè
            for (size_t i = 0; i < Nm; ++i) {
                size_t rand_idx = static_cast<size_t>(my_random_i(Nm - 1));


                int flag = vec_flag[my_random_i(1)];

                long double Am = std::exp(-fitness_m[rand_idx] / (fitness_m[i] + eps));

                for (size_t d = 0; d < dim; ++d) {
                    long double r = dist(gen);
                    X_new_m[i][d] = Xm[rand_idx][d] +
                        flag * C2 * Am * ((ub - lb) * r + lb);
                }
            }

            // Æåíñêèå îñîáè
            for (size_t i = 0; i < Nf; ++i) {

                size_t rand_idx = static_cast<size_t>(my_random_i(Nf - 1));
                //size_t rand_idx = static_cast<size_t>(my_random_d(lb, ub));

                int flag = vec_flag[my_random_i(1)];

                long double Af = std::exp(-fitness_f[rand_idx] / (fitness_f[i] + eps));

                for (size_t d = 0; d < dim; ++d) {
                    long double r = dist(gen);
                    X_new_f[i][d] = Xf[rand_idx][d] +
                        flag * C2 * Af * ((ub - lb) * r + lb);
                }
            }
        }
        else {
            // ========== ÔÀÇÀ ÝÊÑÏËÓÀÒÀÖÈÈ (åñòü ïèùà) ==========

            if (temp > threshold2) {
                // Ãîðÿ÷àÿ ñðåäà
                for (size_t i = 0; i < Nm; ++i) {
                    int flag = vec_flag[my_random_i(0, 1000) % 2];
                    for (size_t d = 0; d < dim; ++d) {
                        long double r = dist(gen);
                        X_new_m[i][d] = X_food[d] +
                            C3 * flag * temp * r * (X_food[d] - Xm[i][d]);
                    }
                }

                for (size_t i = 0; i < Nf; ++i) {
                    int flag = vec_flag[my_random_i(0, 1000) % 2];
                    for (size_t d = 0; d < dim; ++d) {
                        long double r = dist(gen);
                        X_new_f[i][d] = X_food[d] +
                            flag * C3 * temp * r * (X_food[d] - Xf[i][d]);
                    }
                }
            }
            else {
                // Õîëîäíàÿ ñðåäà
                if (dist(gen) > 0.6) {
                    // Ðåæèì áîðüáû
                    for (size_t i = 0; i < Nm; ++i) {
                        for (size_t d = 0; d < dim; ++d) {
                            long double r = dist(gen);
                            long double FM = std::exp(-best_f_fit / (fitness_m[i] + eps));
                            X_new_m[i][d] = Xm[i][d] +
                                C3 * FM * r * (Q * X_best_f[d] - Xm[i][d]);
                        }
                    }

                    for (size_t i = 0; i < Nf; ++i) {
                        for (size_t d = 0; d < dim; ++d) {
                            long double r = dist(gen);
                            long double FF = std::exp(-best_m_fit / (fitness_f[i] + eps));
                            X_new_f[i][d] = Xf[i][d] +
                                C3 * FF * r * (Q * X_best_m[d] - Xf[i][d]);
                        }
                    }
                }
                else {
                    // Ðåæèì ñïàðèâàíèÿ
                    for (size_t i = 0; i < Nm; ++i) {
                        for (size_t d = 0; d < dim; ++d) {
                            long double r = dist(gen);
                            // Áåçîïàñíûé äîñòóï ê æåíñêèì îñîáÿì
                            size_t f_idx = (i < Nf) ? i : Nf - 1;
                            long double Mm = std::exp(-fitness_f[f_idx] / (fitness_m[i] + eps));
                            X_new_m[i][d] = Xm[i][d] +
                                C3 * r * Mm * (Q * Xf[f_idx][d] - Xm[i][d]);
                        }
                    }

                    for (size_t i = 0; i < Nf; ++i) {
                        for (size_t d = 0; d < dim; ++d) {
                            long double r = dist(gen);
                            // Áåçîïàñíûé äîñòóï ê ìóæñêèì îñîáÿì
                            size_t m_idx = (i < Nm) ? i : Nm - 1;
                            int flag = vec_flag[my_random_i(0, 1000) % 2];
                            long double Mf = std::exp(-fitness_m[m_idx] / (fitness_f[i] + eps));
                            X_new_f[i][d] = Xf[i][d] +
                                C3 * flag * r * Mf * (Q * Xm[m_idx][d] - Xf[i][d]);
                        }
                    }

                    // Âîçìîæíîñòü îòêëàäêè ÿèö (ñáðîñ õóäøèõ îñîáåé)
                    if (vec_flag[my_random_i(0, 1000) % 2] == 1) {
                        for (size_t d = 0; d < dim; ++d) {
                            long double r = dist(gen);
                            X_new_m[worst_m_idx][d] = lb + r * (ub - lb);
                            X_new_f[worst_f_idx][d] = lb + r * (ub - lb);
                        }
                    }
                }
            }
        }

        // 4. Ïðîâåðêà ãðàíèö è îáíîâëåíèå ïîïóëÿöèé

        // Ìóæñêèå îñîáè
        for (size_t i = 0; i < Nm; ++i) {
            // Ïðîâåðêà ãðàíèö
            for (size_t d = 0; d < dim; ++d) {
                if (X_new_m[i][d] > ub) X_new_m[i][d] = ub;
                if (X_new_m[i][d] < lb) X_new_m[i][d] = lb;
            }

            // Îöåíêà íîâîé ïîçèöèè
            long double new_fit = objective_func(X_new_m[i]);
            total_evaluations++;

            if (new_fit < fitness_m[i]) {
                fitness_m[i] = new_fit;
                Xm[i] = X_new_m[i];

                // Îáíîâëåíèå ëó÷øåé ìóæñêîé îñîáè
                if (new_fit < best_m_fit) {
                    best_m_fit = new_fit;
                    best_m_idx = i;
                    X_best_m = Xm[i];
                }

                // Îáíîâëåíèå õóäøåé ìóæñêîé îñîáè
                if (new_fit > worst_m_fit) {
                    worst_m_fit = new_fit;
                    worst_m_idx = i;
                }
            }
        }

        // Æåíñêèå îñîáè
        for (size_t i = 0; i < Nf; ++i) {
            // Ïðîâåðêà ãðàíèö
            for (size_t d = 0; d < dim; ++d) {
                if (X_new_f[i][d] > ub) X_new_f[i][d] = ub;
                if (X_new_f[i][d] < lb) X_new_f[i][d] = lb;
            }

            // Îöåíêà íîâîé ïîçèöèè
            long double new_fit = objective_func(X_new_f[i]);
            total_evaluations++;

            if (new_fit < fitness_f[i]) {
                fitness_f[i] = new_fit;
                Xf[i] = X_new_f[i];

                // Îáíîâëåíèå ëó÷øåé æåíñêîé îñîáè
                if (new_fit < best_f_fit) {
                    best_f_fit = new_fit;
                    best_f_idx = i;
                    X_best_f = Xf[i];
                }

                // Îáíîâëåíèå õóäøåé æåíñêîé îñîáè
                if (new_fit > worst_f_fit) {
                    worst_f_fit = new_fit;
                    worst_f_idx = i;
                }
            }
        }

        // 5. Îáíîâëåíèå ãëîáàëüíîãî ëó÷øåãî ðåøåíèÿ
        if (best_m_fit < best_f_fit) {
            global_best_fit = best_m_fit;
            X_food = X_best_m;
        }
        else {
            global_best_fit = best_f_fit;
            X_food = X_best_f;
        }

        // Ñîõðàíåíèå èñòîðèè

        best_history[t] = global_best_fit;
        if (time_flag) {
            int now = clock() - start;
            iteration_counts.push_back(now);
        }
        else {
            iteration_counts.push_back(total_evaluations);
        }
        Record.push_back(X_food);
    }
    /* RecordResult(Record);
     // Âûâîä ðåçóëüòàòà
     std::cout << "Ëó÷øåå ðåøåíèå: ";
     for (size_t i = 0; i < dim; ++i) {
         std::cout << X_food[i] << " ";
     }
     std::cout << "\nËó÷øåå çíà÷åíèå ôóíêöèè: " << global_best_fit << std::endl;
     std::cout << "Âñåãî îöåíîê ôóíêöèè: " << total_evaluations << std::endl;*/
    for (int i = 0; i < X_food.getDimension(); ++i) {
        std::cout << X_food[i] << " ";
    }
    std::cout << std::endl;

    return best_history;
}

// Функция для многократного запуска и усреднения
void Middle(
    size_t N,
    size_t T,
    long double (*objective_func)(Point),
    size_t dim,
    double lb,
    double ub,
    const std::string& plot_name,
    bool time_flag = false,
    size_t num_trials = 10)
{
    std::vector<int> iteration_counts;
    std::vector<double> average_history(T, 0.0);
    std::vector<std::vector<long double>> Table;
    for (size_t trial = 0; trial < num_trials; ++trial) {
        std::vector<int> current_iterations;
        std::vector<long double> history = SnakeOptimization(
            N, T, objective_func, dim, lb, ub, current_iterations, time_flag);
        Table.push_back(history);
        // Если это первый запуск, сохраняем счетчик итераций
        if (trial == 0) {
            iteration_counts = current_iterations;
        }

        // Суммируем историю
        for (size_t i = 0; i < T; ++i) {
            average_history[i] += log(history[i] + 1);
        }
    }

    // Усреднение
    for (size_t i = 0; i < T; ++i) {
        average_history[i] /= num_trials;
    }
    saveSetVectorsToCSV(Table, "C:\\Users\\troit\\OneDrive\\Рабочий стол\\Science_result\\SO.csv", "SO");
    // Построение графика
    plt::named_plot(plot_name, iteration_counts, average_history);
    plt::xlabel("number of calls Fitness function");
    plt::ylabel("Best_Fitness");
    //plt::title("Сходимость Snake Optimization");
    plt::legend();
    plt::grid(true);
}

void Middle_parallel(
    size_t N,
    size_t T,
    long double (*objective_func)(Point),
    size_t dim,
    double lb,
    double ub,
    const std::string& plot_name,
    bool time_flag,
    size_t num_trials = 10)
{
    std::vector<int> iteration_counts;
    std::vector<double> average_history(T, 0.0);
    std::vector<std::vector<long double>> Table(num_trials);
    std::vector<std::vector<int>> all_iterations(num_trials);

    std::cout << "Запуск " << num_trials << " параллельных испытаний SO..." << std::endl;

    // Счётчик завершённых задач
    volatile int completed = 0;

#pragma omp parallel for
    for (int trial = 0; trial < static_cast<int>(num_trials); ++trial) {
        std::vector<int> current_iterations;
        std::vector<long double> history = SnakeOptimization(
            N, T, objective_func, dim, lb, ub, current_iterations/*,time_flag*/);

#pragma omp critical
        {
            Table[trial] = history;
            all_iterations[trial] = current_iterations;
            completed++;

            // Простой прогресс-бар
            int percent = (completed * 100) / num_trials;
            std::cout << "\rПрогресс: " << percent << "% (" << completed
                << "/" << num_trials << ")" << std::flush;
        }
    }
    std::cout << std::endl << "Все испытания завершены!" << std::endl;

    // Сохраняем итерации первого запуска
    if (num_trials > 0) {
        iteration_counts = all_iterations[0];
    }

    // Суммирование и усреднение
    for (size_t trial = 0; trial < num_trials; ++trial) {
        for (size_t i = 0; i < T && i < Table[trial].size(); ++i) {
            average_history[i] += Table[trial][i];
        }
    }

    for (size_t i = 0; i < T; ++i) {
        average_history[i] /= num_trials;
    }

    // Сохранение
    std::string filename = "C:\\Users\\troit\\OneDrive\\Рабочий стол\\Science_result\\SO_parallel.csv";
    saveSetVectorsToCSV(Table, filename.c_str(), "SO");
    std::cout << "Результаты сохранены в " << filename << std::endl;

    // График
    plt::named_plot(plot_name, iteration_counts, average_history);
    plt::xlabel("number of calls Fitness function");
    plt::ylabel("Best_Fitness");
    plt::legend();
    plt::grid(true);
}



std::vector<long double> SnakeOptimization2(
    size_t N,
    size_t T,
    long double (*objective_func)(Point),
    size_t dim,
    long double lb,
    long double ub,
    std::vector<int>& iteration_counts)
{
    const double eps = 1e-12;
    const std::vector<int> vec_flag = { 1, -1 };
    const double threshold = 0.25;
    //const double threshold2 = 0.6; //Статья
    const double threshold2 = 0.6;
    const double C1 = 0.5;
    const double C2 = 0.05;
    const double C3 = 2.0;

    int total_evaluations = 0;

    // 1. Инициализация популяции
    Population X(N, dim);
    X.GeneratePopulation(lb, ub);

    std::vector<long double> fitness(N);
    for (size_t i = 0; i < N; ++i) {
        fitness[i] = objective_func(X[i]);
    }
    total_evaluations += N;

    // Нахождение глобального лучшего
    size_t global_best_idx = 0;
    long double global_best_fit = fitness[0];
    for (size_t i = 1; i < N; ++i) {
        if (fitness[i] < global_best_fit) {
            global_best_fit = fitness[i];
            global_best_idx = i;
        }
    }
    Point X_food = X[global_best_idx];

    // 2. Разделение на мужские и женские особи
    size_t Nm = N / 2;
    size_t Nf = N - Nm;

    Population Xm(Nm, dim);
    std::vector<long double> fitness_m(Nm);
    for (size_t i = 0; i < Nm; ++i) {
        Xm[i] = X[i];
        fitness_m[i] = fitness[i];
    }

    Population Xf(Nf, dim);
    std::vector<long double> fitness_f(Nf);
    for (size_t i = 0; i < Nf; ++i) {
        Xf[i] = X[Nm + i];
        fitness_f[i] = fitness[Nm + i];
    }

    // Нахождение лучших и худших
    size_t best_m_idx = 0, worst_m_idx = 0;
    long double best_m_fit = fitness_m[0], worst_m_fit = fitness_m[0];
    for (size_t i = 1; i < Nm; ++i) {
        if (fitness_m[i] < best_m_fit) { best_m_fit = fitness_m[i]; best_m_idx = i; }
        if (fitness_m[i] > worst_m_fit) { worst_m_fit = fitness_m[i]; worst_m_idx = i; }
    }

    size_t best_f_idx = 0, worst_f_idx = 0;
    long double best_f_fit = fitness_f[0], worst_f_fit = fitness_f[0];
    for (size_t i = 1; i < Nf; ++i) {
        if (fitness_f[i] < best_f_fit) { best_f_fit = fitness_f[i]; best_f_idx = i; }
        if (fitness_f[i] > worst_f_fit) { worst_f_fit = fitness_f[i]; worst_f_idx = i; }
    }

    Point X_best_m = Xm[best_m_idx];
    Point X_best_f = Xf[best_f_idx];

    std::vector<long double> best_history(T);

    // 3. Основной цикл
    for (size_t t = 0; t < T; ++t) {
        long double temp = std::exp(-static_cast<long double>(t) / T);
        long double Q = C1 * std::exp((static_cast<long double>(t) - T) / T);
        if (Q > 1.0) Q = 1.0;

        Population X_new_m(Nm, dim);
        Population X_new_f(Nf, dim);

        if (Q < threshold) {
            // ========== ФАЗА ИССЛЕДОВАНИЯ ==========
            for (size_t i = 0; i < Nm; ++i) {
                size_t rand_idx = my_random_i(0, Nm - 1);

                int flag = vec_flag[my_random_i(0, 2)];
                long double Am = std::exp(-fitness_m[rand_idx] / (fitness_m[i] + eps));

                for (size_t d = 0; d < dim; ++d) {
                    long double r = my_random_d(0, 1);
                    long double val = Xm[rand_idx][d] + flag * C2 * Am * ((ub - lb) * r + lb);
                    X_new_m[i][d] = std::clamp(val, lb, ub);
                }
            }

            for (size_t i = 0; i < Nf; ++i) {
                size_t rand_idx = my_random_i(0, Nf - 1);
                int flag = vec_flag[my_random_i(0, 2)];
                long double Af = std::exp(-fitness_f[rand_idx] / (fitness_f[i] + eps));

                for (size_t d = 0; d < dim; ++d) {
                    long double r = my_random_d(0, 1);
                    long double val = Xf[rand_idx][d] + flag * C2 * Af * ((ub - lb) * r + lb);
                    X_new_f[i][d] = std::clamp(val, lb, ub);
                }
            }
        }
        else {
            // ========== ФАЗА ЭКСПЛУАТАЦИИ ==========
            if (temp > threshold2) {
                // Горячая среда
                for (size_t i = 0; i < Nm; ++i) {
                    int flag = vec_flag[my_random_i(0, 2)];
                    for (size_t d = 0; d < dim; ++d) {
                        long double r = my_random_d(0, 1);
                        long double val = X_food[d] + flag * C3 * temp * r * (X_food[d] - Xm[i][d]);
                        X_new_m[i][d] = std::clamp(val, lb, ub);
                    }
                }

                for (size_t i = 0; i < Nf; ++i) {
                    int flag = vec_flag[my_random_i(0, 2)];
                    for (size_t d = 0; d < dim; ++d) {
                        long double r = my_random_d(0, 1);
                        long double val = X_food[d] + flag * C3 * temp * r * (X_food[d] - Xf[i][d]);
                        X_new_f[i][d] = std::clamp(val, lb, ub);
                    }
                }
            }
            else {
                // Холодная среда
                if (my_random_d(0, 1) > 0.6) {
                    // Режим борьбы
                    for (size_t i = 0; i < Nm; ++i) {
                        for (size_t d = 0; d < dim; ++d) {
                            long double r = my_random_d(0, 1);
                            long double FM = std::exp(-best_f_fit / (fitness_m[i] + eps));
                            long double val = Xm[i][d] + C3 * FM * r * (Q * X_best_f[d] - Xm[i][d]);
                            X_new_m[i][d] = std::clamp(val, lb, ub);
                        }
                    }

                    for (size_t i = 0; i < Nf; ++i) {
                        for (size_t d = 0; d < dim; ++d) {
                            long double r = my_random_d(0, 1);
                            long double FF = std::exp(-best_m_fit / (fitness_f[i] + eps));
                            long double val = Xf[i][d] + C3 * FF * r * (Q * X_best_m[d] - Xf[i][d]);
                            X_new_f[i][d] = std::clamp(val, lb, ub);
                        }
                    }
                }
                else {
                    // Режим спаривания
                    for (size_t i = 0; i < Nm; ++i) {
                        size_t f_idx = (i < Nf) ? i : my_random_i(0, Nf);
                        for (size_t d = 0; d < dim; ++d) {
                            long double r = my_random_d(0, 1);
                            long double Mm = std::exp(-fitness_f[f_idx] / (fitness_m[i] + eps));
                            long double val = Xm[i][d] + C3 * r * Mm * (Q * Xf[f_idx][d] - Xm[i][d]);
                            X_new_m[i][d] = std::clamp(val, lb, ub);
                        }
                    }

                    for (size_t i = 0; i < Nf; ++i) {
                        size_t m_idx = (i < Nm) ? i : my_random_i(0, Nm);
                        for (size_t d = 0; d < dim; ++d) {
                            long double r = my_random_d(0, 1);
                            long double Mf = std::exp(-fitness_m[m_idx] / (fitness_f[i] + eps));
                            long double val = Xf[i][d] + C3 * r * Mf * (Q * Xm[m_idx][d] - Xf[i][d]);
                            X_new_f[i][d] = std::clamp(val, lb, ub);
                        }
                    }

                    // Откладка яиц
                    if (my_random_i(0, 2) == 1) {
                        for (size_t d = 0; d < dim; ++d) {
                            X_new_m[worst_m_idx][d] = lb + my_random_d(0, 1) * (ub - lb);
                            X_new_f[worst_f_idx][d] = lb + my_random_d(0, 1) * (ub - lb);
                        }
                    }
                }
            }
        }

        // Обновление мужской популяции
        for (size_t i = 0; i < Nm; ++i) {
            long double new_fit = objective_func(X_new_m[i]);
            total_evaluations++;

            if (new_fit < fitness_m[i]) {
                fitness_m[i] = new_fit;
                Xm[i] = X_new_m[i];

                if (new_fit < best_m_fit) {
                    best_m_fit = new_fit;
                    best_m_idx = i;
                    X_best_m = Xm[i];
                }
                if (new_fit > worst_m_fit) {
                    worst_m_fit = new_fit;
                    worst_m_idx = i;
                }
            }
        }

        // Обновление женской популяции
        for (size_t i = 0; i < Nf; ++i) {
            long double new_fit = objective_func(X_new_f[i]);
            total_evaluations++;

            if (new_fit < fitness_f[i]) {
                fitness_f[i] = new_fit;
                Xf[i] = X_new_f[i];

                if (new_fit < best_f_fit) {
                    best_f_fit = new_fit;
                    best_f_idx = i;
                    X_best_f = Xf[i];
                }
                if (new_fit > worst_f_fit) {
                    worst_f_fit = new_fit;
                    worst_f_idx = i;
                }
            }
        }

        // Обновление глобального лучшего
        if (best_m_fit < best_f_fit) {
            global_best_fit = best_m_fit;
            X_food = X_best_m;
        }
        else {
            global_best_fit = best_f_fit;
            X_food = X_best_f;
        }

        best_history[t] = global_best_fit;
        iteration_counts.push_back(total_evaluations);
    }

    return best_history;
}

void Middle2_parallel(
    size_t N,
    size_t T,
    long double (*objective_func)(Point),
    size_t dim,
    double lb,
    double ub,
    const std::string& plot_name,
    size_t num_trials = 10)
{
    std::vector<int> iteration_counts;
    std::vector<double> average_history(T, 0.0);
    std::vector<std::vector<long double>> Table(num_trials);
    std::vector<std::vector<int>> all_iterations(num_trials);

    std::cout << "Запуск " << num_trials << " параллельных испытаний SO..." << std::endl;

    // Счётчик завершённых задач
    volatile int completed = 0;

#pragma omp parallel for
    for (int trial = 0; trial < static_cast<int>(num_trials); ++trial) {
        std::vector<int> current_iterations;
        std::vector<long double> history = SnakeOptimization2(
            N, T, objective_func, dim, lb, ub, current_iterations);

#pragma omp critical
        {
            Table[trial] = history;
            all_iterations[trial] = current_iterations;
            completed++;

            // Простой прогресс-бар
            int percent = (completed * 100) / num_trials;
            std::cout << "\rПрогресс: " << percent << "% (" << completed
                << "/" << num_trials << ")" << std::flush;
        }
    }
    std::cout << std::endl << "Все испытания завершены!" << std::endl;

    // Сохраняем итерации первого запуска
    if (num_trials > 0) {
        iteration_counts = all_iterations[0];
    }

    // Суммирование и усреднение
    for (size_t trial = 0; trial < num_trials; ++trial) {
        for (size_t i = 0; i < T && i < Table[trial].size(); ++i) {
            average_history[i] += Table[trial][i];
        }
    }

    for (size_t i = 0; i < T; ++i) {
        average_history[i] /= num_trials;
    }

    // Сохранение
    std::string filename = "C:\\Users\\troit\\OneDrive\\Рабочий стол\\Science_result\\SO_parallel.csv";
    saveSetVectorsToCSV(Table, filename.c_str(), "SO");
    std::cout << "Результаты сохранены в " << filename << std::endl;

    // График
    plt::named_plot(plot_name, iteration_counts, average_history);
    plt::xlabel("number of calls Fitness function");
    plt::ylabel("Best_Fitness");
    plt::legend();
    plt::grid(true);
}

void Middle2(
    size_t N,
    size_t T,
    long double (*objective_func)(Point),
    size_t dim,
    double lb,
    double ub,

    const std::string& plot_name,
    size_t num_trials = 10)
{
    std::vector<int> iteration_counts;
    std::vector<double> average_history(T, 0.0);
    std::vector<std::vector<long double>> Table;
    for (size_t trial = 0; trial < num_trials; ++trial) {
        std::vector<int> current_iterations;
        std::vector<long double> history = SnakeOptimization2(
            N, T, objective_func, dim, lb, ub, current_iterations);
        Table.push_back(history);
        // Если это первый запуск, сохраняем счетчик итераций
        if (trial == 0) {
            iteration_counts = current_iterations;
        }

        // Суммируем историю
        for (size_t i = 0; i < T; ++i) {
            average_history[i] += history[i];
        }
    }

    // Усреднение
    for (size_t i = 0; i < T; ++i) {
        average_history[i] /= num_trials;
    }
    saveSetVectorsToCSV(Table, "C:\\Users\\troit\\OneDrive\\Рабочий стол\\Science_result\\SO.csv", "SO");
    // Построение графика
    plt::named_plot(plot_name, iteration_counts, average_history);
    plt::xlabel("number of calls Fitness function");
    plt::ylabel("Best_Fitness");
    //plt::title("Сходимость Snake Optimization");
    plt::legend();
    plt::grid(true);
}