#pragma once
#include <omp.h>
#include "Mutation.h"
#include "Function.h"
#include "LocalSearch.h"
//#include "DiffEval.h"
namespace plt = matplotlibcpp;
int ReferencePoint(double (*Phi)(size_t), std::vector<long double> Fitness, int k) {
	double r = my_random_d();
	double sum = 0;
	int m = 0;
	long double ref = Fitness[0];
	std::vector<long double> poweredF;
	for (int i = 0; i < Fitness.size(); ++i) {
		poweredF.push_back(pow(Fitness[i], Phi(k)));
		sum += poweredF[i];
	}
	double sum_p = 0;
	for (int i = 0; i < Fitness.size(); ++i) {
		sum_p += (poweredF[i]) / sum;
		if (sum_p > r) {
			m = i;
			break;
		}

	}
	return  m;
}

std::vector<long double> SoFA(int N, int T, long double (*J)(Point), double (*Phi)(size_t), double (*Sigma)(size_t), int dim, double a, double b, double Q, size_t M, std::vector<int>& X, long double (*Distr)(long double, double, double, double)) {
	double Variance = Sigma(0);
	//int m = 0;
	//double P;
	int min_index = 0;
	std::vector<double> ref_numb;
	int numb_iteration = 0;
	std::vector<double> FirstPoint;
	std::vector<double> SecondPoint;
	std::vector<long double> Fitness;
	std::vector<long double> exp_Fitness;

	//std::vector<std::vector<double>> Population=CreatePopulation(N,dim,a,b);
	Population Pop = Population(N, dim);

	numb_iteration += N;
	long double Fitness_min = std::numeric_limits<long double>::max();
	Pop.GeneratePopulation(a, b);
	for (int i = 0; i < N; ++i) {
		Fitness.push_back(J(Pop[i]));
		exp_Fitness.push_back(exponential_Fitness(Fitness[i], dim));
		if (Fitness[i] < Fitness_min) {
			Fitness_min = Fitness[i];
			min_index = i;
		}
	}
	std::vector<long double> Y;
	int Ref_index;
	//Ref_index = ReferencePoint(Phi, exp_Fitness, 0);
	long double Fitness_candidate;
	std::vector<int> Ref_list;
	long double noize;
	long double noize_middle;
	std::vector<int> noize_iter;
	std::vector<long double> noize_list;
	std::vector<long double> sigma_list;
	for (int i = 0; i < T; i++) {

		Ref_index = ReferencePoint(Phi, exp_Fitness, i);
		//std::cout << m << std::endl;
		ref_numb.push_back(Ref_index);
		//double q = my_random_d(0, 1);
		Point Further(dim);

		noize = 0.0;
		for (int j = 0; j < dim; ++j) {

			double q = my_random_d(0, 1);
			if (q < Q) {
				Variance = 2 * Sigma(i);
				Further[j] = Distr(Pop[Ref_index][j], a, b, Variance);
				if (!((i + 9) % 5000))
					noize += fabs(Further[j] - Pop[Ref_index][j]);
			}
			else {
				Variance = 2 * Sigma(i + M);
				Further[j] = Distr(Pop[Ref_index][j], a, b, Variance);
				if (!((i + 9) % 5000))
					noize += fabs(Further[j] - Pop[Ref_index][j]);
			}
		}

		if (!((i + 9) % 5000)) {
			sigma_list.push_back(Variance);
			noize_middle = noize / dim;
			noize_list.push_back(noize_middle);
			noize_iter.push_back(i);

		}
		Fitness_candidate = J(Further);
		numb_iteration++;
		if (Fitness_candidate < Fitness[Ref_index]) {
			Pop[Ref_index] = Further;
			Fitness[Ref_index] = Fitness_candidate;
			exp_Fitness[Ref_index] = exponential_Fitness(Fitness_candidate, dim);

			if (Fitness_candidate < Fitness_min) {
				Fitness_min = Fitness_candidate;
				min_index = Ref_index;
			}
		}
		X.push_back(numb_iteration);
		//Y.push_back(exponential_Fitness(Fitness_min,dim));
		Y.push_back(Fitness_min);
	}
	/*plt::plot(sigma_list, noize_list, ".");

		plt::xlabel("Sigma");
		plt::ylabel("Middle_Noize");
		plt::show();*/
		/*std::vector<double> colors;
		for (size_t j = 0; j < FirstPoint.size(); ++j) {
			colors.push_back(j); // или любое другое значение для цвета
		}
		const char* plot_name = "Go";
		std::vector<double> color
		plt::scatter(FirstPoint, SecondPoint,50); // 50 - размер точек

		plt::colorbar();
		plt::show();*/
		/*plt::plot(ref_numb, X, ".");
		plt::xlabel("Reference point index");
		plt::ylabel("Iteration");
		plt::show();
		plt::close();*/
	return Y;
}
std::vector<long double> SoFA_Time(int N, int T, long double (*J)(Point), double (*Phi)(size_t), double (*Sigma)(size_t), int dim, double a, double b, double Q, size_t M, std::vector<int>& X, long double (*Distr)(long double, double, double, double)) {
	double Variance = Sigma(0);
	//int m = 0;
	//double P;
	int min_index = 0;
	std::vector<double> ref_numb;
	//int numb_iteration = 0;
	std::vector<double> FirstPoint;
	std::vector<double> SecondPoint;
	std::vector<long double> Fitness;
	std::vector<long double> exp_Fitness;
	int start = clock();
	//std::vector<std::vector<double>> Population=CreatePopulation(N,dim,a,b);
	Population Pop = Population(N, dim);
	long double Fitness_min = std::numeric_limits<long double>::max();
	Pop.GeneratePopulation(a, b);
	for (int i = 0; i < N; ++i) {
		Fitness.push_back(J(Pop[i]));
		exp_Fitness.push_back(exponential_Fitness(Fitness[i], dim));
		if (Fitness[i] < Fitness_min) {
			Fitness_min = Fitness[i];
			min_index = i;
		}
	}
	std::vector<long double> Y;
	int Ref_index;
	//Ref_index = ReferencePoint(Phi, exp_Fitness, 0);
	long double Fitness_candidate;
	std::vector<int> Ref_list;
	long double noize;
	long double noize_middle;
	std::vector<int> noize_iter;
	std::vector<long double> noize_list;
	std::vector<long double> sigma_list;
	for (int i = 0; i < T; i++) {

		Ref_index = ReferencePoint(Phi, exp_Fitness, i);
		//std::cout << m << std::endl;
		ref_numb.push_back(Ref_index);
		//double q = my_random_d(0, 1);
		Point Further(dim);

		noize = 0.0;
		for (int j = 0; j < dim; ++j) {

			double q = my_random_d(0, 1);
			if (q < Q) {
				Variance = 2 * Sigma(i);
				Further[j] = Distr(Pop[Ref_index][j], a, b, Variance);
				if (!((i + 9) % 5000))
					noize += fabs(Further[j] - Pop[Ref_index][j]);
			}
			else {
				Variance = 2 * Sigma(i + M);
				Further[j] = Distr(Pop[Ref_index][j], a, b, Variance);
				if (!((i + 9) % 5000))
					noize += fabs(Further[j] - Pop[Ref_index][j]);
			}
		}

		if (!((i + 9) % 5000)) {
			sigma_list.push_back(Variance);
			noize_middle = noize / dim;
			noize_list.push_back(noize_middle);
			noize_iter.push_back(i);

		}
		Fitness_candidate = J(Further);
		if (Fitness_candidate < Fitness[Ref_index]) {
			Pop[Ref_index] = Further;
			Fitness[Ref_index] = Fitness_candidate;
			exp_Fitness[Ref_index] = exponential_Fitness(Fitness_candidate, dim);

			if (Fitness_candidate < Fitness_min) {
				Fitness_min = Fitness_candidate;
				min_index = Ref_index;
			}
		}
		int now = clock() - start;
		X.push_back(now);
		Y.push_back(Fitness_min);
	}

	return Y;
}

void MiddleSF(int N, int T, long double (*J)(Point), double (*Phi)(size_t), double (*Sigma)(size_t),
	int dim, double a, double b, double Q, int M, long double(*Distr)(long double, double, double, double),
	const char* plot_name, bool time_flag = false) {

	std::vector<long double> averaged(T, 0.0);
	std::vector<std::vector<long double>> Table;
	std::vector<int> dummy;
	std::vector<int> last_iterations;

	for (size_t run = 0; run < 10; ++run) {
		dummy.clear();

		std::vector<long double> history;
		if (!time_flag) {
			history = SoFA(N, T, J, Phi, Sigma, dim, a, b, Q, M, dummy, Distr);
		}
		else {
			history = SoFA_Time(N, T, J, Phi, Sigma, dim, a, b, Q, M, dummy, Distr);
		}

		// Сохраняем итерации последнего запуска для графика
		if (run == 9) {
			last_iterations = dummy;
		}

		// Суммируем историю для усреднения
		for (int i = 0; i < T; ++i) {
			averaged[i] += log(1 + history[i]);
		}

		Table.push_back(history);
	}

	// Усреднение
	for (int i = 0; i < T; ++i) {
		averaged[i] /= 10.0;
	}

	// Сохраняем результаты
	saveSetVectorsToCSV(Table, "C:\\Users\\troit\\OneDrive\\Рабочий стол\\Science_result\\SoFA.csv", "SoFA");

	// Строим график
	plt::named_plot(plot_name, last_iterations, averaged);
	plt::legend();
}

void MiddleSF_parallel(int N, int T, long double (*J)(Point), double (*Phi)(size_t), double (*Sigma)(size_t),
	int dim, double a, double b, double Q, int M, long double(*Distr)(long double, double, double, double),
	const char* plot_name, bool time_flag = false) {

	const int NUM_RUNS = 10;
	std::vector<std::vector<long double>> Table(NUM_RUNS);  // для хранения результатов
	std::vector<std::vector<int>> all_iterations(NUM_RUNS); // для хранения итераций

	// Параллельные запуски
#pragma omp parallel for
	for (int i = 0; i < NUM_RUNS; ++i) {
		std::vector<int> current_iter;
		//std::vector<double> curr_it;
		Point F;
		if (!time_flag) {
			F = SoFA(N, T, J, Phi, Sigma, dim, a, b, Q, M, current_iter, Distr);
		}
		else {
			F = SoFA_Time(N, T, J, Phi, Sigma, dim, a, b, Q, M, current_iter, Distr);
		}
		// Критическая секция для сохранения результатов
#pragma omp critical
		{
			Table[i] = F.getPoint();  // предполагаю, что Point можно преобразовать в vector
			all_iterations[i] = current_iter;
			std::cout << "Запуск " << i + 1 << " завершён" << std::endl;
		}
	}

	// Усреднение результатов (последовательно)
	std::vector<long double> averaged(T, 0.0);
	for (int run = 0; run < NUM_RUNS; ++run) {
		for (int t = 0; t < T; ++t) {
			averaged[t] += Table[run][t];
		}
	}
	for (int t = 0; t < T; ++t) {
		averaged[t] /= NUM_RUNS;
	}

	// Сохраняем результаты
	saveSetVectorsToCSV(Table, "C:\\Users\\troit\\OneDrive\\Рабочий стол\\Science_result\\SoFA_parallel.csv", "SoFA");
	std::cout << std::endl;
	std::cout << "SoFA " << exponential_Fitness(averaged[T - 1], dim);
	std::cout << std::endl;
	// Для графика используем итерации последнего запуска
	plt::named_plot(plot_name, all_iterations[NUM_RUNS - 1], averaged);
	plt::legend();
}

std::vector<long double> SoFAM(int N, int T, long double (*J)(Point), double (*Phi)(size_t), double (*Sigma)(size_t), double (*H)(size_t), int dim, double a, double b, double Q, size_t M, std::vector<int>& X, long double (*Distr)(long double, double, double, double)) {
	double Variance = Sigma(0);
	std::vector<double> ref_numb;
	int numb_iteration = 0;
	double Fitness_candidate;
	std::vector<long double> exp_Fitness;
	std::vector<long double> Fitness;
	Population Pop = Population(N, dim);
	long double Fitness_min = INT_MAX;
	int min_index = 0;
	numb_iteration += N;
	Pop.GeneratePopulation(a, b);
	for (int i = 0; i < N; ++i) {
		Fitness.push_back(J(Pop[i]));
		exp_Fitness.push_back(exponential_Fitness(Fitness[i], dim));
		if (Fitness[i] < Fitness_min) {
			Fitness_min = Fitness[i];
			min_index = i;
		}
	}
	std::vector<long double> Y;
	int Ref_index;
	for (size_t i = 0; i < T; i++) {

		Ref_index = ReferencePoint(Phi, exp_Fitness, i);
		//std::cout << m << std::endl;
		ref_numb.push_back(Ref_index);

		Point Further(dim);
		Point Further_without_mutation(dim);
		double q;
		//Variance = 2 * Sigma(i);
		std::vector<long double> New1 = Zero_Copy(exp_Fitness, Ref_index);
		int  Ref_index1 = ReferencePoint(Phi, New1, i);
		std::vector<long double> New2 = Zero_Copy(New1, Ref_index1);
		int  Ref_index2 = ReferencePoint(Phi, New2, i);
		for (int j = 0; j < dim; j++) {

			q = my_random_d();
			if (q < Q) {
				Variance = 2 * Sigma(i);
				Further_without_mutation[j] = std::min(b * 1.0L, (std::max(a * 1.0L, Pop[Ref_index][j] + (Pop[min_index][j] - Pop[Ref_index][j] + Pop[Ref_index1][j] - Pop[Ref_index2][j]) * H(i))));
				Further[j] = Distr(Further_without_mutation[j], a, b, Variance);
			}
			else {
				Variance = 2 * Sigma(i + M);
				Further[j] = Distr(Pop[Ref_index][j], a, b, Variance);
			}
		}


		Fitness_candidate = J(Further);
		numb_iteration++;
		if (Fitness_candidate < Fitness[Ref_index]) {
			Pop[Ref_index] = Further;
			Fitness[Ref_index] = Fitness_candidate;
			exp_Fitness[Ref_index] = exponential_Fitness(Fitness_candidate, dim);
			if (Fitness_candidate < Fitness_min) {
				Fitness_min = Fitness_candidate;
				min_index = Ref_index;
			}
		}

		X.push_back(numb_iteration);
		//Y.push_back(exponential_Fitness(Fitness_min,dim));
		Y.push_back(Fitness_min);
	}
	/*for (int i = 0; i < dim; ++i) {
		std::cout<<Pop[min_index][i]<<" ";
	}
	std::cout << "\n";*/
	/*plt::plot(ref_numb, X, ".");
	plt::xlabel("Reference point index");
	plt::ylabel("Iteration");
	plt::show();*/
	return Y;
}

void MiddleSM(int N, int T, long double (*J)(Point), double (*Phi)(size_t), double (*Sigma)(size_t), double (*H)(size_t),
	int dim, double a, double b, double Q, int M, long double(*Distr)(long double, double, double, double),
	const char* plot_name) {

	std::vector<int>Iteration;
	Point middle(T);
	Point F(T);
	std::vector<std::vector<long double>> Table;
	for (size_t i = 0; i < 10; ++i) {
		Iteration = {};
		F = SoFAM(N, T, J, Phi, Sigma, H, dim, a, b, Q, M, Iteration, Distr);
		std::vector<long double> F1 = F;
		for (int j = 0; j < F.getDimension(); ++j) {
			F[j] = log(F[j] + 1);
		}
		middle = middle + F;
		Table.push_back(F1);
	}

	saveSetVectorsToCSV(Table, "C:\\Users\\troit\\OneDrive\\Рабочий стол\\Science_result\\SO.csv", "SO");
	Point middle1 = middle * 0.1;
	plt::named_plot(plot_name, Iteration, middle1.getPoint());
	plt::legend();

}

void MiddleSM_parallel(int N, int T, long double (*J)(Point), double (*Phi)(size_t),
	double (*Sigma)(size_t), double (*H)(size_t),
	int dim, double a, double b, double Q, int M,
	long double(*Distr)(long double, double, double, double),
	const char* plot_name) {

	const int NUM_RUNS = 10;
	std::vector<Point> results(NUM_RUNS, Point(T));  // для хранения результатов каждого запуска
	std::vector<std::vector<int>> iterations(NUM_RUNS);  // для хранения итераций каждого запуска

	// Параллельные запуски
#pragma omp parallel for num_threads(10)
	for (int run = 0; run < NUM_RUNS; ++run) {
		std::vector<int> current_iter;
		Point F = SoFAM(N, T, J, Phi, Sigma, H, dim, a, b, Q, M, current_iter, Distr);

		// Критическая секция для сохранения результатов
#pragma omp critical
		{
			results[run] = F;
			iterations[run] = current_iter;
			std::cout << "Запуск " << run + 1 << " завершён" << std::endl;
		}
	}

	// Усреднение результатов (последовательно, это быстро)
	Point middle(T);
	for (int run = 0; run < NUM_RUNS; ++run) {
		middle = middle + results[run];
	}
	Point middle1 = middle * (1.0 / NUM_RUNS);  // 0.1 заменено на 1.0/10 для ясности

	// Для графика используем итерации последнего запуска
	// (можно также усреднить, но это сложнее)
	std::cout << std::endl;
	std::cout << "SoFAM " << exponential_Fitness(middle1[T - 1], dim);
	std::cout << std::endl;
	plt::named_plot(plot_name, iterations[NUM_RUNS - 1], middle1.getPoint());
	plt::legend();
}

std::vector<long double> SoFA_Var_Mod(int N, int T, long double (*J)(Point), double (*Phi)(size_t), double (*Sigma)(size_t), int dim, double a, double b, double Q, size_t M, std::vector<int>& X, long double (*Distr)(long double, double, double, double)) {
	double Variance = Sigma(0);
	//int m = 0;
	//double P;
	int min_index = 0;
	std::vector<double> ref_numb;
	int numb_iteration = 0;
	std::vector<double> FirstPoint;
	std::vector<double> SecondPoint;
	std::vector<long double> Fitness;
	std::vector<long double> exp_Fitness;

	//std::vector<std::vector<double>> Population=CreatePopulation(N,dim,a,b);
	Population Pop = Population(N, dim);

	numb_iteration += N;
	long double Fitness_min = std::numeric_limits<long double>::max();
	Pop.GeneratePopulation(a, b);
	for (int i = 0; i < N; ++i) {
		Fitness.push_back(J(Pop[i]));
		exp_Fitness.push_back(exponential_Fitness(Fitness[i], dim));
		if (Fitness[i] < Fitness_min) {
			Fitness_min = Fitness[i];
			min_index = i;
		}
	}
	std::vector<long double> Y;
	int Ref_index;
	//Ref_index = ReferencePoint(Phi, exp_Fitness, 0);
	long double Fitness_candidate;
	std::vector<int> Ref_list;
	long double noize;
	long double noize_middle;
	std::vector<int> noize_iter;
	std::vector<long double> noize_list;
	std::vector<long double> sigma_list;
	for (int i = 0; i < T; i++) {
		Ref_index = ReferencePoint(Phi, exp_Fitness, i);
		ref_numb.push_back(Ref_index);
		//double q = my_random_d(0, 1);
		Point Further(dim);
		Point Ref_Grad = Gradient(Pop[Ref_index], J);
		Ref_Grad.Normalization();
		numb_iteration += dim;
		noize = 0.0;
		for (int j = 0; j < dim; ++j) {
			double q = my_random_d(0, 1);
			if (q < Q) {
				Variance = 2 * Sigma(i) * Ref_Grad[j];
				Further[j] = Distr(Pop[Ref_index][j], a, b, Variance);
				//if (!((i + 9) % 5000))
					//noize += fabs(Further[j] - Pop[Ref_index][j]);
			}
			else {
				Variance = 2 * Sigma(i + M) * Ref_Grad[j];
				Further[j] = Distr(Pop[Ref_index][j], a, b, Variance);
				//if (!((i + 9) % 5000))
					//noize += fabs(Further[j] - Pop[Ref_index][j]);
			}
		}

		/*if (!((i + 9) % 5000)) {
			sigma_list.push_back(Variance);
			noize_middle = noize / dim;
			noize_list.push_back(noize_middle);
			noize_iter.push_back(i);

		}*/
		Fitness_candidate = J(Further);
		numb_iteration++;
		if (Fitness_candidate < Fitness[Ref_index]) {
			Pop[Ref_index] = Further;
			Fitness[Ref_index] = Fitness_candidate;
			exp_Fitness[Ref_index] = exponential_Fitness(Fitness_candidate, dim);

			if (Fitness_candidate < Fitness_min) {
				Fitness_min = Fitness_candidate;
				min_index = Ref_index;
			}
		}
		X.push_back(numb_iteration);
		//Y.push_back(exponential_Fitness(Fitness_min,dim));
		Y.push_back(Fitness_min);
	}
	/*plt::plot(sigma_list, noize_list, ".");

		plt::xlabel("Sigma");
		plt::ylabel("Middle_Noize");
		plt::show();*/
		//std::vector<double> colors;
	return Y;
}

std::vector<long double> SoFA_Var_Mod_Time(int N, int T, long double (*J)(Point), double (*Phi)(size_t), double (*Sigma)(size_t), int dim, double a, double b, double Q, size_t M, std::vector<int>& X, long double (*Distr)(long double, double, double, double)) {
	double Variance = Sigma(0);
	//int m = 0;
	//double P;
	int min_index = 0;
	std::vector<double> ref_numb;
	int numb_iteration = 0;
	std::vector<double> FirstPoint;
	std::vector<double> SecondPoint;
	std::vector<long double> Fitness;
	std::vector<long double> exp_Fitness;

	//std::vector<std::vector<double>> Population=CreatePopulation(N,dim,a,b);
	Population Pop = Population(N, dim);

	numb_iteration += N;
	long double Fitness_min = std::numeric_limits<long double>::max();
	Pop.GeneratePopulation(a, b);
	for (int i = 0; i < N; ++i) {
		Fitness.push_back(J(Pop[i]));
		exp_Fitness.push_back(exponential_Fitness(Fitness[i], dim));
		if (Fitness[i] < Fitness_min) {
			Fitness_min = Fitness[i];
			min_index = i;
		}
	}
	std::vector<long double> Y;
	int Ref_index;
	//Ref_index = ReferencePoint(Phi, exp_Fitness, 0);
	long double Fitness_candidate;
	std::vector<int> Ref_list;
	long double noize;
	long double noize_middle;
	std::vector<int> noize_iter;
	std::vector<long double> noize_list;
	std::vector<long double> sigma_list;
	int start = clock();
	for (int i = 0; i < T; i++) {
		Ref_index = ReferencePoint(Phi, exp_Fitness, i);
		//std::cout << m << std::endl;
		ref_numb.push_back(Ref_index);
		//double q = my_random_d(0, 1);
		Point Further(dim);
		Point Ref_Grad = Gradient(Pop[Ref_index], J);
		Ref_Grad.Normalization();
		numb_iteration += dim;
		noize = 0.0;
		for (int j = 0; j < dim; ++j) {

			double q = my_random_d(0, 1);
			if (q < Q) {
				Variance = 2 * Sigma(i) * Ref_Grad[j];
				Further[j] = Distr(Pop[Ref_index][j], a, b, Variance);
				//if (!((i + 9) % 5000))
					//noize += fabs(Further[j] - Pop[Ref_index][j]);
			}
			else {
				Variance = 2 * Sigma(i + M) * Ref_Grad[j];
				Further[j] = Distr(Pop[Ref_index][j], a, b, Variance);
				//if (!((i + 9) % 5000))
					//noize += fabs(Further[j] - Pop[Ref_index][j]);
			}
		}

		/*if (!((i + 9) % 5000)) {
			sigma_list.push_back(Variance);
			noize_middle = noize / dim;
			noize_list.push_back(noize_middle);
			noize_iter.push_back(i);

		}*/
		Fitness_candidate = J(Further);
		numb_iteration++;
		if (Fitness_candidate < Fitness[Ref_index]) {
			Pop[Ref_index] = Further;
			Fitness[Ref_index] = Fitness_candidate;
			exp_Fitness[Ref_index] = exponential_Fitness(Fitness_candidate, dim);

			if (Fitness_candidate < Fitness_min) {
				Fitness_min = Fitness_candidate;
				min_index = Ref_index;
			}
		}
		int now = clock() - start;
		X.push_back(now);
		//Y.push_back(exponential_Fitness(Fitness_min,dim));
		Y.push_back(Fitness_min);
	}
	/*plt::plot(sigma_list, noize_list, ".");

		plt::xlabel("Sigma");
		plt::ylabel("Middle_Noize");
		plt::show();*/
		//std::vector<double> colors;
	return Y;
}

std::vector<long double> SoFAM_Grad(int N, int T, long double (*J)(Point), double (*Phi)(size_t), double (*Sigma)(size_t), int dim, double a, double b, double Q, size_t M, std::vector<int>& X, long double (*Distr)(long double, double, double, double)) {
	std::vector<double> Variance(dim, Sigma(0));
	//int m = 0;
	//double P;
	int min_index = 0;
	std::vector<double> ref_numb;
	int numb_iteration = 0;
	std::vector<double> FirstPoint;
	std::vector<double> SecondPoint;
	std::vector<long double> Fitness;
	std::vector<long double> exp_Fitness;

	//std::vector<std::vector<double>> Population=CreatePopulation(N,dim,a,b);
	Population Pop = Population(N, dim);

	numb_iteration += N;
	long double Fitness_min = std::numeric_limits<long double>::max();
	Pop.GeneratePopulation(a, b);
	for (int i = 0; i < N; ++i) {
		Fitness.push_back(J(Pop[i]));
		exp_Fitness.push_back(exponential_Fitness(Fitness[i], dim));
		if (Fitness[i] < Fitness_min) {
			Fitness_min = Fitness[i];
			min_index = i;
		}
	}
	numb_iteration += N;
	std::vector<long double> Y;
	int Ref_index;
	//Ref_index = ReferencePoint(Phi, exp_Fitness, 0);
	long double Fitness_candidate;
	std::vector<int> Ref_list;
	long double noize;
	long double noize_middle;
	std::vector<int> noize_iter;
	std::vector<long double> noize_list;
	std::vector<long double> sigma_list;
	for (int i = 0; i < T; i++) {
		Ref_index = ReferencePoint(Phi, exp_Fitness, i);
		//std::cout << m << std::endl;
		ref_numb.push_back(Ref_index);
		//double q = my_random_d(0, 1);
		Point Further(dim);
		Point Further_Swarm(dim);
		Point Ref_Grad = Gradient(Pop[Ref_index], J);
		Ref_Grad.Normalization();
		numb_iteration += dim;

		noize = 0.0;
		double k = 0.1;
		while (true) {
			for (int j = 0; j < dim; ++j) {

				double q = my_random_d(0, 1);

				if (q < Q) {

					Variance[j] = 2 * Sigma(i);
					Further_Swarm[j] = Pop[Ref_index][j] - k * Ref_Grad[j];
					//Further[j] = Distr(Further_Swarm[j], a, b, Variance);
					//if (!((i + 9) % 5000))
						//noize += fabs(Further[j] - Pop[Ref_index][j]);
				}
				else {
					Further_Swarm[j] = Pop[Ref_index][j];
					Variance[j] = 2 * Sigma(i + M);
					//Further[j] = Distr(Further_Swarm[j], a, b, Variance);
					//if (!((i + 9) % 5000))
						//noize += fabs(Further[j] - Pop[Ref_index][j]);
				}
			}
			if (J(Further_Swarm) > Fitness[Ref_index]) {
				k *= 0.5;
				numb_iteration++;
				continue;

			}
			else {
				break;
			}
		}
		for (int j = 0; j < dim; ++j) {
			Further[j] = Distr(Further_Swarm[j], a, b, Variance[j]);
		}
		/*if (!((i + 9) % 5000)) {
			sigma_list.push_back(Variance);
			noize_middle = noize / dim;
			noize_list.push_back(noize_middle);
			noize_iter.push_back(i);

		}*/
		Fitness_candidate = J(Further);
		numb_iteration++;
		if (Fitness_candidate < Fitness[Ref_index]) {
			Pop[Ref_index] = Further;
			Fitness[Ref_index] = Fitness_candidate;
			exp_Fitness[Ref_index] = exponential_Fitness(Fitness_candidate, dim);

			if (Fitness_candidate < Fitness_min) {
				Fitness_min = Fitness_candidate;
				min_index = Ref_index;
			}
		}
		X.push_back(numb_iteration);
		//Y.push_back(exponential_Fitness(Fitness_min,dim));
		Y.push_back(Fitness_min);
	}
	/*plt::plot(sigma_list, noize_list, ".");

		plt::xlabel("Sigma");
		plt::ylabel("Middle_Noize");
		plt::show();*/
		//std::vector<double> colors;
	return Y;
}

std::vector<long double> SoFAM_Grad_Time(int N, int T, long double (*J)(Point), double (*Phi)(size_t), double (*Sigma)(size_t), int dim, double a, double b, double Q, size_t M, std::vector<int>& X, long double (*Distr)(long double, double, double, double)) {
	double Variance = Sigma(0);
	//int m = 0;
	//double P;
	int min_index = 0;
	std::vector<double> ref_numb;
	int numb_iteration = 0;
	std::vector<double> FirstPoint;
	std::vector<double> SecondPoint;
	std::vector<long double> Fitness;
	std::vector<long double> exp_Fitness;

	//std::vector<std::vector<double>> Population=CreatePopulation(N,dim,a,b);
	Population Pop = Population(N, dim);

	numb_iteration += N;
	long double Fitness_min = std::numeric_limits<long double>::max();
	Pop.GeneratePopulation(a, b);
	for (int i = 0; i < N; ++i) {
		Fitness.push_back(J(Pop[i]));
		exp_Fitness.push_back(exponential_Fitness(Fitness[i], dim));
		if (Fitness[i] < Fitness_min) {
			Fitness_min = Fitness[i];
			min_index = i;
		}
	}
	std::vector<long double> Y;
	int Ref_index;
	//Ref_index = ReferencePoint(Phi, exp_Fitness, 0);
	long double Fitness_candidate;
	std::vector<int> Ref_list;
	long double noize;
	long double noize_middle;
	std::vector<int> noize_iter;
	std::vector<long double> noize_list;
	std::vector<long double> sigma_list;
	int start = clock();
	for (int i = 0; i < T; i++) {
		Ref_index = ReferencePoint(Phi, exp_Fitness, i);
		//std::cout << m << std::endl;
		ref_numb.push_back(Ref_index);
		//double q = my_random_d(0, 1);
		Point Further(dim);
		Point Ref_Grad = Gradient(Pop[Ref_index], J);
		Ref_Grad.Normalization();
		numb_iteration += dim;
		noize = 0.0;
		for (int j = 0; j < dim; ++j) {

			double q = my_random_d(0, 1);
			if (q < Q) {
				Variance = 2 * Sigma(i) * Ref_Grad[j];
				Further[j] = Distr(Pop[Ref_index][j], a, b, Variance);
				//if (!((i + 9) % 5000))
					//noize += fabs(Further[j] - Pop[Ref_index][j]);
			}
			else {
				Variance = 2 * Sigma(i + M);
				Further[j] = Distr(Pop[Ref_index][j], a, b, Variance);
				//if (!((i + 9) % 5000))
					//noize += fabs(Further[j] - Pop[Ref_index][j]);
			}
		}

		/*if (!((i + 9) % 5000)) {
			sigma_list.push_back(Variance);
			noize_middle = noize / dim;
			noize_list.push_back(noize_middle);
			noize_iter.push_back(i);

		}*/
		Fitness_candidate = J(Further);
		numb_iteration++;
		if (Fitness_candidate < Fitness[Ref_index]) {
			Pop[Ref_index] = Further;
			Fitness[Ref_index] = Fitness_candidate;
			exp_Fitness[Ref_index] = exponential_Fitness(Fitness_candidate, dim);

			if (Fitness_candidate < Fitness_min) {
				Fitness_min = Fitness_candidate;
				min_index = Ref_index;
			}
		}
		int now = clock() - start;
		X.push_back(now);
		//Y.push_back(exponential_Fitness(Fitness_min,dim));
		Y.push_back(Fitness_min);
	}
	/*plt::plot(sigma_list, noize_list, ".");

		plt::xlabel("Sigma");
		plt::ylabel("Middle_Noize");
		plt::show();*/
		//std::vector<double> colors;
	return Y;
}
void MiddleSM_Grad(int N, int T, long double (*J)(Point), double (*Phi)(size_t), double (*Sigma)(size_t),
	int dim, double a, double b, double Q, int M, long double(*Distr)(long double, double, double, double),
	const char* plot_name, bool time_flag = false) {

	std::vector<int>Iteration;
	Point middle(T);
	Point F(T);
	std::vector<std::vector<long double>> Table;
	for (size_t i = 0; i < 10; ++i) {
		Iteration = {};
		if (!time_flag) {
			F = SoFAM_Grad(N, T, J, Phi, Sigma, dim, a, b, Q, M, Iteration, Distr);
		}
		else {
			F = SoFAM_Grad_Time(N, T, J, Phi, Sigma, dim, a, b, Q, M, Iteration, Distr);
		}
		for (int j = 0; j < F.getDimension(); ++j) {
			F[j] = log(F[j] + 1);
		}
		middle = middle + F;
		Table.push_back(F);
	}
	saveSetVectorsToCSV(Table, "C:\\Users\\troit\\OneDrive\\Рабочий стол\\Science_result\\SoFA_Var_Mod.csv", "SoFA_Var_Mod");
	Point middle1 = middle * 0.1;
	plt::named_plot(plot_name, Iteration, middle1.getPoint());
	plt::legend();
}


/*std::vector<long double> SoFA_Var_Mod_Time(int N, int T, long double (*J)(Point), double (*Phi)(size_t), double (*Sigma)(size_t), int dim, double a, double b, double Q, size_t M, std::vector<int>& X, long double (*Distr)(long double, double, double, double)) {
	double Variance = Sigma(0);
	//int m = 0;
	//double P;
	int min_index = 0;
	std::vector<double> ref_numb;
	int numb_iteration = 0;
	std::vector<double> FirstPoint;
	std::vector<double> SecondPoint;
	std::vector<long double> Fitness;
	std::vector<long double> exp_Fitness;

	//std::vector<std::vector<double>> Population=CreatePopulation(N,dim,a,b);
	Population Pop = Population(N, dim);

	numb_iteration += N;
	long double Fitness_min = std::numeric_limits<long double>::max();
	Pop.GeneratePopulation(a, b);
	for (int i = 0; i < N; ++i) {
		Fitness.push_back(J(Pop[i]));
		exp_Fitness.push_back(exponential_Fitness(Fitness[i], dim));
		if (Fitness[i] < Fitness_min) {
			Fitness_min = Fitness[i];
			min_index = i;
		}
	}
	std::vector<long double> Y;
	int Ref_index;
	//Ref_index = ReferencePoint(Phi, exp_Fitness, 0);
	long double Fitness_candidate;
	std::vector<int> Ref_list;
	long double noize;
	long double noize_middle;
	std::vector<int> noize_iter;
	std::vector<long double> noize_list;
	std::vector<long double> sigma_list;
	int start = clock();
	for (int i = 0; i < T; i++) {
		Ref_index = ReferencePoint(Phi, exp_Fitness, i);
		//std::cout << m << std::endl;
		ref_numb.push_back(Ref_index);
		//double q = my_random_d(0, 1);
		Point Further(dim);
		Point Ref_Grad = Gradient(Pop[Ref_index], J);
		Ref_Grad.Normalization();
		numb_iteration += dim;
		noize = 0.0;
		for (int j = 0; j < dim; ++j) {

			double q = my_random_d(0, 1);
			if (q < Q) {
				Variance = 2 * Sigma(i) * Ref_Grad[j];
				Further[j] = Distr(Pop[Ref_index][j], a, b, Variance);
				//if (!((i + 9) % 5000))
					//noize += fabs(Further[j] - Pop[Ref_index][j]);
			}
			else {
				Variance = 2 * Sigma(i + M);
				Further[j] = Distr(Pop[Ref_index][j], a, b, Variance);
				//if (!((i + 9) % 5000))
					//noize += fabs(Further[j] - Pop[Ref_index][j]);
			}
		}

		if (!((i + 9) % 5000)) {
			sigma_list.push_back(Variance);
			noize_middle = noize / dim;
			noize_list.push_back(noize_middle);
			noize_iter.push_back(i);

		}
		Fitness_candidate = J(Further);
		numb_iteration++;
		if (Fitness_candidate < Fitness[Ref_index]) {
			Pop[Ref_index] = Further;
			Fitness[Ref_index] = Fitness_candidate;
			exp_Fitness[Ref_index] = exponential_Fitness(Fitness_candidate, dim);

			if (Fitness_candidate < Fitness_min) {
				Fitness_min = Fitness_candidate;
				min_index = Ref_index;
			}
		}
		int now = clock() - start;
		X.push_back(now);
		//Y.push_back(exponential_Fitness(Fitness_min,dim));
		Y.push_back(Fitness_min);
	}
	plt::plot(sigma_list, noize_list, ".");

		plt::xlabel("Sigma");
		plt::ylabel("Middle_Noize");
		plt::show();
		//std::vector<double> colors;
	return Y;
}*/

void MiddleSFVar_parallel(int N, int T, long double (*J)(Point), double (*Phi)(size_t), double (*Sigma)(size_t),
	int dim, double a, double b, double Q, int M, long double(*Distr)(long double, double, double, double),
	const char* plot_name, bool time_flag = false) {

	const int NUM_RUNS = 10;
	std::vector<std::vector<long double>> Table(NUM_RUNS);  // для хранения результатов
	std::vector<std::vector<int>> all_iterations(NUM_RUNS); // для хранения итераций

	// Параллельные запуски
#pragma omp parallel for
	for (int i = 0; i < NUM_RUNS; ++i) {
		std::vector<int> current_iter;
		Point F;
		if (time_flag) {
			F = SoFA_Var_Mod_Time(N, T, J, Phi, Sigma, dim, a, b, Q, M, current_iter, Distr);
		}
		else {
			F = SoFA_Var_Mod(N, T, J, Phi, Sigma, dim, a, b, Q, M, current_iter, Distr);
		}

		// Критическая секция для сохранения результатов
#pragma omp critical
		{
			Table[i] = F.getPoint();  // предполагаю, что Point можно преобразовать в vector
			all_iterations[i] = current_iter;
			std::cout << "Запуск " << i + 1 << " завершён" << std::endl;
		}
	}

	// Усреднение результатов (последовательно)
	std::vector<long double> averaged(T, 0.0);
	for (int run = 0; run < NUM_RUNS; ++run) {
		for (int t = 0; t < T; ++t) {
			averaged[t] += Table[run][t];
		}
	}
	for (int t = 0; t < T; ++t) {
		averaged[t] /= NUM_RUNS;
	}

	// Сохраняем результаты
	saveSetVectorsToCSV(Table, "C:\\Users\\troit\\OneDrive\\Рабочий стол\\Science_result\\SoFA_parallel.csv", "SoFA");
	std::cout << std::endl;
	std::cout << "SoFA_Var_Mod " << exponential_Fitness(averaged[T - 1], dim);
	std::cout << std::endl;
	// Для графика используем итерации последнего запуска
	plt::named_plot(plot_name, all_iterations[NUM_RUNS - 1], averaged);
	plt::legend();
}

void MiddleSF_Var(int N, int T, long double (*J)(Point), double (*Phi)(size_t), double (*Sigma)(size_t),
	int dim, double a, double b, double Q, int M, long double(*Distr)(long double, double, double, double),
	const char* plot_name, bool time_flag = false) {

	std::vector<int>Iteration;
	Point middle(T);
	Point F(T);
	std::vector<std::vector<long double>> Table;
	for (size_t i = 0; i < 10; ++i) {
		Iteration = {};
		if (time_flag) {
			F = SoFA_Var_Mod_Time(N, T, J, Phi, Sigma, dim, a, b, Q, M, Iteration, Distr);
		}
		else {
			F = SoFA_Var_Mod(N, T, J, Phi, Sigma, dim, a, b, Q, M, Iteration, Distr);
		}
		for (int j = 0; j < F.getDimension(); ++j) {
			F[j] = log(F[j] + 1);
		}
		middle = middle + F;
		Table.push_back(F);
	}
	saveSetVectorsToCSV(Table, "C:\\Users\\troit\\OneDrive\\Рабочий стол\\Science_result\\SoFA_Var_Mod.csv", "SoFA_Var_Mod");
	Point middle1 = middle * 0.1;
	plt::named_plot(plot_name, Iteration, middle1.getPoint());
	plt::legend();
}

std::vector<long double> SoFA_Var_Mod_Numb(int N, int T, long double (*J)(Point), double (*Phi)(size_t), double (*Sigma)(size_t), int dim, double a, double b, double Q, size_t M, std::vector<int>& X, long double (*Distr)(long double, double, double, double)) {
	double Variance = Sigma(0);
	//int m = 0;
	//double P;
	int min_index = 0;
	std::vector<double> ref_numb;
	int numb_iteration = 0;
	std::vector<double> FirstPoint;
	std::vector<double> SecondPoint;
	std::vector<long double> Fitness;
	std::vector<long double> exp_Fitness;

	//std::vector<std::vector<double>> Population=CreatePopulation(N,dim,a,b);
	Population Pop = Population(N, dim);

	numb_iteration += N;
	long double Fitness_min = std::numeric_limits<long double>::max();
	Pop.GeneratePopulation(a, b);
	std::vector<std::vector<DistIndex>> DI = CreatePopulationSet(Pop);
	for (int i = 0; i < N; ++i) {
		Fitness.push_back(J(Pop[i]));
		exp_Fitness.push_back(exponential_Fitness(Fitness[i], dim));
		if (Fitness[i] < Fitness_min) {
			Fitness_min = Fitness[i];
			min_index = i;
		}
	}
	std::vector<long double> Y;
	int Ref_index;
	//Ref_index = ReferencePoint(Phi, exp_Fitness, 0);
	long double Fitness_candidate;
	std::vector<int> Ref_list;
	long double noize;
	long double noize_middle;
	std::vector<int> noize_iter;
	std::vector<long double> noize_list;
	std::vector<long double> sigma_list;
	for (int i = 0; i < T; i++) {
		Ref_index = ReferencePoint(Phi, exp_Fitness, i);
		//std::cout << m << std::endl;
		ref_numb.push_back(Ref_index);
		//double q = my_random_d(0, 1);
		Point Further(dim);
		std::vector<DistIndex> D = DI[Ref_index];
		std::vector<DistIndex> DIK = getKNearest(D, dim);
		std::vector <Point> MatrToSolv = MatrixToSolve(Pop[Ref_index], D, Pop);
		Point SystPoint;
		std::vector<long double> B;
		for (int k = 0; k < DIK.size(); ++k) {
			SystPoint = Pop[Ref_index] - Pop[DIK[k].getIndex()];
			long double Norm = SystPoint.norm();
			B.push_back(abs(Fitness[Ref_index] - Fitness[DIK[k].getIndex()]) / Norm);
		}
		Point Ref_Grad = solveLinearSystem(MatrToSolv, B);
		Ref_Grad.Normalization();
		noize = 0.0;
		for (int j = 0; j < dim; ++j) {

			double q = my_random_d(0, 1);
			if (q < Q) {
				Variance = 2 * Sigma(i) * Ref_Grad[j];
				Further[j] = Distr(Pop[Ref_index][j], a, b, Variance);
				//if (!((i + 9) % 5000))
					//noize += fabs(Further[j] - Pop[Ref_index][j]);
			}
			else {
				Variance = 2 * Sigma(i + M) * Ref_Grad[j];
				Further[j] = Distr(Pop[Ref_index][j], a, b, Variance);
				//if (!((i + 9) % 5000))
					//noize += fabs(Further[j] - Pop[Ref_index][j]);
			}
		}

		/*if (!((i + 9) % 5000)) {
			sigma_list.push_back(Variance);
			noize_middle = noize / dim;
			noize_list.push_back(noize_middle);
			noize_iter.push_back(i);

		}*/
		Fitness_candidate = J(Further);
		numb_iteration++;
		if (Fitness_candidate < Fitness[Ref_index]) {
			Pop[Ref_index] = Further;
			Fitness[Ref_index] = Fitness_candidate;
			exp_Fitness[Ref_index] = exponential_Fitness(Fitness_candidate, dim);

			if (Fitness_candidate < Fitness_min) {
				Fitness_min = Fitness_candidate;
				min_index = Ref_index;
			}
		}
		X.push_back(numb_iteration);
		//Y.push_back(exponential_Fitness(Fitness_min,dim));
		Y.push_back(Fitness_min);
	}
	/*plt::plot(sigma_list, noize_list, ".");

		plt::xlabel("Sigma");
		plt::ylabel("Middle_Noize");
		plt::show();*/
		//std::vector<double> colors;
	return Y;
}

std::vector<long double> SoFA_Var_Mod_Numb_Time(int N, int T, long double (*J)(Point), double (*Phi)(size_t), double (*Sigma)(size_t), int dim, double a, double b, double Q, size_t M, std::vector<int>& X, long double (*Distr)(long double, double, double, double)) {
	double Variance = Sigma(0);
	//int m = 0;
	//double P;
	int min_index = 0;
	std::vector<double> ref_numb;
	int numb_iteration = 0;
	std::vector<double> FirstPoint;
	std::vector<double> SecondPoint;
	std::vector<long double> Fitness;
	std::vector<long double> exp_Fitness;

	//std::vector<std::vector<double>> Population=CreatePopulation(N,dim,a,b);
	Population Pop = Population(N, dim);

	numb_iteration += N;
	long double Fitness_min = std::numeric_limits<long double>::max();
	Pop.GeneratePopulation(a, b);
	std::vector<std::vector<DistIndex>> DI = CreatePopulationSet(Pop);
	for (int i = 0; i < N; ++i) {
		Fitness.push_back(J(Pop[i]));
		exp_Fitness.push_back(exponential_Fitness(Fitness[i], dim));
		if (Fitness[i] < Fitness_min) {
			Fitness_min = Fitness[i];
			min_index = i;
		}
	}
	std::vector<long double> Y;
	int Ref_index;
	//Ref_index = ReferencePoint(Phi, exp_Fitness, 0);
	long double Fitness_candidate;
	std::vector<int> Ref_list;
	long double noize;
	long double noize_middle;
	std::vector<int> noize_iter;
	std::vector<long double> noize_list;
	std::vector<long double> sigma_list;
	int start = clock();
	for (int i = 0; i < T; i++) {
		Ref_index = ReferencePoint(Phi, exp_Fitness, i);
		//std::cout << m << std::endl;
		ref_numb.push_back(Ref_index);
		//double q = my_random_d(0, 1);
		Point Further(dim);
		std::vector<DistIndex> D = DI[Ref_index];
		std::vector<DistIndex> DIK = getKNearest(D, dim);
		std::vector <Point> MatrToSolv = MatrixToSolve(Pop[Ref_index], D, Pop);
		Point SystPoint;
		std::vector<long double> B;
		for (int k = 0; k < DIK.size(); ++k) {
			SystPoint = Pop[Ref_index] - Pop[DIK[k].getIndex()];
			long double Norm = SystPoint.norm();
			B.push_back(abs(Fitness[Ref_index] - Fitness[DIK[k].getIndex()]) / Norm);
		}
		Point Ref_Grad = solveLinearSystem(MatrToSolv, B);
		Ref_Grad.Normalization();
		noize = 0.0;
		for (int j = 0; j < dim; ++j) {

			double q = my_random_d(0, 1);
			if (q < Q) {
				Variance = 2 * Sigma(i) * Ref_Grad[j];
				Further[j] = Distr(Pop[Ref_index][j], a, b, Variance);
				//if (!((i + 9) % 5000))
					//noize += fabs(Further[j] - Pop[Ref_index][j]);
			}
			else {
				Variance = 2 * Sigma(i + M) * Ref_Grad[j];
				Further[j] = Distr(Pop[Ref_index][j], a, b, Variance);
				//if (!((i + 9) % 5000))
					//noize += fabs(Further[j] - Pop[Ref_index][j]);
			}
		}

		/*if (!((i + 9) % 5000)) {
			sigma_list.push_back(Variance);
			noize_middle = noize / dim;
			noize_list.push_back(noize_middle);
			noize_iter.push_back(i);

		}*/
		Fitness_candidate = J(Further);
		numb_iteration++;
		if (Fitness_candidate < Fitness[Ref_index]) {
			Pop[Ref_index] = Further;
			Fitness[Ref_index] = Fitness_candidate;
			exp_Fitness[Ref_index] = exponential_Fitness(Fitness_candidate, dim);

			if (Fitness_candidate < Fitness_min) {
				Fitness_min = Fitness_candidate;
				min_index = Ref_index;
			}
		}
		int now = clock() - start;
		X.push_back(now);
		//Y.push_back(exponential_Fitness(Fitness_min,dim));
		Y.push_back(Fitness_min);
	}
	/*plt::plot(sigma_list, noize_list, ".");

		plt::xlabel("Sigma");
		plt::ylabel("Middle_Noize");
		plt::show();*/
		//std::vector<double> colors;
	return Y;
}
void MiddleSF_Var_Numb_parallel(int N, int T, long double (*J)(Point), double (*Phi)(size_t), double (*Sigma)(size_t),
	int dim, double a, double b, double Q, int M, long double(*Distr)(long double, double, double, double),
	const char* plot_name, bool time_flag = false) {

	const int NUM_RUNS = 10;
	std::vector<std::vector<long double>> Table(NUM_RUNS);  // для хранения результатов
	std::vector<std::vector<int>> all_iterations(NUM_RUNS); // для хранения итераций

	// Параллельные запуски
#pragma omp parallel for
	for (int i = 0; i < NUM_RUNS; ++i) {
		std::vector<int> current_iter;
		Point F;
		if (time_flag) {
			F = SoFA_Var_Mod_Numb_Time(N, T, J, Phi, Sigma, dim, a, b, Q, M, current_iter, Distr);
		}
		else {
			F = SoFA_Var_Mod_Numb(N, T, J, Phi, Sigma, dim, a, b, Q, M, current_iter, Distr);
		}

		// Критическая секция для сохранения результатов
#pragma omp critical
		{
			Table[i] = F.getPoint();  // предполагаю, что Point можно преобразовать в vector
			all_iterations[i] = current_iter;
			std::cout << "Запуск " << i + 1 << " завершён" << std::endl;
		}
	}

	// Усреднение результатов (последовательно)
	std::vector<long double> averaged(T, 0.0);
	for (int run = 0; run < NUM_RUNS; ++run) {
		for (int t = 0; t < T; ++t) {
			averaged[t] += Table[run][t];
		}
	}
	for (int t = 0; t < T; ++t) {
		averaged[t] /= NUM_RUNS;
	}

	// Сохраняем результаты
	saveSetVectorsToCSV(Table, "C:\\Users\\troit\\OneDrive\\Рабочий стол\\Science_result\\SoFA_parallel.csv", "SoFA");
	std::cout << std::endl;
	std::cout << "SoFA_Var_Mod " << exponential_Fitness(averaged[T - 1], dim);
	std::cout << std::endl;
	// Для графика используем итерации последнего запуска
	plt::named_plot(plot_name, all_iterations[NUM_RUNS - 1], averaged);
	plt::legend();
}
void MiddleSF_Var_Numb(int N, int T, long double (*J)(Point), double (*Phi)(size_t), double (*Sigma)(size_t),
	int dim, double a, double b, double Q, int M, long double(*Distr)(long double, double, double, double),
	const char* plot_name, bool time_flag = false) {

	std::vector<int>Iteration;
	Point middle(T);
	Point F(T);
	std::vector<std::vector<long double>> Table;
	for (size_t i = 0; i < 10; ++i) {
		Iteration = {};
		if (time_flag) {
			F = SoFA_Var_Mod_Numb_Time(N, T, J, Phi, Sigma, dim, a, b, Q, M, Iteration, Distr);
		}
		else {
			F = SoFA_Var_Mod_Numb(N, T, J, Phi, Sigma, dim, a, b, Q, M, Iteration, Distr);
		}
		for (int j = 0; j < F.getDimension(); ++j) {
			F[j] = log(F[j] + 1);
		}
		middle = middle + F;
		Table.push_back(F);
	}
	saveSetVectorsToCSV(Table, "C:\\Users\\troit\\OneDrive\\Рабочий стол\\Science_result\\SoFA_Var_Mod.csv", "SoFA_Var_Mod");
	Point middle1 = middle * 0.1;
	plt::named_plot(plot_name, Iteration, middle1.getPoint());
	plt::legend();
}
std::vector<long double> SoFAM_Grad_Numb(int N, int T, long double (*J)(Point), double (*Phi)(size_t), double (*Sigma)(size_t), int dim, double a, double b, double Q, size_t M, std::vector<int>& X, long double (*Distr)(long double, double, double, double)) {
	std::vector<long double> Variance(dim, Sigma(0));
	//int m = 0;
	//double P;
	int min_index = 0;
	std::vector<double> ref_numb;
	int numb_iteration = 0;
	std::vector<double> FirstPoint;
	std::vector<double> SecondPoint;
	std::vector<long double> Fitness;
	std::vector<long double> exp_Fitness;

	//std::vector<std::vector<double>> Population=CreatePopulation(N,dim,a,b);
	Population Pop = Population(N, dim);

	numb_iteration += N;
	long double Fitness_min = std::numeric_limits<long double>::max();
	Pop.GeneratePopulation(a, b);
	std::vector<std::vector<DistIndex>> DI = CreatePopulationSet(Pop);
	for (int i = 0; i < N; ++i) {
		Fitness.push_back(J(Pop[i]));
		exp_Fitness.push_back(exponential_Fitness(Fitness[i], dim));
		if (Fitness[i] < Fitness_min) {
			Fitness_min = Fitness[i];
			min_index = i;
		}
	}
	std::vector<long double> Y;
	int Ref_index;
	//Ref_index = ReferencePoint(Phi, exp_Fitness, 0);
	long double Fitness_candidate;
	std::vector<int> Ref_list;
	long double noize;
	long double noize_middle;
	std::vector<int> noize_iter;
	std::vector<long double> noize_list;
	std::vector<long double> sigma_list;
	for (int i = 0; i < T; i++) {
		Ref_index = ReferencePoint(Phi, exp_Fitness, i);
		//std::cout << m << std::endl;
		ref_numb.push_back(Ref_index);
		//double q = my_random_d(0, 1);
		Point Further(dim);
		Point Further_Swarm(dim);
		std::vector<DistIndex> D = DI[Ref_index];
		std::vector<DistIndex> DIK = getKNearest(D, dim);
		std::vector <Point> MatrToSolv = MatrixToSolve(Pop[Ref_index], D, Pop);
		Point SystPoint;
		std::vector<long double> B;
		for (int k = 0; k < DIK.size(); ++k) {
			SystPoint = Pop[Ref_index] - Pop[DIK[k].getIndex()];
			long double Norm = SystPoint.norm();
			B.push_back(abs(Fitness[Ref_index] - Fitness[DIK[k].getIndex()]) / Norm);
		}
		Point Ref_Grad = solveLinearSystem(MatrToSolv, B);
		Ref_Grad.Normalization();
		noize = 0.0;
		double k = 0.1;
		//while (true) {
		for (int j = 0; j < dim; ++j) {

			double q = my_random_d(0, 1);

			if (q < Q) {

				Variance[j] = 2 * Sigma(i);
				Further_Swarm[j] = Pop[Ref_index][j] - k * Ref_Grad[j];
				//Further[j] = Distr(Further_Swarm[j], a, b, Variance);
				//if (!((i + 9) % 5000))
					//noize += fabs(Further[j] - Pop[Ref_index][j]);
			}
			else {
				Further_Swarm[j] = Pop[Ref_index][j];
				Variance[j] = 2 * Sigma(i + M);
				//Further[j] = Distr(Further_Swarm[j], a, b, Variance);
				//if (!((i + 9) % 5000))
					//noize += fabs(Further[j] - Pop[Ref_index][j]);
			}
		}
		/*if (J(Further_Swarm) > Fitness[Ref_index]) {
			k *= 0.5;
			numb_iteration++;
			continue;

		}
		else {
			break;
		}*/
		//}
		for (int j = 0; j < dim; ++j) {
			Further[j] = Distr(Further_Swarm[j], a, b, Variance[j]);
		}
		/*if (!((i + 9) % 5000)) {
			sigma_list.push_back(Variance);
			noize_middle = noize / dim;
			noize_list.push_back(noize_middle);
			noize_iter.push_back(i);

		}*/
		Fitness_candidate = J(Further);
		numb_iteration++;
		if (Fitness_candidate < Fitness[Ref_index]) {
			Pop[Ref_index] = Further;
			Fitness[Ref_index] = Fitness_candidate;
			exp_Fitness[Ref_index] = exponential_Fitness(Fitness_candidate, dim);

			if (Fitness_candidate < Fitness_min) {
				Fitness_min = Fitness_candidate;
				min_index = Ref_index;
			}
		}
		X.push_back(numb_iteration);
		//Y.push_back(exponential_Fitness(Fitness_min,dim));
		Y.push_back(Fitness_min);
	}
	/*plt::plot(sigma_list, noize_list, ".");

		plt::xlabel("Sigma");
		plt::ylabel("Middle_Noize");
		plt::show();*/
		//std::vector<double> colors;
	return Y;
}
void MiddleSM_Grad_Numb(int N, int T, long double (*J)(Point), double (*Phi)(size_t), double (*Sigma)(size_t),
	int dim, double a, double b, double Q, int M, long double(*Distr)(long double, double, double, double),
	const char* plot_name, bool time_flag = false) {

	std::vector<int>Iteration;
	Point middle(T);
	Point F(T);
	std::vector<std::vector<long double>> Table;
	for (size_t i = 0; i < 10; ++i) {
		Iteration = {};
		/*if (time_flag) {
			F = SoFA_Var_Mod_Numb_Time(N, T, J, Phi, Sigma, dim, a, b, Q, M, Iteration, Distr);
		}
		else {*/
		F = SoFAM_Grad_Numb(N, T, J, Phi, Sigma, dim, a, b, Q, M, Iteration, Distr);
		//}
		for (int j = 0; j < F.getDimension(); ++j) {
			F[j] = log(F[j] + 1);
		}
		middle = middle + F;

		Table.push_back(F);
	}
	saveSetVectorsToCSV(Table, "C:\\Users\\troit\\OneDrive\\Рабочий стол\\Science_result\\SoFA_Var_Mod.csv", "SoFA_Var_Mod");
	Point middle1 = middle * 0.1;
	plt::named_plot(plot_name, Iteration, middle1.getPoint());
	plt::legend();
}

std::vector<long double> SoFA2(int N, int T, long double (*J)(Point), double (*Phi)(size_t), double (*Sigma)(size_t), double a, double b, double Q, size_t M, long double (*Distr)(long double, double, double, double)) {
	double Variance = Sigma(0);
	//int m = 0;
	//double P;
	int min_index = 0;
	std::vector<double> ref_numb;
	int numb_iteration = 0;
	std::vector<double> FirstPoint;
	std::vector<double> SecondPoint;
	std::vector<long double> Fitness;
	std::vector<long double> exp_Fitness;
	std::vector<long double> X1;
	std::vector<long double> X2;
	std::vector<long double> R1(1);
	std::vector<long double> R2(1);

	//std::vector<std::vector<double>> Population=CreatePopulation(N,dim,a,b);
	Population Pop = Population(N, 2);

	numb_iteration += N;
	long double Fitness_min = std::numeric_limits<long double>::max();
	Pop.GeneratePopulation(a, b);
	for (int i = 0; i < N; ++i) {

		Fitness.push_back(J(Pop[i]));
		exp_Fitness.push_back(exponential_Fitness(Fitness[i], 2));
		if (Fitness[i] < Fitness_min) {
			Fitness_min = Fitness[i];
			min_index = i;
		}

	}
	std::vector<long double> Y;
	int Ref_index;
	//Ref_index = ReferencePoint(Phi, exp_Fitness, 0);
	long double Fitness_candidate;
	std::vector<int> Ref_list;
	long double noize;
	long double noize_middle;
	std::vector<int> noize_iter;
	std::vector<long double> noize_list;
	std::vector<long double> sigma_list;
	for (int i = 0; i < T; i++) {
		X1.clear();
		X2.clear();
		R1.clear();
		R2.clear();
		Ref_index = ReferencePoint(Phi, exp_Fitness, i);
		//std::cout << m << std::endl;
		ref_numb.push_back(Ref_index);
		//double q = my_random_d(0, 1);
		Point Further(2);

		noize = 0.0;
		for (int j = 0; j < 2; ++j) {

			double q = my_random_d(0, 1);
			if (q < Q) {
				Variance = 2 * Sigma(i);
				Further[j] = Distr(Pop[Ref_index][j], a, b, Variance);
				if (!((i + 9) % 5000))
					noize += fabs(Further[j] - Pop[Ref_index][j]);
			}
			else {
				Variance = 2 * Sigma(i + M);
				Further[j] = Distr(Pop[Ref_index][j], a, b, Variance);
				if (!((i + 9) % 5000))
					noize += fabs(Further[j] - Pop[Ref_index][j]);
			}
		}

		if (!((i + 9) % 5000)) {
			sigma_list.push_back(Variance);
			noize_middle = noize / 2;
			noize_list.push_back(noize_middle);
			noize_iter.push_back(i);

		}
		Fitness_candidate = J(Further);
		numb_iteration++;
		if (Fitness_candidate < Fitness[Ref_index])
		{
			R1.push_back(Pop[Ref_index][0]);
			R2.push_back(Pop[Ref_index][1]);
			plt::plot(X1, X2, ".");
			plt::plot(R1, R2, "r.");
			R1.clear();
			R2.clear();

			//plt::pause(1);
			Pop[Ref_index] = Further;
			R1.push_back(Pop[Ref_index][0]);
			R2.push_back(Pop[Ref_index][1]);
			//plt::plot(R1, R2, ".");
			Fitness[Ref_index] = Fitness_candidate;
			exp_Fitness[Ref_index] = exponential_Fitness(Fitness_candidate, 2);
			for (int j = 0; j < N; ++j) {
				X1.push_back(Pop[j][0]);
				X2.push_back(Pop[j][1]);

			}
			std::cout << i << std::endl;
			plt::pause(1);
			plt::cla();
			plt::plot(X1, X2, ".");
			//plt::plot( R1,R2, "." );
			plt::pause(1);
			//plt::show();
			//plt::plot(X1, X2, ".");
			//plt::show();
			if (Fitness_candidate < Fitness_min) {
				Fitness_min = Fitness_candidate;
				min_index = Ref_index;
			}
		}
		//X.push_back(numb_iteration);
		//Y.push_back(exponential_Fitness(Fitness_min,dim));
		Y.push_back(Fitness_min);


	}
	/*plt::plot(sigma_list, noize_list, ".");

		plt::xlabel("Sigma");
		plt::ylabel("Middle_Noize");
		plt::show();*/
		/*std::vector<double> colors;
		for (size_t j = 0; j < FirstPoint.size(); ++j) {
			colors.push_back(j); // или любое другое значение для цвета
		}
		const char* plot_name = "Go";
		std::vector<double> color
		plt::scatter(FirstPoint, SecondPoint,50); // 50 - размер точек

		plt::colorbar();
		plt::show();*/
		/*plt::plot(ref_numb, X, ".");
		plt::xlabel("Reference point index");
		plt::ylabel("Iteration");
		plt::show();
		plt::close();*/
	return Y;
}