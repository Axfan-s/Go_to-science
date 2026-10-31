#pragma once
#include "Structures.h"
#include "matplotlibcpp.h"
#include <algorithm>
namespace plt = matplotlibcpp;
Point Mutation(Population Pop, const char* MutName, double F_i, long double a, long double b, int cur_index, int best_index = 0) {
	int dim = Pop[0].getDimension();
	Point Candidate(dim);
	if (MutName == "DE_Rand1") {
		int first_index = my_random_i(Pop.getSize() - 1);
		int second_index = my_random_i(Pop.getSize() - 1);
		int third_index = my_random_i(Pop.getSize() - 1);
		while (first_index == cur_index) {
			first_index = my_random_i(Pop.getSize() - 1);
		}
		while ((second_index == first_index) || (second_index == cur_index)) {
			second_index = my_random_i(Pop.getSize() - 1);
		}
		while ((third_index == second_index) || (third_index == first_index) || (third_index == cur_index)) {
			third_index = my_random_i(Pop.getSize() - 1);
		}

		for (int i = 0; i < dim; ++i) {
			double UpdateF = my_random_d();
			if (UpdateF <= 0.1) {
				F_i = my_random_d();
			}
			//std::cout << first_index << " " << second_index << " " << third_index << std::endl;;
			Candidate[i] = Pop[first_index][i] + F_i * (Pop[second_index][i] - Pop[third_index][i]);
			Candidate[i] = std::clamp(Candidate[i], a, b);
		}
	}
	else if (MutName == "DE_Best1") {
		int first_index = my_random_i(Pop.getSize() - 1);
		int second_index = my_random_i(Pop.getSize() - 1);
		while (first_index == cur_index) {
			first_index = my_random_i(Pop.getSize() - 1);
		}
		while ((second_index == first_index) || (second_index == cur_index)) {
			second_index = my_random_i(Pop.getSize() - 1);
		}
		for (int i = 0; i < dim; ++i) {
			double UpdateF = my_random_d();
			if (UpdateF <= 0.1) {
				F_i = my_random_d();
			}
			//std::cout << first_index << " " << second_index << " " << third_index << std::endl;;
			Candidate[i] = Pop[best_index][i] + F_i * (Pop[first_index][i] - Pop[second_index][i]);
			//std::cout << Pop[first_index][i]  <<" "<< Pop[second_index][i]<< Pop[best_index][i]<< std::endl;
			Candidate[i] = std::clamp(Candidate[i], a, b);
		}
	}
	return Candidate;
}

Point Crossingover(Point Cand, Point Cur, double CR) {
	for (int i = 0; i < Cand.getDimension(); ++i) {
		double ran = my_random_d();
		int jr = my_random_i(0, Cur.getDimension());
		if (ran > CR && i != jr) {
			Cand[i] = Cur[i];
		}
		double UpdateCR = my_random_d();
		if (UpdateCR <= 0.1) {
			CR = my_random_d();
		}
	}

	/*for (int i = 0; i < Cur.getDimension(); ++i) {
		if ( i != jr) {
			Cand[i] = Cur[i];
		}
	}*/
	return Cand;
}

std::vector<long double> DifferencialEvolution(long double(*J)(Point), int T, int dim, int N, double a, double b, const char* MutName) {
	Population Pop(N, dim);
	Pop.GeneratePopulation(a, b);
	double F_i = 0.5;
	double CR = 0.9;
	int best_index = 0;
	std::vector<long double> Fitness_Best;
	std::vector<long double> Fitness(N);
	for (int i = 0; i < N; ++i) {
		Fitness[i] = J(Pop[i]);
		if (Fitness[i] < Fitness[best_index]) {
			best_index = i;
		}
	}
	for (int t = 1; t < T; ++t) {
		for (int i = 0; i < N; ++i) {
			Point Mut = Mutation(Pop, MutName, F_i, a, b, i, best_index);
			Point Cand = Crossingover(Mut, Pop[i], CR);
			long double Fitness_Cand = J(Cand);
			if (Fitness_Cand < Fitness[i]) {
				Fitness[i] = Fitness_Cand;
				Pop[i] = Cand;
				if (Fitness_Cand < Fitness[best_index]) {
					best_index = i;
				}
			}
		}
		Fitness_Best.push_back(Fitness[best_index]);
	}
	for (int i = 0; i < dim; ++i) {
		//std::cout << Pop[best_index][i]<<" ";
	}
	//std::cout<<std::endl;
	return Fitness_Best;
}

void MiddleDE(
	long double (*J)(Point),
	int T,           // число поколений
	int dim,         // размерность пространства
	int N,           // размер популяции
	double a, double b,
	const char* MutName
) {
	const int NUM_RUNS = 10;
	std::vector<double> Fitness_Total(T - 1, 0.0);

	for (int run = 0; run < NUM_RUNS; ++run) {
		std::vector<long double> Fitness_cur = DifferencialEvolution(
			J, T, dim, N, a, b, MutName
		);

		// Проверка размера
		size_t hist_size = Fitness_cur.size();
		for (int gen = 0; gen < T - 1 && gen < hist_size; ++gen) {
			Fitness_Total[gen] += static_cast<double>(log(Fitness_cur[gen] + 1));
		}
	}

	// Усреднение
	for (int gen = 0; gen < T - 1; ++gen) {
		Fitness_Total[gen] /= NUM_RUNS;
	}

	// Формирование оси X (номера поколений)
	std::vector<int> generations(T - 1);
	for (int gen = 0; gen < T - 1; ++gen) {
		generations[gen] = gen * N;
	}

	// Построение графика
	plt::named_plot(MutName, generations, Fitness_Total);
	plt::legend();
	plt::grid(true);
}