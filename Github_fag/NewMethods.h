#pragma once
#include "Mutation.h"
#include "SoFA.h"
#include "matplotlibcpp.h"
//#include "Local_Search.h"
#include <algorithm>
#include <omp.h>
#include "DiffEval.h"
namespace plt = matplotlibcpp;
long double Give_Sum_Probability(double (*Phi)(size_t), std::vector<long double> Fitness, int k, std::vector<long double>& poweredF) {
	double r = my_random_d();
	long double sum = 0;
	int m = 0;
	long double ref = Fitness[0];
	//std::vector<long double> poweredF;
	for (int i = 0; i < Fitness.size(); ++i) {
		poweredF[i] = pow(Fitness[i], Phi(k));
		sum += poweredF[i];
	}

	return  sum;
}
int ReferencePoint_New(double sum, std::vector<long double> poweredF) {
	double sum_p = 0;
	double r = my_random_d();
	int m = 0;
	for (int i = 0; i < poweredF.size(); ++i) {
		sum_p += (poweredF[i]) / sum;
		if (sum_p > r) {
			m = i;
			break;
		}

	}
	return  m;
}
std::vector<long double> SoFAMS(int N, int T, long double (*J)(Point), double (*Phi)(size_t), double (*Sigma)(size_t), int dim, double a, double b, std::vector<int>& iteration_counts, long double(*Distr)(long double, double, double, double)) {
	double Variance = Sigma(0);
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
	X.GeneratePopulation(a, b);

	std::vector<long double> fitness(N);
	for (size_t i = 0; i < N; ++i) {
		fitness[i] = J(X[i]);
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

	// Íàõîæäåíèå ëó÷øèõ â êàæäîé ãðóïïå
	size_t best_m_idx = 0;
	long double best_m_fit = fitness_m[0];
	size_t worst_m_idx = 0;
	long double worst_m_fit = fitness_m[0];

	for (size_t i = 0; i < Nm; ++i) {
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

	for (size_t i = 0; i < Nf; ++i) {
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
		std::vector<long double> poweredF;
		if (Q < threshold) {
			// ========== ÔÀÇÀ ÈÑÑËÅÄÎÂÀÍÈß (áåç ïèùè) ==========

			// Ìóæñêèå îñîáè
			//for (size_t i = 0; i < Nm; ++i) {

			std::vector<long double> Fitness_Copy_m = Zero_Copy(fitness_m, fitness_m.size());
			for (size_t k = 0; k < Nm; ++k) {
				long double sum_prob = Give_Sum_Probability(Phi, Fitness_Copy_m, t, poweredF);
				int i = ReferencePoint_New(sum_prob, poweredF);
				size_t rand_idx = static_cast<size_t>(my_random_i(Nm - 1));
				int flag = vec_flag[my_random_i(1)];

				long double Am = std::exp(-fitness_m[rand_idx] / (fitness_m[i] + eps));

				for (size_t d = 0; d < dim; ++d) {
					long double r = my_random_d();
					X_new_m[i][d] = Xm[rand_idx][d] +
						flag * C2 * Am * ((b - a) * r + a);
				}
				sum_prob -= poweredF[i];
				poweredF = Zero_Copy(poweredF, i);
				Fitness_Copy_m = Zero_Copy(Fitness_Copy_m, i);

				//}

				// Æåíñêèå îñîáè
			}

			poweredF.clear();
			std::vector<long double> Fitness_Copy_f = Zero_Copy(fitness_f, fitness_f.size());
			for (size_t k = 0; k < Nf; ++k) {
				long double sum_prob = Give_Sum_Probability(Phi, Fitness_Copy_f, t, poweredF);
				int i = ReferencePoint_New(sum_prob, poweredF);
				size_t rand_idx = static_cast<size_t>(my_random_i(Nf - 1));
				int flag = vec_flag[my_random_i(1)];

				long double Af = std::exp(-fitness_f[rand_idx] / (fitness_f[i] + eps));
				//std::cout << i << ":" << "   ";
				for (size_t d = 0; d < dim; ++d) {
					long double r = my_random_d();
					X_new_f[i][d] = Xf[rand_idx][d] +
						flag * C2 * Af * ((b - a) * r + a);

					//std::cout << X_new_f[i][d] << " ";
				}
				sum_prob -= poweredF[i];
				poweredF = Zero_Copy(poweredF, i);
				Fitness_Copy_f = Zero_Copy(Fitness_Copy_f, i);
				//std::cout << std::endl;
			}
			poweredF.clear();
		}
		else {
			// ========== ÔÀÇÀ ÝÊÑÏËÓÀÒÀÖÈÈ (åñòü ïèùà) ==========

			if (temp > threshold2) {
				// Ãîðÿ÷àÿ ñðåäà
				std::vector<long double> Fitness_Copy_m = Zero_Copy(fitness_m, fitness_m.size());
				for (size_t k = 0; k < Nm; ++k) {
					long double sum_prob = Give_Sum_Probability(Phi, Fitness_Copy_m, t, poweredF);
					int i = ReferencePoint_New(sum_prob, poweredF);
					int flag = vec_flag[my_random_i(0, 1000) % 2];
					for (size_t d = 0; d < dim; ++d) {
						long double r = my_random_d();
						X_new_m[i][d] = X_food[d] +
							C3 * flag * temp * r * (X_food[d] - Xm[i][d]);
					}
					sum_prob -= poweredF[i];
					poweredF = Zero_Copy(poweredF, i);
					Fitness_Copy_m = Zero_Copy(Fitness_Copy_m, i);
				}
				poweredF.clear();

				std::vector<long double> Fitness_Copy_f = Zero_Copy(fitness_f, fitness_f.size());
				for (size_t k = 0; k < Nf; ++k) {
					long double sum_prob = Give_Sum_Probability(Phi, Fitness_Copy_f, t, poweredF);
					int i = ReferencePoint_New(sum_prob, poweredF);
					int flag = vec_flag[my_random_i(0, 1000) % 2];
					for (size_t d = 0; d < dim; ++d) {
						long double r = my_random_d();
						X_new_f[i][d] = X_food[d] +
							flag * C3 * temp * r * (X_food[d] - Xf[i][d]);
					}
					sum_prob -= poweredF[i];
					poweredF = Zero_Copy(poweredF, i);
					Fitness_Copy_f = Zero_Copy(Fitness_Copy_f, i);
				}
				poweredF.clear();
			}
			else {
				// Õîëîäíàÿ ñðåäà
				if (my_random_d() > 0.6) {
					// Ðåæèì áîðüáû

					std::vector<long double> Fitness_Copy_m = Zero_Copy(fitness_m, fitness_m.size());
					for (size_t k = 0; k < Nm; ++k) {
						long double sum_prob = Give_Sum_Probability(Phi, Fitness_Copy_m, t, poweredF);
						int i = ReferencePoint_New(sum_prob, poweredF);
						for (size_t d = 0; d < dim; ++d) {
							long double r = my_random_d();
							long double FM = std::exp(-best_f_fit / (fitness_m[i] + eps));
							X_new_m[i][d] = Xm[i][d] +
								C3 * FM * r * (Q * X_best_f[d] - Xm[i][d]);
						}
						sum_prob -= poweredF[i];
						poweredF = Zero_Copy(poweredF, i);
						Fitness_Copy_m = Zero_Copy(Fitness_Copy_m, i);
					}
					poweredF.clear();

					std::vector<long double> Fitness_Copy_f = Zero_Copy(fitness_f, fitness_f.size());
					for (size_t k = 0; k < Nf; ++k) {
						long double sum_prob = Give_Sum_Probability(Phi, Fitness_Copy_f, t, poweredF);
						int i = ReferencePoint_New(sum_prob, poweredF);
						for (size_t d = 0; d < dim; ++d) {
							long double r = my_random_d();
							long double FF = std::exp(-best_m_fit / (fitness_f[i] + eps));
							X_new_f[i][d] = Xf[i][d] +
								C3 * FF * r * (Q * X_best_m[d] - Xf[i][d]);
						}
						sum_prob -= poweredF[i];
						poweredF = Zero_Copy(poweredF, i);
						Fitness_Copy_f = Zero_Copy(Fitness_Copy_f, i);
					}
					poweredF.clear();
				}
				else {
					// Ðåæèì ñïàðèâàíèÿ
					std::vector<long double> Fitness_Copy_m = Zero_Copy(fitness_m, fitness_m.size());
					for (size_t k = 0; k < Nm; ++k) {
						long double sum_prob = Give_Sum_Probability(Phi, Fitness_Copy_m, t, poweredF);
						int i = ReferencePoint_New(sum_prob, poweredF);
						for (size_t d = 0; d < dim; ++d) {
							long double r = my_random_d();
							// Áåçîïàñíûé äîñòóï ê æåíñêèì îñîáÿì
							size_t f_idx = (i < Nf) ? i : Nf - 1;
							long double Mm = std::exp(-fitness_f[f_idx] / (fitness_m[i] + eps));
							X_new_m[i][d] = Xm[i][d] +
								C3 * r * Mm * (Q * Xf[f_idx][d] - Xm[i][d]);
						}
						sum_prob -= poweredF[i];
						poweredF = Zero_Copy(poweredF, i);
						Fitness_Copy_m = Zero_Copy(Fitness_Copy_m, i);
					}
					poweredF.clear();
					std::vector<long double> Fitness_Copy_f = Zero_Copy(fitness_f, fitness_f.size());
					for (size_t k = 0; k < Nf; ++k) {
						long double sum_prob = Give_Sum_Probability(Phi, Fitness_Copy_f, t, poweredF);
						int i = ReferencePoint_New(sum_prob, poweredF);
						for (size_t d = 0; d < dim; ++d) {
							long double r = my_random_d();
							// Áåçîïàñíûé äîñòóï ê ìóæñêèì îñîáÿì
							size_t m_idx = (i < Nm) ? i : Nm - 1;
							int flag = vec_flag[my_random_i(0, 1000) % 2];
							long double Mf = std::exp(-fitness_m[m_idx] / (fitness_f[i] + eps));
							X_new_f[i][d] = Xf[i][d] +
								C3 * flag * r * Mf * (Q * Xm[m_idx][d] - Xf[i][d]);
						}
						sum_prob -= poweredF[i];
						poweredF = Zero_Copy(poweredF, i);
						Fitness_Copy_f = Zero_Copy(Fitness_Copy_f, i);
					}

					// Âîçìîæíîñòü îòêëàäêè ÿèö (ñáðîñ õóäøèõ îñîáåé)
					/*if (vec_flag[my_random_i(0, 1000) % 2] == 1) {
						for (size_t d = 0; d < dim; ++d) {
							long double r = my_random_d();
							X_new_m[worst_m_idx][d] = a + r * (b - a);
							X_new_f[worst_f_idx][d] = a + r * (b - a);
						}
					}*/
				}
			}
		}

		// 4. Ïðîâåðêà ãðàíèö è îáíîâëåíèå ïîïóëÿöèé

		// Ìóæñêèå îñîáè
		for (size_t i = 0; i < Nm; ++i) {
			// Ïðîâåðêà ãðàíèö
			for (size_t d = 0; d < dim; ++d) {

				if (X_new_m[i][d] > b) X_new_m[i][d] = b;
				if (X_new_m[i][d] < a) X_new_m[i][d] = a;
				X_new_m[i][d] = Distr(X_new_m[i][d], a, b, Sigma(t));
			}

			// Îöåíêà íîâîé ïîçèöèè
			long double new_fit = J(X_new_m[i]);
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
		//std::cout<<"Mutation" << std::endl;
		for (size_t i = 0; i < Nf; ++i) {
			// Ïðîâåðêà ãðàíèö

			for (size_t d = 0; d < dim; ++d) {
				if (X_new_f[i][d] > b) X_new_f[i][d] = b;
				if (X_new_f[i][d] < a) X_new_f[i][d] = a;

				X_new_f[i][d] = Distr(X_new_f[i][d], a, b, Sigma(t));
				//std::cout << X_new_f[i][d] << " ";
			}
			//std::cout << std::endl;

			// Îöåíêà íîâîé ïîçèöèè
			long double new_fit = J(X_new_f[i]);
			total_evaluations++;

			if (new_fit < fitness_f[i]) {
				fitness_f[i] = new_fit;
				Xf[i] = X_new_f[i];
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

		best_history[t] = global_best_fit;

		iteration_counts.push_back(total_evaluations);
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
void MiddleSMS(
	size_t N,
	size_t T,
	long double (*objective_func)(Point),
	double (*Phi)(size_t),
	double (*Sigma)(size_t),
	size_t dim,
	double lb,
	double ub,


	double Q,
	int M,
	long double(*Distr)(long double, double, double, double),
	const std::string& plot_name,
	bool time_flag = false,

	size_t num_trials = 10)
{
	std::vector<int> iteration_counts;
	std::vector<double> average_history(T, 0.0);
	std::vector<std::vector<long double>> Table;
	for (size_t trial = 0; trial < num_trials; ++trial) {
		std::vector<int> current_iterations;
		std::vector<long double> history = SoFAMS(
			N, T, objective_func, Phi, Sigma, dim, lb, ub, current_iterations, Distr);
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


std::vector<long double> SoFASM(int N, int T, long double (*J)(Point), double (*Phi)(size_t), double (*Sigma)(size_t), int dim, double a, double b, double Q, size_t M, std::vector<int>& X, long double(*Distr)(long double, double, double, double)) {
	double Variance = Sigma(0);
	size_t m = 0;
	//double P;
	std::vector<long double> ref_numb;
	int numb_iteration = 0;
	double Threshold = 0.25;
	double Thresold2 = 0.6;
	double C1 = 0.5;
	double C2 = .05;
	double C3 = 2;
	double Temp;
	double tt;
	double A;
	int flag_index;
	int Flag;
	double ran2;
	double L;
	std::vector<long double> exp_Fitness;
	long double Fitness_candidate;
	long double Fitness_min = INT_MAX;
	std::vector<long double> Fitness;
	std::vector<int> vec_flag = { 1, -1 };
	int rand_leader_index;
	int min_index = 0;
	//std::vector<std::vector<double>> Population=CreatePopulation(N,dim,a,b);
	Population Pop = Population(N, dim);
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

	for (int i = 0; i < T; i++) {
		Temp = std::exp(-i * 1.0 / T);
		tt = (T - i) * 1.0 / T;  // Будет от 1.0 до 0.0
		L = C1 * std::exp(-tt);
		if (L > 1)
			L = 1;
		Ref_index = ReferencePoint(Phi, exp_Fitness, i);
		Point Further_Swarm(dim);

		//std::cout << m << std::endl;
		ref_numb.push_back(Ref_index);


		Point Further(dim);
		if (i == 1000) {
			int r = 3;
		}

		/*for (size_t j = 0; j < Ref.getDimension(); j++) {
			double eta = my_rand() % 1001 * 1.0 / 1000;
			double Expr = eta * (boost::math::erf((b - Ref[j]) / sqrt(2 * Variance))) + (1 - eta) * boost::math::erf((a - Ref[j]) / sqrt(2 * Variance));
			double R;
			if (Expr >= 1) {
				R = b;
			}
			else if (Expr <= -1) {
				R = a;
			}
			else {
				R = sqrt(2 * Variance) * (boost::math::erf_inv(Expr)) + Ref[j];
			}
			//double R = Mutation(Ref[j], a, b, Variance);
			Further[j]=R;
		}*/
		for (int j = 0; j < dim; ++j) {
			double q = my_random_d(0, 1);
			if (q < Q) {
				Variance = 2 * Sigma(i);
				if (L < Threshold) {
					rand_leader_index = my_random_i(N);
					double denominator = std::max(Fitness[Ref_index], 1e-300L);
					A = std::exp(-Fitness[rand_leader_index] / denominator);

					flag_index = (my_random_i(2));
					Flag = vec_flag[flag_index];
					ran2 = my_random_d(0, 1);
					Further_Swarm[j] = std::min(b * 1.0L, (std::max(a * 1.0L, Pop[Ref_index][j] + H(i) * (Flag * C2 * A * ((b - a) * ran2 + a)))));
					Further[j] = Distr(Further_Swarm[j], a, b, Variance);
				}
				else {
					if (Temp > Thresold2) {
						rand_leader_index = my_random_i(N);
						flag_index = my_random_i(2);
						Flag = vec_flag[flag_index];
						ran2 = my_random_d();
						Further_Swarm[j] = std::min(b * 1.0L, (std::max(a * 1.0L, Pop[min_index][j] + ((Pop[min_index][j] - Pop[Ref_index][j]) * C3 * Flag * Temp * ran2) * H(i))));
						Further[j] = Distr(Further_Swarm[j], a, b, Variance);

					}
					else {
						Further[j] = Distr(Pop[Ref_index][j], a, b, Variance);

					}
				}
			}
			else {
				Variance = 2 * Sigma(i + M);
				Further[j] = Distr(Pop[Ref_index][j], a, b, Variance);
			}
		}
		Fitness_candidate = J(Further);
		numb_iteration++;
		if (!(i % 7000000)) {
			std::vector<std::vector<long double>> Table = Pop.PopulationToMatrix();
			saveSetVectorsToCSV(Table, "C:\\Users\\troit\\OneDrive\\Рабочий стол\\Science_result\\SoFASM_Pop.csv", "SoFASM_Pop");
		}
		if (Fitness_candidate < Fitness[Ref_index]) {
			Fitness[Ref_index] = Fitness_candidate;
			Pop[Ref_index] = Further;

			exp_Fitness[Ref_index] = exponential_Fitness(Fitness_candidate, dim);
			if (Fitness_candidate < Fitness[min_index]) {
				min_index = Ref_index;

			}
		}

		X.push_back(numb_iteration);
		Y.push_back(Fitness[min_index]);
	}
	return Y;
}

std::vector<long double> SoFASM_Time(int N, int T, long double (*J)(Point), double (*Phi)(size_t), double (*Sigma)(size_t), double (*H)(size_t), int dim, double a, double b, double Q, size_t M, std::vector<int>& X, long double(*Distr)(long double, double, double, double)) {
	double Variance = Sigma(0);
	size_t m = 0;
	//double P;
	std::vector<long double> ref_numb;
	int numb_iteration = 0;
	double Threshold = 0.25;
	double Thresold2 = 0.6;
	double C1 = 0.5;
	double C2 = .05;
	double C3 = 2;
	double Temp;
	double tt;
	double A;
	int flag_index;
	int Flag;
	double ran2;
	double L;
	int start = clock();
	std::vector<long double> exp_Fitness;
	long double Fitness_candidate;
	long double Fitness_min = INT_MAX;
	std::vector<long double> Fitness;
	std::vector<int> vec_flag = { 1, -1 };
	int rand_leader_index;
	int min_index = 0;
	//std::vector<std::vector<double>> Population=CreatePopulation(N,dim,a,b);
	Population Pop = Population(N, dim);
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

	for (int i = 0; i < T; i++) {
		Temp = std::exp(-i * 1.0 / T);
		tt = (T - i) * 1.0 / T;  // Будет от 1.0 до 0.0

		L = C1 * std::exp(-tt);
		if (L > 1)
			L = 1;
		Ref_index = ReferencePoint(Phi, exp_Fitness, i);
		Point Further_Swarm(dim);

		//std::cout << m << std::endl;
		ref_numb.push_back(Ref_index);


		Point Further(dim);
		if (i == 1000) {
			int r = 3;
		}

		/*for (size_t j = 0; j < Ref.getDimension(); j++) {
			double eta = my_rand() % 1001 * 1.0 / 1000;
			double Expr = eta * (boost::math::erf((b - Ref[j]) / sqrt(2 * Variance))) + (1 - eta) * boost::math::erf((a - Ref[j]) / sqrt(2 * Variance));
			double R;
			if (Expr >= 1) {
				R = b;
			}
			else if (Expr <= -1) {
				R = a;
			}
			else {
				R = sqrt(2 * Variance) * (boost::math::erf_inv(Expr)) + Ref[j];
			}
			//double R = Mutation(Ref[j], a, b, Variance);
			Further[j]=R;
		}*/
		for (int j = 0; j < dim; ++j) {
			double q = my_random_d(0, 1);
			if (q < Q) {
				Variance = 2 * Sigma(i);
				if (L < Threshold) {
					rand_leader_index = my_random_i(N);
					double denominator = std::max(Fitness[Ref_index], 1e-300L);
					A = std::exp(-Fitness[rand_leader_index] / denominator);

					flag_index = (my_random_i(2));
					Flag = vec_flag[flag_index];
					ran2 = my_random_d(0, 1);
					Further_Swarm[j] = std::min(b * 1.0L, (std::max(a * 1.0L, Pop[Ref_index][j] + H(i) * (Flag * C2 * A * ((b - a) * ran2 + a)))));
					Further[j] = Distr(Further_Swarm[j], a, b, Variance);
				}
				else {
					if (Temp > Thresold2) {
						rand_leader_index = my_random_i(N);
						flag_index = my_random_i(2);
						Flag = vec_flag[flag_index];
						ran2 = my_random_d();
						Further_Swarm[j] = std::min(b * 1.0L, (std::max(a * 1.0L, Pop[min_index][j] + ((Pop[min_index][j] - Pop[Ref_index][j]) * C3 * Flag * Temp * ran2) * H(i))));
						Further[j] = Distr(Further_Swarm[j], a, b, Variance);

					}
					else {
						Further[j] = Distr(Pop[Ref_index][j], a, b, Variance);

					}
				}
			}
			else {
				Variance = 2 * Sigma(i + M);
				Further[j] = Distr(Pop[Ref_index][j], a, b, Variance);
			}
		}
		Fitness_candidate = J(Further);
		numb_iteration++;
		if (!(i % 7000000)) {
			std::vector<std::vector<long double>> Table = Pop.PopulationToMatrix();
			saveSetVectorsToCSV(Table, "C:\\Users\\troit\\OneDrive\\Рабочий стол\\Science_result\\SoFASM_Pop.csv", "SoFASM_Pop");
		}
		if (Fitness_candidate < Fitness[Ref_index]) {
			Fitness[Ref_index] = Fitness_candidate;
			Pop[Ref_index] = Further;

			exp_Fitness[Ref_index] = exponential_Fitness(Fitness_candidate, dim);
			if (Fitness_candidate < Fitness[min_index]) {
				min_index = Ref_index;

			}
		}

		X.push_back(numb_iteration);
		Y.push_back(Fitness[min_index]);
	}
	return Y;
}

void MiddleSSM(int N, int T, long double (*J)(Point), double (*Phi)(size_t), double (*Sigma)(size_t),
	int dim, double a, double b, double Q, int M, long double(*Distr)(long double, double, double, double),
	const char* plot_name, bool time_flag = false) {

	std::vector<int>Iteration;
	Point middle(T);
	Point F(T);
	for (size_t i = 0; i < 10; ++i) {
		Iteration = {};
		if (time_flag) {
			F = SoFASM(N, T, J, Phi, Sigma, dim, a, b, Q, M, Iteration, Distr);
		}
		else {
			F = SoFASM_Time(N, T, J, Phi, Sigma, H, dim, a, b, Q, M, Iteration, Distr);
		}
		for (int j = 0; j < F.getDimension(); ++j) {
			F[j] = log(F[j] + 1);
		}
		middle = middle + F;
	}
	Point middle1 = middle * 0.1;
	plt::named_plot(plot_name, Iteration, middle1.getPoint());
	plt::legend();

}


void MiddleSSM_parallel(int N, int T, long double (*J)(Point), double (*Phi)(size_t),
	double (*Sigma)(size_t), double (*H)(size_t),
	int dim, double a, double b, double Q, int M,
	long double(*Distr)(long double, double, double, double),
	const char* plot_name, bool time_flag = false) {

	const int NUM_RUNS = 10;
	std::vector<Point> results(NUM_RUNS, Point(T));  // для хранения результатов каждого запуска
	std::vector<std::vector<int>> iterations(NUM_RUNS);  // для хранения итераций каждого запуска

	// Параллельные запуски
#pragma omp parallel for num_threads(10)
	for (int run = 0; run < NUM_RUNS; ++run) {
		std::vector<int> current_iter;
		Point F;
		if (time_flag) {
			SoFASM_Time(N, T, J, Phi, Sigma, H, dim, a, b, Q, M, current_iter, Distr);
		}
		else {
			SoFASM(N, T, J, Phi, Sigma, dim, a, b, Q, M, current_iter, Distr);

		}

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
	std::cout << "\n";
	std::cout << "SoFASM " << exponential_Fitness(middle1[T - 1], dim);
	plt::named_plot(plot_name, iterations[NUM_RUNS - 1], middle1.getPoint());
	std::cout << std::endl;
	plt::legend();
}


std::vector<long double> SoFASM2(int N, int T, long double (*J)(Point), double (*Phi)(size_t), double (*Sigma)(size_t), double (*H)(size_t), int dim, double a, double b, double Q, size_t M, std::vector<int>& X, long double(*Distr)(long double, double, double, double)) {
	double Variance = Sigma(0);
	size_t m = 0;
	//double P;
	std::vector<long double> ref_numb;
	int numb_iteration = 0;
	double Threshold = 0.25;
	double Thresold2 = 0.6;
	double C1 = 0.5;
	double C2 = .05;
	double C3 = 2;
	double Temp;
	double tt;
	double A;
	int flag_index;
	int Flag;
	double ran2;
	double L;
	std::vector<long double> exp_Fitness;
	long double Fitness_candidate;
	long double Fitness_min = INT_MAX;
	std::vector<long double> Fitness;
	std::vector<int> vec_flag = { 1, -1 };
	int rand_leader_index;
	int min_index = 0;
	//std::vector<std::vector<double>> Population=CreatePopulation(N,dim,a,b);
	Population Pop = Population(N, dim);
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

	for (int i = 0; i < T; i++) {
		Temp = std::exp(-i * 1.0 / T);
		tt = (T - i) * 1.0 / T;  // Будет от 1.0 до 0.0
		L = C1 * std::exp(-tt);
		if (L > 1)
			L = 1;
		Ref_index = ReferencePoint(Phi, exp_Fitness, i);
		Point Further_Swarm(dim);

		//std::cout << m << std::endl;
		ref_numb.push_back(Ref_index);


		Point Further(dim);
		if (i == 1000) {
			int r = 3;
		}

		/*for (size_t j = 0; j < Ref.getDimension(); j++) {
			double eta = my_rand() % 1001 * 1.0 / 1000;
			double Expr = eta * (boost::math::erf((b - Ref[j]) / sqrt(2 * Variance))) + (1 - eta) * boost::math::erf((a - Ref[j]) / sqrt(2 * Variance));
			double R;
			if (Expr >= 1) {
				R = b;
			}
			else if (Expr <= -1) {
				R = a;
			}
			else {
				R = sqrt(2 * Variance) * (boost::math::erf_inv(Expr)) + Ref[j];
			}
			//double R = Mutation(Ref[j], a, b, Variance);
			Further[j]=R;
		}*/
		for (int j = 0; j < dim; ++j) {
			double q = my_random_d(0, 1);
			if (q < Q) {
				Variance = 2 * Sigma(i);
				if (L < Threshold) {
					rand_leader_index = my_random_i(N);
					double denominator = std::max(Fitness[Ref_index], 1e-300L);
					A = std::exp(-Fitness[rand_leader_index] / denominator);

					flag_index = (my_random_i(2));
					Flag = vec_flag[flag_index];
					ran2 = my_random_d(0, 1);
					Further_Swarm[j] = std::min(b * 1.0L, (std::max(a * 1.0L, Pop[min_index][j] + ((Pop[min_index][j] - Pop[Ref_index][j]) * C3 * Flag * Temp * ran2) * H(i))));

					Further[j] = Distr(Further_Swarm[j], a, b, Variance);
				}
				else {
					if (Temp > Thresold2) {
						rand_leader_index = my_random_i(N);
						flag_index = my_random_i(2);
						Flag = vec_flag[flag_index];
						ran2 = my_random_d();
						Further_Swarm[j] = std::min(b * 1.0L, (std::max(a * 1.0L, Pop[Ref_index][j] + H(i) * (Flag * C2 * A * ((b - a) * ran2 + a)))));
						Further[j] = Distr(Further_Swarm[j], a, b, Variance);

					}
					else {
						Further[j] = Distr(Pop[Ref_index][j], a, b, Variance);

					}
				}
			}
			else {
				Variance = 2 * Sigma(i + M);
				Further[j] = Distr(Pop[Ref_index][j], a, b, Variance);
			}
		}
		Fitness_candidate = J(Further);
		numb_iteration++;
		if (!(i % 7000000)) {
			std::vector<std::vector<long double>> Table = Pop.PopulationToMatrix();
			saveSetVectorsToCSV(Table, "C:\\Users\\troit\\OneDrive\\Рабочий стол\\Science_result\\SoFASM_Pop.csv", "SoFASM_Pop");
		}
		if (Fitness_candidate < Fitness[Ref_index]) {
			Fitness[Ref_index] = Fitness_candidate;
			Pop[Ref_index] = Further;

			exp_Fitness[Ref_index] = exponential_Fitness(Fitness_candidate, dim);
			if (Fitness_candidate < Fitness[min_index]) {
				min_index = Ref_index;

			}
		}

		X.push_back(numb_iteration);
		Y.push_back(Fitness[min_index]);
	}
	return Y;
}

void MiddleSSM2(int N, int T, long double (*J)(Point), double (*Phi)(size_t), double (*Sigma)(size_t), double (*H)(size_t),
	int dim, double a, double b, double Q, int M, long double(*Distr)(long double, double, double, double),
	const char* plot_name) {

	std::vector<int>Iteration;
	Point middle(T);
	Point F(T);
	for (size_t i = 0; i < 10; ++i) {
		Iteration = {};
		F = SoFASM2(N, T, J, Phi, Sigma, H, dim, a, b, Q, M, Iteration, Distr);
		middle = middle + F;
	}
	Point middle1 = middle * 0.1;
	plt::named_plot(plot_name, Iteration, middle1.getPoint());
	plt::legend();

}

void MiddleSSM2_parallel(int N, int T, long double (*J)(Point), double (*Phi)(size_t),
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
		Point F = SoFASM2(N, T, J, Phi, Sigma, H, dim, a, b, Q, M, current_iter, Distr);

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
	std::cout << "SoFASM2 " << exponential_Fitness(middle1[T - 1], dim);
	std::cout << std::endl;
	plt::named_plot(plot_name, iterations[NUM_RUNS - 1], middle1.getPoint());
	plt::legend();
}



std::vector<long double> SoFASM_without_H_and_Rand(int N, int T, long double (*J)(Point), double (*Phi)(size_t), double (*Sigma)(size_t), double (*H)(size_t), int dim, double a, double b, double Q, size_t M, std::vector<int>& X, long double(*Distr)(long double, double, double, double)) {
	double Variance = Sigma(0);
	size_t m = 0;
	//double P;
	std::vector<long double> ref_numb;
	int numb_iteration = 0;
	double Threshold = 0.25;
	double Thresold2 = 0.6;
	double C1 = 0.5;
	double C2 = .05;
	double C3 = 2;
	double Temp;
	double tt;
	double A;
	int flag_index;
	int Flag;
	double ran2;
	double L;
	std::vector<long double> exp_Fitness;
	long double Fitness_candidate;
	long double Fitness_min = INT_MAX;
	std::vector<long double> Fitness;
	std::vector<int> vec_flag = { 1, -1 };
	int rand_leader_index;
	int min_index = 0;
	//std::vector<std::vector<double>> Population=CreatePopulation(N,dim,a,b);
	Population Pop = Population(N, dim);
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

	for (int i = 0; i < T; i++) {
		Temp = std::exp(-i * 1.0 / T);
		tt = (T - i) * 1.0 / T;  // Будет от 1.0 до 0.0
		L = C1 * std::exp(-tt);
		if (L > 1)
			L = 1;
		Ref_index = ReferencePoint(Phi, exp_Fitness, i);
		Point Further_Swarm(dim);

		//std::cout << m << std::endl;
		ref_numb.push_back(Ref_index);


		Point Further(dim);


		/*for (size_t j = 0; j < Ref.getDimension(); j++) {
			double eta = my_rand() % 1001 * 1.0 / 1000;
			double Expr = eta * (boost::math::erf((b - Ref[j]) / sqrt(2 * Variance))) + (1 - eta) * boost::math::erf((a - Ref[j]) / sqrt(2 * Variance));
			double R;
			if (Expr >= 1) {
				R = b;
			}
			else if (Expr <= -1) {
				R = a;
			}
			else {
				R = sqrt(2 * Variance) * (boost::math::erf_inv(Expr)) + Ref[j];
			}
			//double R = Mutation(Ref[j], a, b, Variance);
			Further[j]=R;
		}*/
		for (int j = 0; j < dim; ++j) {
			double q = my_random_d(0, 1);
			if (q < Q) {
				Variance = 2 * Sigma(i);
				if (L < Threshold) {
					rand_leader_index = my_random_i(N);
					double denominator = std::max(Fitness[Ref_index], 1e-300L);
					A = std::exp(-Fitness[rand_leader_index] / denominator);

					flag_index = (my_random_i(2));
					Flag = vec_flag[flag_index];
					ran2 = my_random_d(0, 1);
					Further_Swarm[j] = std::min(b * 1.0L, (std::max(a * 1.0L, Pop[Ref_index][j] + (Flag * C2 * A * ((b - a) + a)))));
					Further[j] = Distr(Further_Swarm[j], a, b, Variance);
				}
				else {
					if (Temp > Thresold2) {
						rand_leader_index = my_random_i(N);
						//std::cout << "???????????" << std::endl;
						flag_index = my_random_i(2);
						Flag = vec_flag[flag_index];
						ran2 = my_random_d();
						Further_Swarm[j] = std::min(b * 1.0L, (std::max(a * 1.0L, Pop[min_index][j] + ((Pop[min_index][j] - Pop[Ref_index][j]) * C3 * Flag * Temp * ran2) * H(i))));
						Further[j] = Distr(Further_Swarm[j], a, b, Variance);

					}
					else {
						Further[j] = Distr(Pop[Ref_index][j], a, b, Variance);

					}
				}
			}
			else {
				Variance = 2 * Sigma(i + M);
				Further[j] = Distr(Pop[Ref_index][j], a, b, Variance);
			}
		}
		/*if (i % 10000 == 0) {
			std::cout << "Variance=" << Variance << std::endl;
			for (int l = 0; l < Further.getDimension(); ++l) {
				std::cout << Further_Swarm[l] << " ";
			}

			std::cout << std::endl;
			for (int l = 0; l < Further.getDimension(); ++l) {
				std::cout << Further[l] << " ";
							}
			std::cout << std::endl;
		}*/
		Fitness_candidate = J(Further);
		numb_iteration++;
		if (!(i % 7000000)) {
			std::vector<std::vector<long double>> Table = Pop.PopulationToMatrix();
			saveSetVectorsToCSV(Table, "C:\\Users\\troit\\OneDrive\\Рабочий стол\\Science_result\\SoFASM_Pop.csv", "SoFASM_Pop");
		}
		if (Fitness_candidate < Fitness[Ref_index]) {
			Fitness[Ref_index] = Fitness_candidate;
			Pop[Ref_index] = Further;

			exp_Fitness[Ref_index] = exponential_Fitness(Fitness_candidate, dim);
			if (Fitness_candidate < Fitness[min_index]) {
				min_index = Ref_index;

			}
		}

		X.push_back(numb_iteration);
		Y.push_back(Fitness[min_index]);
	}
	return Y;
}

void MiddleSSM_without_H_and_Rand(int N, int T, long double (*J)(Point), double (*Phi)(size_t), double (*Sigma)(size_t), double (*H)(size_t),
	int dim, double a, double b, double Q, int M, long double(*Distr)(long double, double, double, double),
	const char* plot_name) {

	std::vector<int>Iteration;
	Point middle(T);
	Point F(T);
	for (size_t i = 0; i < 10; ++i) {
		Iteration = {};
		F = SoFASM_without_H_and_Rand(N, T, J, Phi, Sigma, H, dim, a, b, Q, M, Iteration, Distr);
		middle = middle + F;
	}
	Point middle1 = middle * 0.1;
	plt::named_plot(plot_name, Iteration, middle1.getPoint());
	plt::legend();

}
void MiddleSSM_without_H_and_Rand_parallel(int N, int T, long double (*J)(Point), double (*Phi)(size_t),
	double (*Sigma)(size_t),
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
		Point F = SoFASM_without_H_and_Rand(N, T, J, Phi, Sigma, H, dim, a, b, Q, M, current_iter, Distr);

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
	std::cout << "SoFASM without H and rand " << exponential_Fitness(middle1[T - 1], dim);
	std::cout << std::endl;
	plt::named_plot(plot_name, iterations[NUM_RUNS - 1], middle1.getPoint());
	plt::legend();
}

std::vector<long double> SoFASM_without_H(int N, int T, long double (*J)(Point), double (*Phi)(size_t), double (*Sigma)(size_t), double (*H)(size_t), int dim, double a, double b, double Q, size_t M, std::vector<int>& X, long double(*Distr)(long double, double, double, double)) {
	double Variance = Sigma(0);
	size_t m = 0;
	//double P;
	std::vector<long double> ref_numb;
	int numb_iteration = 0;
	double Threshold = 0.25;
	double Thresold2 = 0.6;
	double C1 = 0.5;
	double C2 = .05;
	double C3 = 2;
	double Temp;
	double tt;
	double A;
	int flag_index;
	int Flag;
	double ran2;
	double L;
	std::vector<long double> exp_Fitness;
	long double Fitness_candidate;
	long double Fitness_min = INT_MAX;
	std::vector<long double> Fitness;
	std::vector<int> vec_flag = { 1, -1 };
	int rand_leader_index;
	int min_index = 0;
	//std::vector<std::vector<double>> Population=CreatePopulation(N,dim,a,b);
	Population Pop = Population(N, dim);
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

	for (int i = 0; i < T; i++) {
		Temp = std::exp(-i * 1.0 / T);
		tt = (T - i) * 1.0 / T;  // Будет от 1.0 до 0.0
		L = C1 * std::exp(-tt);
		if (L > 1)
			L = 1;
		Ref_index = ReferencePoint(Phi, exp_Fitness, i);
		Point Further_Swarm(dim);

		//std::cout << m << std::endl;
		ref_numb.push_back(Ref_index);


		Point Further(dim);
		if (i == 1000) {
			int r = 3;
		}

		/*for (size_t j = 0; j < Ref.getDimension(); j++) {
			double eta = my_rand() % 1001 * 1.0 / 1000;
			double Expr = eta * (boost::math::erf((b - Ref[j]) / sqrt(2 * Variance))) + (1 - eta) * boost::math::erf((a - Ref[j]) / sqrt(2 * Variance));
			double R;
			if (Expr >= 1) {
				R = b;
			}
			else if (Expr <= -1) {
				R = a;
			}
			else {
				R = sqrt(2 * Variance) * (boost::math::erf_inv(Expr)) + Ref[j];
			}
			//double R = Mutation(Ref[j], a, b, Variance);
			Further[j]=R;
		}*/
		for (int j = 0; j < dim; ++j) {
			double q = my_random_d(0, 1);
			if (q < Q) {
				Variance = 2 * Sigma(i);
				if (L < Threshold) {
					rand_leader_index = my_random_i(N);
					double denominator = std::max(Fitness[Ref_index], 1e-300L);
					A = std::exp(-Fitness[rand_leader_index] / denominator);

					flag_index = (my_random_i(2));
					Flag = vec_flag[flag_index];
					ran2 = my_random_d(0, 1);
					Further_Swarm[j] = std::min(b * 1.0L, (std::max(a * 1.0L, Pop[Ref_index][j] + (Flag * C2 * A * ((b - a) * ran2 + a)))));
					Further[j] = Distr(Further_Swarm[j], a, b, Variance);
				}
				else {
					if (Temp > Thresold2) {
						rand_leader_index = my_random_i(N);
						flag_index = my_random_i(2);
						Flag = vec_flag[flag_index];
						ran2 = my_random_d();
						Further_Swarm[j] = std::min(b * 1.0L, (std::max(a * 1.0L, Pop[min_index][j] + ((Pop[min_index][j] - Pop[Ref_index][j]) * C3 * Flag * Temp * ran2) * H(i))));
						Further[j] = Distr(Further_Swarm[j], a, b, Variance);

					}
					else {
						Further[j] = Distr(Pop[Ref_index][j], a, b, Variance);

					}
				}
			}
			else {
				Variance = 2 * Sigma(i + M);
				Further[j] = Distr(Pop[Ref_index][j], a, b, Variance);
			}
		}
		Fitness_candidate = J(Further);
		numb_iteration++;
		if (!(i % 7000000)) {
			std::vector<std::vector<long double>> Table = Pop.PopulationToMatrix();
			saveSetVectorsToCSV(Table, "C:\\Users\\troit\\OneDrive\\Рабочий стол\\Science_result\\SoFASM_Pop.csv", "SoFASM_Pop");
		}
		if (Fitness_candidate < Fitness[Ref_index]) {
			Fitness[Ref_index] = Fitness_candidate;
			Pop[Ref_index] = Further;

			exp_Fitness[Ref_index] = exponential_Fitness(Fitness_candidate, dim);
			if (Fitness_candidate < Fitness[min_index]) {
				min_index = Ref_index;

			}
		}

		X.push_back(numb_iteration);
		Y.push_back(Fitness[min_index]);
	}
	return Y;
}

void MiddleSSM_without_H_parallel(int N, int T, long double (*J)(Point), double (*Phi)(size_t),
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
		Point F = SoFASM_without_H(N, T, J, Phi, Sigma, H, dim, a, b, Q, M, current_iter, Distr);

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
	std::cout << "SoFASM without H " << exponential_Fitness(middle1[T - 1], dim);
	std::cout << std::endl;
	plt::named_plot(plot_name, iterations[NUM_RUNS - 1], middle1.getPoint());
	plt::legend();
}
std::vector<long double> SoFASM_without_rand(int N, int T, long double (*J)(Point), double (*Phi)(size_t), double (*Sigma)(size_t), double (*H)(size_t), int dim, double a, double b, double Q, size_t M, std::vector<int>& X, long double(*Distr)(long double, double, double, double)) {
	double Variance = Sigma(0);
	size_t m = 0;
	//double P;
	std::vector<long double> ref_numb;
	int numb_iteration = 0;
	double Threshold = 0.25;
	double Thresold2 = 0.6;
	double C1 = 0.5;
	double C2 = .05;
	double C3 = 2;
	double Temp;
	double tt;
	double A;
	int flag_index;
	int Flag;
	double ran2;
	double L;
	int help;
	std::vector<long double> exp_Fitness;
	long double Fitness_candidate;
	long double Fitness_min = INT_MAX;
	std::vector<long double> Fitness;
	std::vector<int> vec_flag = { 1, -1 };
	int rand_leader_index;
	int min_index = 0;
	//std::vector<std::vector<double>> Population=CreatePopulation(N,dim,a,b);
	Population Pop = Population(N, dim);
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

	for (int i = 0; i < T; i++) {
		if (i == 50000)
			help = 390;
		Temp = std::exp(-i * 1.0 / T);
		tt = (T - i) * 1.0 / T;  // Будет от 1.0 до 0.0
		L = C1 * std::exp(-tt);
		if (L > 1)
			L = 1;
		Ref_index = ReferencePoint(Phi, exp_Fitness, i);
		Point Further_Swarm(dim);

		//std::cout << m << std::endl;
		ref_numb.push_back(Ref_index);


		Point Further(dim);
		if (i == 1000) {
			int r = 3;
		}

		/*for (size_t j = 0; j < Ref.getDimension(); j++) {
			double eta = my_rand() % 1001 * 1.0 / 1000;
			double Expr = eta * (boost::math::erf((b - Ref[j]) / sqrt(2 * Variance))) + (1 - eta) * boost::math::erf((a - Ref[j]) / sqrt(2 * Variance));
			double R;
			if (Expr >= 1) {
				R = b;
			}
			else if (Expr <= -1) {
				R = a;
			}
			else {
				R = sqrt(2 * Variance) * (boost::math::erf_inv(Expr)) + Ref[j];
			}
			//double R = Mutation(Ref[j], a, b, Variance);
			Further[j]=R;
		}*/
		for (int j = 0; j < dim; ++j) {
			double q = my_random_d(0, 1);
			if (q < Q) {
				Variance = 2 * Sigma(i);
				if (L < Threshold) {
					rand_leader_index = my_random_i(N);
					double denominator = std::max(Fitness[Ref_index], 1e-300L);
					A = std::exp(-Fitness[rand_leader_index] / denominator);

					flag_index = (my_random_i(2));
					Flag = vec_flag[flag_index];
					ran2 = my_random_d(0, 1);
					Further_Swarm[j] = std::min(b * 1.0L, (std::max(a * 1.0L, Pop[Ref_index][j] + H(i) * (Flag * C2 * A * ((b - a) + a)))));
					Further[j] = Distr(Further_Swarm[j], a, b, Variance);
				}
				else {
					if (Temp > Thresold2) {
						rand_leader_index = my_random_i(N);
						flag_index = my_random_i(2);
						Flag = vec_flag[flag_index];
						ran2 = my_random_d();
						Further_Swarm[j] = std::min(b * 1.0L, (std::max(a * 1.0L, Pop[min_index][j] + ((Pop[min_index][j] - Pop[Ref_index][j]) * C3 * Flag * Temp * ran2) * H(i))));
						Further[j] = Distr(Further_Swarm[j], a, b, Variance);

					}
					else {
						Further[j] = Distr(Pop[Ref_index][j], a, b, Variance);

					}
				}
			}
			else {
				Variance = 2 * Sigma(i + M);
				Further[j] = Distr(Pop[Ref_index][j], a, b, Variance);
			}
		}
		Fitness_candidate = J(Further);
		numb_iteration++;
		if (!(i % 7000000)) {
			std::vector<std::vector<long double>> Table = Pop.PopulationToMatrix();
			saveSetVectorsToCSV(Table, "C:\\Users\\troit\\OneDrive\\Рабочий стол\\Science_result\\SoFASM_Pop.csv", "SoFASM_Pop");
		}
		if (Fitness_candidate < Fitness[Ref_index]) {
			Fitness[Ref_index] = Fitness_candidate;
			Pop[Ref_index] = Further;

			exp_Fitness[Ref_index] = exponential_Fitness(Fitness_candidate, dim);
			if (Fitness_candidate < Fitness[min_index]) {
				min_index = Ref_index;

			}
		}

		X.push_back(numb_iteration);
		Y.push_back(Fitness[min_index]);
	}
	return Y;
}

void MiddleSSM_without_Rand_parallel(int N, int T, long double (*J)(Point), double (*Phi)(size_t),
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
		Point F = SoFASM_without_rand(N, T, J, Phi, Sigma, H, dim, a, b, Q, M, current_iter, Distr);

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
	std::cout << std::endl;
	std::cout << "SoFASM without Rand " << exponential_Fitness(middle1[T - 1], dim);
	std::cout << std::endl;
	// Для графика используем итерации последнего запуска
	// (можно также усреднить, но это сложнее)
	plt::named_plot(plot_name, iterations[NUM_RUNS - 1], middle1.getPoint());
	plt::legend();
}

/*std::vector<long double> SoFAM_DE(int N, int T, Point(*J)(Point, double, bool), double (*Phi)(size_t), double (*Sigma)(size_t), int dim, Point a, Point b, double Q, size_t M, std::vector<int>& X, long double (*Distr)(long double, double, double, double), const char* MutName) {

	double Variance = Sigma(0);
	std::vector<double> ref_numb;
	int numb_iteration = 0;
	//double Fitness_candidate;
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

	std::vector<long double> Fitness_candidate(N);
	Point Further1(dim);
	std::vector<Point> Further(N, Further1);
	for (int l = 0; l < N; ++l) {
		for (int i = 0; i < dim; ++i) {
			Further[l][i] = my_random_d(a, b);

		}
		Fitness_candidate[l] = J(Further[l]);
	}
	for (size_t i = 0; i < T; i++) {
		double F_i = 0.5;
		std::vector<long double> poweredF(N);
		//std::cout << m << std::endl;
		//ref_numb.push_back(Ref_index);
		//std::vector<long double> exp_Fitness_Copy = exp_Fitness;
		Point Further_without_mutation(dim);
		double q;
		//Variance = 2 * Sigma(i);
		long double sum_prob = Give_Sum_Probability(Phi, exp_Fitness, i, poweredF);
		//long double sum_prob2 = sum_prob;
		std::vector<long double> poweredF_Copy;

		sum_prob -= poweredF[min_index];
		long double sum_prob2 = sum_prob;
		for (int k = 0; k < N - 4; ++k) {
			poweredF=Zero_Copy(poweredF, min_index);
			poweredF_Copy = poweredF;
			int Ref_index = ReferencePoint_New(sum_prob, poweredF_Copy);
			sum_prob -= poweredF[Ref_index];
			poweredF_Copy = Zero_Copy(poweredF_Copy, Ref_index);
			poweredF = Zero_Copy(poweredF, Ref_index);
			//exp_Fitness = Zero_Copy(exp_Fitness, Ref_index);

			sum_prob2 -= poweredF[Ref_index];
			int  Ref_index1 = ReferencePoint_New(sum_prob2, poweredF_Copy);
			int jr = my_random_i(0, N);
			//std::vector<long double> New2 = Zero_Copy(New1, Ref_index1);
			for (int j = 0; j < dim; j++) {


				double UpdateF = my_random_d();
				if (UpdateF <= 0.1) {
					F_i = my_random_d();
				}
				q = my_random_d();
				double UpdateCR = my_random_d();
				if (UpdateCR <= 0.1) {
					Q = my_random_d();
				}
				if (q < Q) {
					Variance = 2 * Sigma(i);
					Further_without_mutation[j] = Pop[min_index][j] + F_i * (Pop[Ref_index][j] - Pop[Ref_index1][j]);
					Further[k][j] = Distr(Further_without_mutation[j], a, b, Variance);
				}
				else {
					Variance = 2 * Sigma(i + M);
					Further[k][j] = Distr(Pop[Ref_index][j], a, b, Variance);
				}
			}


		}
		//std::vector<long double> Fitness_candidate(N);
		for (int l = 0; l < N; ++l) {
			long double Jf = J(Further[l]);
			Fitness_candidate[l] = Jf;

			if (Fitness_candidate[l] < Fitness[l]) {
				Pop[l] = Further[l];
				Fitness[l] = Fitness_candidate[l];
				exp_Fitness[l] = exponential_Fitness(Fitness_candidate[l], dim);
				if (Fitness_candidate[l] < Fitness_min) {
					Fitness_min = Fitness_candidate[l];
					min_index = l;
				}
			}
		}
		numb_iteration += N;
		X.push_back(numb_iteration);
		//Y.push_back(exponential_Fitness(Fitness_min,dim));
		Y.push_back(Fitness_min);
	}*/
	/*for (int i = 0; i < dim; ++i) {
		std::cout<<Pop[min_index][i]<<" ";
	}
	std::cout << "\n";*/
	/*plt::plot(ref_numb, X, ".");
	plt::xlabel("Reference point index");
	plt::ylabel("Iteration");
	plt::show();
	return Y;
}
void MiddleSM_DE(int N, int T, long double (*J)(Point), double (*Phi)(size_t), double (*Sigma)(size_t),
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
		else {
		F = SoFAM_DE(N, T, J, Phi, Sigma, dim, a, b, Q, M, Iteration, Distr, "DE_best1");
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
}*/

std::vector<long double> SoFAM_DE2(int N, int T, long double (*J)(Point), double (*Phi)(size_t), double (*Sigma)(size_t), int dim, double a, double b, double Q, size_t M, std::vector<int>& X, long double (*Distr)(long double, double, double, double), const char* MutName) {
	double Variance = Sigma(0);
	std::vector<double> ref_numb;
	int numb_iteration = 0;
	//double Fitness_candidate;
	std::vector<long double> exp_Fitness;
	std::vector<long double> Fitness;
	Population Pop = Population(N, dim);
	long double Fitness_min = INT_MAX;
	double CR = 0.9;
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
		double F_i = 0.5;
		std::vector<long double> poweredF(N);
		//std::cout << m << std::endl;
		//ref_numb.push_back(Ref_index);
		//std::vector<long double> exp_Fitness_Copy = exp_Fitness;
		Point Further_without_mutation(dim);
		double q;
		//Variance = 2 * Sigma(i);
		Ref_index = ReferencePoint(Phi, exp_Fitness, i);

		Point Further(dim);
		std::vector<Point> Further_Pop(N, Further);

		for (int j = 0; j < dim; ++j) {

			double q = my_random_d(0, 1);
			if (q < Q) {
				Variance = 2 * Sigma(i);
				Further_Pop[Ref_index][j] = Distr(Pop[Ref_index][j], a, b, Variance);
			}
			else {
				Variance = 2 * Sigma(i + M);
				Further_Pop[Ref_index][j] = Distr(Pop[Ref_index][j], a, b, Variance);

			}

		}
		for (int l = 0; l < N; ++l) {
			if (l == Ref_index) {
				continue;
			}
			Point Mut = Mutation(Pop, MutName, F_i, a, b, l, min_index);

			Further_Pop[l] = Crossingover(Mut, Pop[l], CR);
		}
		for (int l = 0; l < N; ++l) {
			long double Fitness_Cand = J(Further_Pop[l]);
			if (Fitness_Cand < Fitness[l]) {
				Fitness[l] = Fitness_Cand;
				Pop[l] = Further_Pop[l];
				exp_Fitness[l] = exponential_Fitness(Fitness[l], dim);
				if (Fitness_Cand < Fitness[min_index]) {
					min_index = l;
				}
			}
		}
		std::cout << i << std::endl;
		numb_iteration += N;
		X.push_back(numb_iteration);
		//Y.push_back(exponential_Fitness(Fitness_min,dim));
		Y.push_back(Fitness[min_index]);

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
void MiddleSM_DE2(int N, int T, long double (*J)(Point), double (*Phi)(size_t), double (*Sigma)(size_t),
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
			F = SoFAM_DE2(N, T, J, Phi, Sigma, dim, a, b, Q, M, Iteration, Distr, "DE_best1");
			//}
			for (int j = 0; j < F.getDimension(); ++j) {
				F[j] = log(F[j] + 1);
			}
			middle = middle + F;

			Table.push_back(F);
		}
	}
	saveSetVectorsToCSV(Table, "C:\\Users\\troit\\OneDrive\\Рабочий стол\\Science_result\\SoFA_Var_Mod.csv", "SoFA_Var_Mod");
	Point middle1 = middle * 0.1;
	plt::named_plot(plot_name, Iteration, middle1.getPoint());
	plt::legend();
	Iteration.clear();

}

std::vector<long double> SoFDE2(int N, int T, long double (*J)(Point), double (*Phi)(size_t), double (*Sigma)(size_t), double (*H)(size_t), int dim, double a, double b, double Q, size_t M, std::vector<int>& X, long double (*Distr)(long double, double, double, double), const char* MutName) {
	double Variance = Sigma(0);
	std::vector<double> ref_numb;
	int numb_iteration = 0;
	//double Fitness_candidate;
	std::vector<long double> exp_Fitness;
	std::vector<long double> Fitness;
	Population Pop = Population(N, dim);
	long double Fitness_min = INT_MAX;
	double CR = 0.9;
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
		double F_i = 0.5;
		std::vector<long double> poweredF(N);
		//std::cout << m << std::endl;
		//ref_numb.push_back(Ref_index);
		//std::vector<long double> exp_Fitness_Copy = exp_Fitness;
		Point Further_without_mutation(dim);
		double q;
		//Variance = 2 * Sigma(i);
		Ref_index = ReferencePoint(Phi, exp_Fitness, i);

		Point Further(dim);
		std::vector<Point> Further_Pop(N, Further);

		for (int j = 0; j < dim; ++j) {

			double q = my_random_d(0, 1);
			if (q < Q) {
				Variance = 2 * Sigma(i);
				Further_Pop[Ref_index][j] = Distr(Pop[Ref_index][j], a, b, Variance);
			}
			else {
				Variance = 2 * Sigma(i + M);
				Further_Pop[Ref_index][j] = Distr(Pop[Ref_index][j], a, b, Variance);

			}

		}
		double s = my_random_d();
		if (H(i) > s) {
			for (int l = 0; l < N; ++l) {
				if (l == Ref_index) {
					continue;
				}
				Point Mut = Mutation(Pop, MutName, F_i, a, b, l, min_index);

				Further_Pop[l] = Crossingover(Mut, Pop[l], CR);
			}
			for (int l = 0; l < N; ++l) {
				long double Fitness_Cand = J(Further_Pop[l]);
				if (Fitness_Cand < Fitness[l]) {
					Fitness[l] = Fitness_Cand;
					Pop[l] = Further_Pop[l];
					exp_Fitness[l] = exponential_Fitness(Fitness[l], dim);
					if (Fitness_Cand < Fitness[min_index]) {
						min_index = l;
					}
				}
			}
			numb_iteration += N;
		}
		else {
			long double FCand = J(Further_Pop[Ref_index]);
			if (Fitness[Ref_index] > FCand) {
				Fitness[Ref_index] = FCand;
				Pop[Ref_index] = Further_Pop[Ref_index];
				exp_Fitness[Ref_index] = exponential_Fitness(Fitness[Ref_index], dim);
				if (FCand < Fitness[min_index]) {
					min_index = Ref_index;
				}
			}
			numb_iteration++;
		}
		//std::cout << i << std::endl;
		X.push_back(numb_iteration);
		//Y.push_back(exponential_Fitness(Fitness_min,dim));
		Y.push_back(Fitness[min_index]);

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

void MiddleSF_DE2(int N, int T, long double (*J)(Point), double (*Phi)(size_t), double (*Sigma)(size_t), double (*H)(size_t),
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
			F = SoFDE2(N, T, J, Phi, Sigma, H, dim, a, b, Q, M, Iteration, Distr, "DE_best1");
			//}
			for (int j = 0; j < F.getDimension(); ++j) {
				F[j] = log(F[j] + 1);
			}
			middle = middle + F;

			Table.push_back(F);
		}
	}
	saveSetVectorsToCSV(Table, "C:\\Users\\troit\\OneDrive\\Рабочий стол\\Science_result\\SoFA_Var_Mod.csv", "SoFA_Var_Mod");
	Point middle1 = middle * 0.1;
	plt::named_plot(plot_name, Iteration, middle1.getPoint());
	plt::legend();
	Iteration.clear();

}

std::vector<long double> SoFDE3(int N, int T, long double (*J)(Point), double (*Phi)(size_t), double (*Sigma)(size_t), double (*H)(size_t), int dim, double a, double b, double Q, size_t M, std::vector<int>& X, long double (*Distr)(long double, double, double, double), const char* MutName) {
	double Variance = Sigma(0);
	std::vector<double> ref_numb;
	int numb_iteration = 0;
	//double Fitness_candidate;
	std::vector<long double> exp_Fitness;
	std::vector<long double> Fitness;
	Population Pop = Population(N, dim);
	long double Fitness_min = INT_MAX;
	double CR = 0.9;
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
		double F_i = 0.5;
		std::vector<long double> poweredF(N);
		//std::cout << m << std::endl;
		//ref_numb.push_back(Ref_index);
		//std::vector<long double> exp_Fitness_Copy = exp_Fitness;
		Point Further_without_mutation(dim);
		double q;
		//Variance = 2 * Sigma(i);
		Ref_index = ReferencePoint(Phi, exp_Fitness, i);

		Point Further(dim);
		std::vector<Point> Further_Pop(N, Further);

		for (int j = 0; j < dim; ++j) {

			double q = my_random_d(0, 1);
			if (q < Q) {
				Variance = 2 * Sigma(i);
				Further_Pop[Ref_index][j] = Distr(Pop[Ref_index][j], a, b, Variance);
			}
			else {
				Variance = 2 * Sigma(i + M);
				Further_Pop[Ref_index][j] = Distr(Pop[Ref_index][j], a, b, Variance);

			}

		}
		if (i == 10000) {
			int z = 78;
		}
		long double h = H(i);
		numb_iteration++;
		bool f = true;
		for (int l = 0; l < N; ++l) {
			if (l == Ref_index) {
				continue;
			}
			double s = my_random_d();
			

			if (h > s) {
				/*if (f) {
					std::cout << i << std::endl;
				}*/
				numb_iteration++;
				Point Mut = Mutation(Pop, MutName, F_i, a, b, l, min_index);
				Further_Pop[l] = Crossingover(Mut, Pop[l], CR);
				f = false;
			}
			else {
				Further_Pop[l] = Pop[l];
			}

		}
		for (int l = 0; l < N; ++l) {
			long double Fitness_Cand = J(Further_Pop[l]);
			if (Fitness_Cand < Fitness[l]) {
				Fitness[l] = Fitness_Cand;
				Pop[l] = Further_Pop[l];
				exp_Fitness[l] = exponential_Fitness(Fitness[l], dim);
				if (Fitness_Cand < Fitness[min_index]) {
					min_index = l;
				}
			}
		}





		//std::cout << i << std::endl;
		X.push_back(numb_iteration);
		//Y.push_back(exponential_Fitness(Fitness_min,dim));
		Y.push_back(Fitness[min_index]);
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

void MiddleSF_DE3(int N, int T, long double (*J)(Point), double (*Phi)(size_t), double (*Sigma)(size_t), double (*H)(size_t),
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
			F = SoFDE3(N, T, J, Phi, Sigma, H, dim, a, b, Q, M, Iteration, Distr, "DE_best1");
			//}
			for (int j = 0; j < F.getDimension(); ++j) {
				F[j] = log(F[j] + 1);
			}
			middle = middle + F;

			Table.push_back(F);
		}
	}
	saveSetVectorsToCSV(Table, "C:\\Users\\troit\\OneDrive\\Рабочий стол\\Science_result\\SoFA_Var_Mod.csv", "SoFA_Var_Mod");
	Point middle1 = middle * 0.1;
	plt::named_plot(plot_name, Iteration, middle1.getPoint());
	plt::legend();
	Iteration.clear();

}

std::vector<long double> SoFDE_Grad(int N, int T, long double (*J)(Point), double (*Phi)(size_t), double (*Sigma)(size_t), double (*H)(size_t), int dim, double a, double b, double Q, size_t M, std::vector<int>& X, long double (*Distr)(long double, double, double, double), const char* MutName) {
	double Variance = Sigma(0);
	std::vector<double> ref_numb;
	int numb_iteration = 0;
	//double Fitness_candidate;
	std::vector<long double> exp_Fitness;
	std::vector<long double> Fitness;
	Population Pop = Population(N, dim);
	long double Fitness_min = INT_MAX;
	double CR = 0.9;
	int min_index = 0;
	numb_iteration += N;
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
	for (size_t i = 0; i < T; i++) {
		double F_i = 0.5;
		std::vector<long double> poweredF(N);
		//std::cout << m << std::endl;
		//ref_numb.push_back(Ref_index);
		//std::vector<long double> exp_Fitness_Copy = exp_Fitness;
		Point Further_without_mutation(dim);
		double q;
		//Variance = 2 * Sigma(i);
		Ref_index = ReferencePoint(Phi, exp_Fitness, i);

		Point Further(dim);
		std::vector<Point> Further_Pop(N, Further);
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
		for (int j = 0; j < dim; ++j) {

			double q = my_random_d(0, 1);
			if (q < Q) {
				Variance = 2 * Sigma(i) * Ref_Grad[j];
				Further_Pop[Ref_index][j] = Distr(Pop[Ref_index][j], a, b, Variance);
			}
			else {
				Variance = 2 * Sigma(i + M) * Ref_Grad[j];
				Further_Pop[Ref_index][j] = Distr(Pop[Ref_index][j], a, b, Variance);

			}

		}
		numb_iteration++;
		bool f = true;
		for (int l = 0; l < N; ++l) {
			if (l == Ref_index) {
				continue;
			}
			double s = my_random_d();
			long double h = H(i);

			if (h > s) {
				if (f) {
				}
				numb_iteration++;
				Point Mut = Mutation(Pop, MutName, F_i, a, b, l, min_index);
				Further_Pop[l] = Crossingover(Mut, Pop[l], CR);
				f = false;
			}
			else {
				Further_Pop[l] = Pop[l];
			}

		}
		for (int l = 0; l < N; ++l) {
			long double Fitness_Cand = J(Further_Pop[l]);
			if (Fitness_Cand < Fitness[l]) {
				Fitness[l] = Fitness_Cand;
				Pop[l] = Further_Pop[l];
				exp_Fitness[l] = exponential_Fitness(Fitness[l], dim);
				if (Fitness_Cand < Fitness[min_index]) {
					min_index = l;
				}
			}
		}





		//std::cout << i << std::endl;
		X.push_back(numb_iteration);
		//Y.push_back(exponential_Fitness(Fitness_min,dim));
		Y.push_back(Fitness[min_index]);
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

void MiddleSFDE_Grad(int N, int T, long double (*J)(Point), double (*Phi)(size_t), double (*Sigma)(size_t), double (*H)(size_t),
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
			F = SoFDE_Grad(N, T, J, Phi, Sigma, H, dim, a, b, Q, M, Iteration, Distr, "DE_best1");
			//}
			for (int j = 0; j < F.getDimension(); ++j) {
				F[j] = log(F[j] + 1);
			}
			middle = middle + F;

			Table.push_back(F);
		}
	}
	saveSetVectorsToCSV(Table, "C:\\Users\\troit\\OneDrive\\Рабочий стол\\Science_result\\SoFA_Var_Mod.csv", "SoFA_Var_Mod");
	Point middle1 = middle * 0.1;
	plt::named_plot(plot_name, Iteration, middle1.getPoint());
	plt::legend();
	Iteration.clear();

}