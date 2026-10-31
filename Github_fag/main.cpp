#include "matplotlibcpp.h"
#include <math.h>
//#include "Convergence.h"
#include "NewMethods.h"
#include "SnakeOptimizer.h"

#include "DiffEval.h"
int main() {
	
	setlocale(LC_ALL, "Russian");
	
	
	MiddleSFDE_Grad(100, 1000000, Rasstrigin, Phi1, Sigma11, H2, 10, -5.12, 5.12, 0.5, 1000, Cauchy, "SoFDEGrad");

	//MiddleSFDE_Grad(100, 1000000, Ackley, Phi, Sigma11, H, 10, -33, 33, 0.5, 1000, Cauchy, "SoFDE_Grad");
	//MiddleSF(100, 1000000, Ackley, Phi, Sigma11, 10, -33, 33, 0.5, 1000, Cauchy, "SoFA");
	
	std::vector<long double> X;

	
	bool B = false;
	MiddleSF(100, 1000000, Rasstrigin, Phi1, Sigma11, 10, -5.12, 5.12, 0.5, 1000, Cauchy, "SoFA", B);
	Middle(100, 10000, Rasstrigin, 10, -5.12, 5.12, "SO");
	MiddleDE(Rasstrigin, 10000, 10, 100, -5.12, 5.12, "DE_Rand1");
	MiddleDE(Rasstrigin, 10000, 10, 100, -5.12, 5.12, "DE_Best1");
	//MiddleSF_parallel(10, 10000000, Rasstrigin, Phi1, Sigma2, 10, -5, 5, 0.5, 1000, Cauchy, "SoFA_parallel", B);
	//MiddleSF_Var(11, 200000, Styblinsky_Tang, Phi1, Sigma2, 10, -5, 5,  0.5, 1000, Cauchy, "SoFA_Var_Mod", B);
	//MiddleSFVar_parallel(10, 10000000, Rasstrigin, Phi1, Sigma2, 10, -5, 5, 0.5, 1000, Cauchy, "SoFA_Var_Mod_parallel", B);
	//MiddleSF_Var_Numb(11, 2000000, Styblinsky_Tang, Phi1, Sigma2, 10, -5, 5, 0.5, 1000, Cauchy, "SoFA_Var_Mod_Numb",B);
	//Middle(100, 20000, Rasstrigin, 10, -5, 5, "SO", B);
	 //MiddleSM_Grad(11, 200000, Styblinsky_Tang, Phi1, Sigma2, 10, -5, 5, 0.5, 1000, Cauchy, "SoFAM_Grad", B);
	 
	// MiddleSM_Grad2(11, 200000, Styblinsky_Tang, Phi1, Sigma2, 10, -5, 5, 0.5, 1000, Cauchy, "SoFAM_Grad", B);
	 //MiddleSoFAM_Grad(11, 200000, Rasstrigin, Phi1, Sigma2, 10, -5, 5, 0.5, 1000, Cauchy, "SoFAM_Grad", B);
	//plt::title("Rasstrigin");
	//MiddleSSM(11, 2000000, Rasstrigin, Phi1, Sigma2, H, 10, -5, 5, 0.5, 1000, Cauchy, "SoFASM", B);
	//MiddleSF_Var_Numb_parallel(11, 1000000, Rasstrigin, Phi1, Sigma2, 10, -5, 5, 0.5, 1000, Cauchy, "SoFA_Var_Mod_parallel", B);
	//MiddleSSM_parallel(12, 1000000, Rasstrigin, Phi1, Sigma2,H, 10, -5, 5, 0.5, 1000, Cauchy, "SoFASM");
	//MiddleSF_parallel(3, 10000000, Rasstrigin, Phi1, Sigma2, 10, -5, 5, 0.5, 1000, Cauchy, "SoFA");
	//MiddleSF_parallel(3, 10000000, Rasstrigin, Phi1, Sigma2, 10, -5, 5, 0.5, 1000, Cauchy, "SoFA");
	//MiddleSSM_without_H_and_Rand_parallel(10, 10000000, Schwefel, Phi1, Sigma2, 10, -500, 500, 0.5, 1000, Cauchy, "SoFASM_without_H_and_rand");
	//MiddleSFVar_parallel(10, 100000, Styblinsky_Tang, Phi1, Sigma2, 10, -5, 5, 0.5, 1000, Cauchy, "SoFA_Var_mod");
	//MiddleSM_parallel(10, 1000000, Styblinsky_Tang, Phi1, Sigma2,H, 10, -5, 5, 0.5, 1000, Cauchy, "SoFAM");
	//Middle_parallel(10, 100000, Styblinsky_Tang, 10, -5, 5, "SO");

//>>>>>>> 968935c38924bbe3f0eab988a9b2ba7c6e3a74c4
	//MiddleSF_parallel(6, 5000000, Ackley, Phi1, Sigma2, 10, -5, 5, 0.5, 1000, Cauchy, "SoFA");
	/*MiddleSSM_without_H_and_Rand_parallel(6, 5000000, Styblinsky_Tang, Phi1, Sigma2, 10, -5, 5, 0.2, 1000, Cauchy, "SoFASM_0.2");
	MiddleSSM_without_H_and_Rand_parallel(6, 5000000, Styblinsky_Tang, Phi1, Sigma2, 10, -5, 5, 0.5, 1000, Cauchy, "SoFASM_0.5");
	MiddleSSM_without_H_and_Rand_parallel(6, 5000000, Styblinsky_Tang, Phi1, Sigma2, 10, -5, 5, 0.7, 1000, Cauchy, "SoFASM_0.7");
	MiddleSSM_without_H_and_Rand_parallel(6, 5000000, Styblinsky_Tang, Phi1, Sigma2, 10, -5, 5, 0.9, 1000, Cauchy, "SoFASM_0.9");*/
	//MiddleSSM(10, 10000000, Styblinsky_Tang, Phi1, Sigma5, H, 10, -5, 5, 0.1, 1000, Cauchy, "Sigma5");
	/*MiddleSSM(6, 4000000, Styblinsky_Tang, Phi1, Sigma2, H, 10, -5, 5, 0.2, 5000, Cauchy, "M=5000 Q=0.2");
	MiddleSSM(6, 4000000, Styblinsky_Tang, Phi1, Sigma2, H, 10, -5, 5, 0.4, 1000, Cauchy, "M=1000 Q=0.4");
	MiddleSSM(6, 4000000, Styblinsky_Tang, Phi1, Sigma2, H, 10, -5, 5, 0.4, 5000, Cauchy, "M=5000 Q=0.4");
	MiddleSSM(6, 4000000, Styblinsky_Tang, Phi1, Sigma2, H, 10, -5, 5, 0.5, 1000, Cauchy, "M=1000 Q=0.5");
	MiddleSSM(6, 4000000, Styblinsky_Tang, Phi1, Sigma2, H, 10, -5, 5, 0.5, 5000, Cauchy, "M=5000 Q=0.5");
	MiddleSSM(6, 4000000, Styblinsky_Tang, Phi1, Sigma2, H, 10, -5, 5, 0.8, 1000, Cauchy, "M=1000 Q=0.8");
	MiddleSSM(6, 4000000, Styblinsky_Tang, Phi1, Sigma2, H, 10, -5, 5, 0.8, 5000, Cauchy, "M=5000 Q=0.8");*/
	/*MiddleSM(6, 3000000, Rasstrigin, Phi1, Sigma2, H, 10, -5, 5, 0.5, 1000, Cauchy, "SoFAM");
	MiddleSF(6, 3000000, Rasstrigin, Phi1, Sigma2, 10, -5, 5, 0.5, 1000, Cauchy, "SoFA");
	Middle_Gradient(600, 20, Rasstrigin, 10, -5, 5, "Gradient_Search");
	Middle(6, 500000, Rasstrigin,10 , -5, 5, "SO");*/
	/*MiddleSM(10, 200000, GKLS_x, Phi1, Sigma_SM, H, 10, -5, 5, 0.5, 1000, Cauchy, "Sequence1");
	MiddleSM(10, 200000, GKLS_x, Phi1, Sigma11, H, 10, -5, 5, 0.5, 1000, Cauchy, "Sequence2");
	MiddleSM(10, 200000, GKLS_x, Phi1, Sigma2, H, 10, -5, 5, 0.5, 1000, Cauchy, "Sequence3");
	MiddleSM(10, 200000, GKLS_x, Phi1, Sigma3, H, 10, -5, 5, 0.5, 1000, Cauchy, "Sequence4");
	MiddleSM(10, 200000, GKLS_x, Phi1, Sigma4, H, 10, -5, 5, 0.5, 1000, Cauchy, "Sequence5");
	MiddleSM(10, 200000, GKLS_x, Phi1, Sigma5, H, 10, -5, 5, 0.5, 1000, Cauchy, "Sequence6");*/


	/*MiddleSSM(6, 2000000, Rasstrigin, Phi, Sigma, H, 10, -5, 5, 0.3, 1000, Cauchy, "Phi Sigma");
	MiddleSSM(6, 2000000, Rasstrigin, Phi1, Sigma, H, 10, -5, 5, 0.3, 1000, Cauchy, "Phi1 Sigma");
	MiddleSSM(6, 2000000, Rasstrigin, Phi2, Sigma, H, 10, -5, 5, 0.3, 1000, Cauchy, "Phi2 Sigma");
	MiddleSO_Conv(100, 500, GKLS_x, 10, -10, 10, "SO");
	MiddleSF_Conv(100, 500, GKLS_x, Phi1, Sigma2, 10, -10, 10, 0.5, 1000, Cauchy, "SoFA");
	MiddleSFVar_Conv(100, 500, GKLS_x, Phi1, Sigma2, 10, -10, 10, 0.5, 1000, Cauchy, "SoFA_Var_Mod");
	MiddleSM_Grad_Conv(100, 500, GKLS_x, Phi1, Sigma2, 10, -10, 10, 0.5, 1000, Cauchy, "SoFAM_Grad");
	MiddleSF_Var_Numb_Conv(100, 500, GKLS_x, Phi1, Sigma2, 10, -10, 10, 0.5, 1000, Cauchy, "SoFA_Var_Mod_Numb");
	MiddleSM_Grad_Numb_Conv(100, 500, GKLS_x, Phi1, Sigma2, 10, -10, 10, 0.5, 1000, Cauchy, "SoFA_Grad_Numb");
	MiddleSM_Conv(100, 500, GKLS_x, Phi1, Sigma2, H, 10, -10, 10, 0.5, 1000, Cauchy, "SoFAM");*/
	//MiddleSSM_Conv(100, 500, GKLS_x, Phi1, Sigma2, 10, -10, 10, 0.5, 1000, Cauchy, "SoFASM");
	plt::grid(true);

	/*if (B) {
	plt::xlabel("Programm Working Time, ms");
	}
	else {*/
	//<<<<<<< HEAD
	plt::title("Rasstrigin");
	plt::xlabel("Number of calls Fitness-Function");

	plt::ylabel("Fitness");
	//=======
	//>>>>>>> 968935c38924bbe3f0eab988a9b2ba7c6e3a74c4
	plt::show();
	//plt::close();
	//MiddleSM(6, 200000, expRasstrigin, Phi, Sigma, H, 7, -5, 5, 0.5, 1000, Cauchy, "SoFAM");
	//Middle(6, 1000000, Rasstrigin, 7, -5, 5, "SO");
	//Middle(6, 100000, Rasstrigin, 7, -5, 5, "SO");




}