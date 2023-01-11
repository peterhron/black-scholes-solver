#include "FDM.hpp"

int main()
{
	double K = 0.5;  // strike
	double r = 0.05;   // risk-free rate
	double v = 0.2;    // volatility of the U/L
	double T = 1.00;    // Time to maturity

	// FDM discretisation parameters
	double x_dom = 1.0;       // Spot price in [0.0, 1.0]
	unsigned long I = 20;     // #of discretized price steps
	double t_dom = T;         // time in [0.0, T]
	unsigned long N = 20;     // #of discretized time steps

	VanillaPut* opt = new VanillaPut(K, T, r, v);

	blackscholes_pde* bs_eqn = new blackscholes_pde(opt);
	FDMEulerExplicit soln_obj(x_dom, I, t_dom, N, bs_eqn);
	soln_obj.time_march();

	delete bs_eqn;
	delete opt;
	return 7;
}