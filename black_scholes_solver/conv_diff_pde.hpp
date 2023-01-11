#pragma once
#include "option.hpp"
class conv_diff_pde {
public:
	virtual double trans_term(double t, double S) = 0;
	virtual double conv_term(double t, double S) = 0;
	virtual double diff_term(double t, double S) = 0;
	virtual double source_term(double t, double S) = 0;

	virtual double bdry_left_cond(double t, double S) = 0;
	virtual double bdry_right_cond(double t, double S) = 0;
	virtual double init_cond(double S) = 0;
};

class blackscholes_pde : public conv_diff_pde {
public:
	blackscholes_pde(VanillaOption* _option) : option(_option) {}
	
	// BS-model specific coefficients
	double trans_term(double t, double S) { return option->r; };
	double conv_term(double t, double S) { return option->r * S; };
	double diff_term(double t, double S) { return 0.5 * option->sigma * option->sigma * S * S; };
	double source_term(double t, double S) { return 0.0; };

	//Option-specific boundary + initial (terminal) conditions
	double bdry_left_cond(double t, double S) { return option->left_bdry_value(t, S); };
	double bdry_right_cond(double t, double S) { return option->right_bdry_value(t, S); };
	double init_cond(double S) { return option->payoff(S); };

	VanillaOption* option;
};