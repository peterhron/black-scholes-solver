#pragma once
class VanillaOption {
public:
	VanillaOption() {};
	virtual ~VanillaOption() {};
	
	virtual double payoff(double S) = 0;
	virtual double left_bdry_value(double t, double S) = 0;
	virtual double right_bdry_value (double t, double S) = 0;
	double K = 0.0, T = 0.0, r = 0.0, sigma = 0.0;
};

class VanillaCall : public VanillaOption {
public:
	VanillaCall(double _K, double _T, double _r, double _sigma)
		: K(_K), T(_T), r(_r), sigma(_sigma) {}

	virtual double payoff(double S) { return std::max(S - K, 0.0); };
	virtual double left_bdry_value(double t, double S) { return 0.0; };
	virtual double right_bdry_value(double t, double S) { return S - K * exp(-r * (T - t)); };
	double K, T, r, sigma;
};

class VanillaPut : public VanillaOption {
public:
	VanillaPut(double _K, double _T, double _r, double _sigma)
		: K(_K), T(_T), r(_r), sigma(_sigma) {}

	virtual double payoff(double S) { return std::max(K - S, 0.0); };
	virtual double left_bdry_value(double t, double S) { return K * exp(-r * (T - t)); };
	virtual double right_bdry_value(double t, double S) { return 0; };
	double K, T, r, sigma;
};


