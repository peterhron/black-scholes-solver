#include<algorithm>

class payoff{
public: 
	payoff() {};
	virtual ~payoff() {};

	virtual double operator() (const double& S) const = 0;
};

class payoff_call : public payoff {
public:
	payoff_call(double _K) { K = _K; }
	virtual ~payoff_call() {};

	double K;
	virtual double operator() (const double& S) const { return std::max(S - K, 0.0); }
};

class payoff_put : public payoff {
public:
	payoff_put(double _K) {K = _K;}
	virtual ~payoff_put() {};

	double K;
	virtual double operator() (const double& S) const { return std::max(K - S, 0.0); }
};

