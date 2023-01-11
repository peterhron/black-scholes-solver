#pragma once
#include<vector>
#include<fstream>
#include"conv_diff_pde.hpp"

// Finite Difference Methods class
class FDMBase {
protected:
    conv_diff_pde* eqn;

    // Space discretisation
    double x_max;      // spatial domain [0.0, x_max]
    unsigned long I;   // Number of discretized spatial steps
    double dx = x_max / I;         // spatial step size
    std::vector<double> x_values;  // coordinates of the x dimension

    // Time discretisation
    double t_max;      // time domain [0.0, t_max]
    unsigned long N;   // Number of discretized time steps
    double dt;         // time step size (calculated from above)

    // Time steps
    double t_prev, t_curr;   // current and previous times

    // Differencing coefficients
    double alpha, beta, gamma;

    // Storage
    std::vector<double> new_result;   // new solution (becomes N+1)
    std::vector<double> old_result;   // old solution (becomes N)

    // Constructor
    FDMBase(double _x_max, int _I, double _t_max, int _N, conv_diff_pde* _eqn)
        : x_max(_x_max), I(_I), t_max(_t_max), N(_N), eqn(_eqn) {};

    // Override these virtual methods in derived classes for specific FDM techniques
    virtual void calculate_step_sizes() = 0;
    virtual void set_init_conds() = 0;
    virtual void calculate_bdry_conds() = 0;
    virtual void calculate_inner_domain() = 0;

public:
    // Carry out the time-steps
    virtual void time_march() = 0;
};


class FDMEulerExplicit : public FDMBase {
protected:
    void calculate_step_sizes() 
    {
        dx = x_max / static_cast<double>(I - 1);
        dt = t_max / static_cast<double>(N - 1);
    }

    void set_init_conds()
    {
        // spatial params
        double x_curr = 0.0;

        old_result.resize(I, 0.0);
        new_result.resize(I, 0.0);
        x_values.resize(I, 0.0);

        for (int i = 0; i < I; i++) 
        {
            x_curr = static_cast<double>(i) * dx;
            old_result[i] = eqn->init_cond(x_curr);
            x_values[i] = x_curr;
        }

        // temporal params
        t_prev = 0.0;
        t_curr = 0.0;
    }

    void calculate_bdry_conds()
    {
        new_result[0] = eqn->bdry_left_cond(t_prev, x_values[0]);
        new_result[I - 1] = eqn->bdry_right_cond(t_prev, x_values[I - 1]);
    }

    void calculate_inner_domain()
    {
        for (int i = 1; i < I - 1; i++) 
        {
            // differencing coefficients
            double dt_sig = dt * eqn->diff_term(t_prev, x_values[i]);
            double dt_sig_2 = dt * dx * 0.5 * eqn->conv_term(t_prev, x_values[i]);

            alpha = dt_sig - dt_sig_2;
            beta = dx * dx - 2.0 * dt_sig + dt * dx * dx * eqn->trans_term(t_prev, x_values[i]);
            gamma = dt_sig + dt_sig_2; 

            // calculate values in the grid
            new_result[i] = (alpha * old_result[i - 1] +
                beta * old_result[i] +
                gamma * old_result[i + 1]) / (dx * dx) -
                dt * eqn->source_term(t_prev, x_values[i]);
        }
    }

public:
    FDMEulerExplicit(double _x_max, unsigned long _I,
        double _t_max, unsigned long _N, conv_diff_pde* _eqn)
        : FDMBase(_x_max, _I, _t_max, _N, _eqn) 
    {
        calculate_step_sizes();
        set_init_conds();
    }

    void time_march()
    {
        std::ofstream soln_out("output_solution.csv");

        while (t_curr < t_max) 
        {
            t_curr = t_prev + dt;
            calculate_bdry_conds();
            calculate_inner_domain();
            for (int j = 0; j < I; j++) {
                soln_out << x_values[j] << " " << t_prev << " " << new_result[j] << std::endl;
            }

            old_result = new_result;
            t_prev = t_curr;
        }

        soln_out.close();
    }
};