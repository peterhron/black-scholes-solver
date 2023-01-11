# black-scholes-solver
A simple solver for PDEs of the convection-diffusion type. 

To use, define the equation as a child of the *conv_diff_pde* class and apply 
an existing FDM method to it, or implement new methods under the *FDMBase* parent class.
Plan is to add new numerical methods as time goes, possibly for comparison with each other and/or with
analytical solutions where available (such as the Black-Scholes model)

Project is inspired by an article of https://www.quantstart.com/, 
and some parts of their code are used here.