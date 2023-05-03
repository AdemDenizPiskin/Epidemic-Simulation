function xNext = RungeKutta4Stages(x,h,f)
% =============================
% Function computes one step of Heun's method
% Inputs
% - x: current state
% - h: step size
% - f: right hand side of the IVP to be simulated (given as function
% handle): https://www.mathworks.com/help/matlab/matlab_prog/pass-a-function-to-another-function.html
% Output
% - xNext: state in the next time step

k1 = f(x);
k2 = f(x+0.5*h*k1);
k3 = f(x+0.5*h*k2);
k4 = f(x+1*h*k3);


xNext = x + h*(1/6*k1+2/6*k2+2/6*k3+1/6*k4);
