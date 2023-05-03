function x_dot = f(x)
mu = 1/(365*76);
omega = 1/365;
gamma = 1/14;
sigma = 1/7;
alpha = 0;
beta = 0.21;

S = x(1);
E = x(2);
I = x(3);
R = x(4);
N = x(1)+x(2)+x(3)+x(4);
S_dot = mu*N-beta*I*S/N+omega*R-mu*S;
E_dot = beta*I*S/N-sigma*E-mu*E;
I_dot = sigma*E-gamma*I-(mu+alpha)*I;
R_dot = gamma*I-omega*R-mu*R;

x_dot = [S_dot;E_dot;I_dot;R_dot];
end
