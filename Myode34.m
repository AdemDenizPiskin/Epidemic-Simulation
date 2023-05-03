% Some figure formatting
clear all;
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Heun's method for the IVP d/dt x = -x, x(0) = 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% right hand side to be simulated defined with the independent variable x
% d/dt x = -x -> f(x) = -x
%f = @(x)(-x); % you could use any other first-order ODE here
% Initial condition
x0 = [0.999; 0.001; 0; 0];
% =======================
% Simulation
% =======================
t_end = 1800; % final time
h0 = 10;% small step size; this is very close to the actual solution
h_min  = 0.01;
MAX_ITER = 3*ceil(t_end/h0);
h = zeros(1,MAX_ITER);
x(:,1) = x0; % state values during the simulation
time(1) = 0; % simulation time
h(1) = h0;
atol = 10e-5;
rtol = 10e-5;
fac0 = 0.2;
fac1 = 5;
beta  = 0.9;
time(1,1) = 0;
%x4(:,1) = x0;
%sigma_arr = 0;
%eta_arr = zeros(4,1);
tic
for kk = 1:MAX_ITER
    if time(end)+h(kk)>t_end
        h(kk) = t_end-time(end);
    end
    x_4Stages = RungeKutta4Stages(x(:,end),h(kk),@f);
    x_3Stages = RungeKutta3Stages(x(:,end),h(kk),@f);
    eta = abs(x_3Stages-x_4Stages);
    %eta_arr(:,kk) =eta;
    sigma = (sqrt(1/4*sum((eta./(atol+rtol*abs(x_4Stages))).^2)));%^(-1/4);
    %sigma_arr(kk) = sigma;
    h(kk+1) = h(kk)*min(fac1,max(fac0,beta*sigma^(-1/4)));
    %x4 = [x4 x_4Stages];
    if sigma<=1 
        time = [time (time(end)+h(kk))];
        x = [x x_4Stages];
    end
    if time(end)>=t_end
            break
    end
end
tt = toc;
figure;
plot(time/365,x(1,:));
hold on
plot(time/365,x(2,:));
hold on
plot(time/365,x(3,:));
hold on
plot(time/365,x(4,:));
xlabel('time [years]')
ylabel('Population percentage normalized')
grid on;
legend({'S','E','I','R'},'FontSize',14)
title('Epidemic Simulation with constant step size h=15 days')
set(findall(gcf,'Type','line'),'LineWidth',2)
set(findall(gcf,'-property','FontSize'),'FontSize',14);
time_shifted = circshift(time,1);
time_shifted(1) =0;
time_dif = time-time_shifted;
figure,
plot(h(1:kk),'LineWidth',2)
grid on
xlabel('iteration')
ylabel('time step (days)')
figure,
plot(time_dif,'LineWidth',2)
grid on
xlabel('iteration')
ylabel('time step (days)')