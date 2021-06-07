%% Example 1: Compute the distribution translation (See the ¡±Test data¡± file).
% clear; 
% clc;
% load('shift_of_Ricker.mat')
mu = fx';
mu = mu/sum(mu);
v = fy';
v = v/sum(v);
% % convert into standard form
% [A,b,c] = standard(C,mu,v);
% 
% % Interior point methods with solver MOSEK
% [x_m_i,sec_m_i] = Mosekint(A,b,c);
% fprintf('\n')
% % Simplex methods with solver MOSEK
% [x_m_s,sec_m_s] = Mosekspx(A,b,c);
% fprintf('\n')
% % Barrier methods with solver GUROBI
% [x_g_b,sec_g_b] = Gurobibar(A,b,c);
% fprintf('\n')
% % Interior point methods with solver GUROBI
% [x_g_s,sec_g_s] = Gurobispx(A,b,c);
% fprintf('\n')
% % the alternating direction method of multipliers for dual problem
% [x,y,s,prim,dual,iter,secs,errs,valps,valds] = admmdual(A,b,c);
% fprintf('\n')
% the sinkohorn 
iter = 2000;
% epsilon_sinkhorn = 0.03;
% [sinkhorn_x,sinkhorn_opt,sinkhorn_opts,secs_s] = sinkhorn(C,mu,v,iter,epsilon_sinkhorn);
% fprintf('\n')
% ipot
epsilon_ipot = 0.7;
[ipot_x,ipot_opt,ipot_opts,secs_i] = ipot(C,mu,v,iter,epsilon_ipot);



%% Example 2: 
% clear;
% clc;
% m = 1000;
% n = 1000;
% C = 1+9*rand(m,n);
% mu = rand(m,1);
% mu = mu/sum(mu);
% v = rand(n,1);
% v = v/sum(v);
% % convert into standard form
% [A,b,c] = standard(C,mu,v);
%
% % Interior point methods with solver MOSEK
% [x_m_i,sec_m_i] = Mosekint(A,b,c);
% fprintf('\n')
% % Simplex methods with solver MOSEK
% [x_m_s,sec_m_s] = Mosekspx(A,b,c);
% fprintf('\n')
% % Barrier methods with solver GUROBI
% [x_g_b,sec_g_b] = Gurobibar(A,b,c);
% fprintf('\n')
% % Interior point methods with solver GUROBI
% [x_g_s,sec_g_s] = Gurobispx(A,b,c);
% fprintf('\n')
% % the alternating direction method of multipliers for dual problem
% [x,y,s,prim,dual,iter,secs,errs,valps,valds] = admmdual(A,b,c);
% fprintf('\n')
% % the alternating direction method of multipliers for dual problem
% [x,y,s,prim,dual,iter,secs,errs,valps,valds] = admmdual(A,b,c);
% fprintf('\n')
% % the sinkohorn 
% iter = 2000;
% epsilon_sinkhorn = 0.03;
% [sinkhorn_x,sinkhorn_opt,sinkhorn_opts,secs] = sinkhorn(C,mu,v,iter,epsilon_sinkhorn);
% fprintf('\n')
% % ipot
% epsilon_ipot = 0.7;
% [ipot_x,ipot_opt,ipot_opts,sec] = ipot(C,mu,v,iter,epsilon_ipot);

