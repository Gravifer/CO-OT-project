function [xx,secs] = Gurobispx(A,b,c)
% Simplex methods with solver GUROBI
% xx: the solution
% secs: the running time
% =====================
clear model;
sense = [];
for ii = 1:size(A,1)
    sense=[sense;'=']; %#ok<AGROW>
end
model.A = sparse(A);
model.obj = c;
model.modelsense = 'Min';
model.rhs = b;
model.sense = sense;
% Options are:0=primal simplex, 1=dual simplex, 2=barrier
model.Method = 2;
% time start
tstart = tic;
% apply gurobi
result = gurobi(model);
% time end
secs = toc(tstart);
% solution
xx = result.x;
end

