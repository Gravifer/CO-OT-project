function [xx,secs] = Mosekspx(A,b,c)
% Simplex methods with solver MOSEK
% xx: the solution
% secs: the running time
% =====================
blx=zeros(length(c),1);
% time start
tstart = tic;
% apply Mosek
[res] = msklpopt(c,A,b,b,blx,[],[],'minimize');
% time end
secs = toc(tstart);
% Basic solution
sol = res.sol;
xx = sol.bas.xx;
end