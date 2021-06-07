function [T,opt,opts,secs] = ipot(C,mu,v,iter,epsilon)
% C is the cost matrix
% mu and v are the distribution
% iter is the maximun number of iteration
% epsilon is the normalize parameter
% opts is the objective value
[m,n] = size(C);
b = ones(n,1);
G = exp((-1/epsilon).*C);
T = ones(m,n);
ts = tic;
for t = 1:iter
    Q = G.*T;
    for i = 1:2
        a = mu./(Q*b);
        b = v./(Q'*a);
    end
    T = diag(a)*Q*diag(b);
    opts(t) = trace(C'*T); %#ok<AGROW>
end
opt = trace(C'*T);
secs = toc(ts);
end