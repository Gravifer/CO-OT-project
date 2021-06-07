function [P,opt,opts,secs] = sinkhorn(C,mu,v,iter,epsilon)
% C is the cost matrix
% mu and v are the distribution
% iter is the maximun number of iteration
% epsilon is the normalize parameter
% opts is the objective value
[~,n] = size(C);
b = ones(n,1);
K = exp((-1/epsilon).*C);
ts = tic;
for t = 1:iter
    a = mu./(K*b);
    b = v./(K'*a);
    opts(t) = trace(C'*(diag(a)*K*diag(b))); %#ok<AGROW>
end
P = diag(a)*K*diag(b);
opt = trace(C'*P);
secs = toc(ts);
end