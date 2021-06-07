function [A,b,c] = standard(C,mu,v)
% convert into standard form
% C: cost matrix; mu and v: distribution
% A, b, c is the standard form in LP
% min c^t x 
% s.t. Ax=b, x>=0
% =================================
[m,n] = size(C);
c = reshape(C,m*n,1);
b = [mu;v];
A = sparse(m+n,m*n);
for i = 1:m
    Aa = zeros(m,n);
    Aa(i,:) = ones(1,n);
    A(i,:) = reshape(Aa,1,m*n);
end
for i = 1:n
    Aa = zeros(m,n);
    Aa(:,i) = ones(m,1);
    A(i+m,:) = reshape(Aa,1,m*n);
end
end
