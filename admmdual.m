function [x,y,s,prim,dual,iter,secs,errs,valps,valds] = admmdual(A,b,c)
% x the prime variable
% y the dual variable
% s the relax variable
% prim the value of prime problem
% dual the value of dual problem
% iter the number of iteration
% secs the running time

% set parameter
max_iter=10000; % max iteration number
tol = 1e-6; % the accuracy
sigma = 0.618; % the step size rate

% check input dimension
[m,n] = size(A);
assert( length(b)==m, 'mismatch rows of A and length of b');
assert( length(c)==n, 'mismatch columns of A and length of c');
normb = norm(b); normc = norm(c);
y = randn(m,1);
x = randn(n,1); 
s = randn(n,1);

% initial matrix setting
AAt = A*A';

%iteration setting
done = 0;  % stopping condition not satisfied
iter = 0;  % iteration count
incsig = 0;% counts how often sigma is increased consecutively
redsig = 0;% counts how often sigma is reduced consecutively

% the print 
fprintf('\n     iter  | time |    dual         prim       sigma |');
fprintf('\n')

% time start
tstart = tic;
% the main process
x_old = x;
while done == 0
    % fix x,s and update y
    y = -AAt\((A*x-b)/sigma+A*(s-c));   
    % fix y,x and update s by projection
    s = max(c-A'*y-x/sigma,0); 
    % fix y,s and update x
    x = x + (A'*y+s-c)*sigma;
    % check the kkt condition
    err(1) = norm(b - A*x)/(1+normb);
    err(2) = norm(c - A'*y -s)/(1+normc);
    normx = norm(x);
    err(3) = norm(x - max(x,0))/(1+normx);
    err(4) = norm(s'*x)/(1+normx+norm(s));
    err(5) = norm(x-x_old)/(1+normx);
    x_old = x;
    
    %updata parameter sigma
    if (err(1) > 10*err(2))
        redsig = redsig + 1;
    elseif (err(1) < 2*err(2))
        incsig = incsig + 1;
    elseif (err(2) < tol*3)
        sigma = sigma*9;
    else
        incsig = 0; redsig = 0;
    end
    
    % update sigma according to redsig and incsig
    if incsig > 9 
        sigma = sigma*2; incsig = 0;
    end
    if redsig > 9
        sigma = sigma/2; redsig=0;
    end
    
    %check stop condition
    iter = iter + 1;
    secs = toc(tstart);
    prim = c'*x; dual = b'*y;
    done = (max(err)<tol) || (iter>= max_iter);
    errs(iter) = max(err(1),err(2)); %#ok<AGROW>
    valps(iter) = prim; %#ok<AGROW>
    valds(iter) = dual; %#ok<AGROW>
    % done = done || (norm(prim-dual)<1e-6);
    
    %print the iteration (mod 50)
    if (mod(iter, 50)==0)
        fprintf( '%8.0d %8.4f %13.6e %13.6e %9.6f \n', ...
        iter, secs, dual, prim, sigma);
    end
    
end

% the output we solve this problem 
if max(err)<tol
  fprintf('required accuracy reached after %8.0d iterations,',iter)
else
  fprintf('max iterations ');
end
fprintf('total time: %8.2f \n', secs);

end