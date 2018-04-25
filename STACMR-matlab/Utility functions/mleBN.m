function f = mleBN (n, x)
% function f = mleBN (n, x)
% returns G-squared derived from maximum likelihood
% n = ncond x 2 matrix of counts (hits, misses) where ncond = no. of conditions
% x = ncond x 1 vector of model means (i.e., proportions)
% works as a likelihood ratio test of difference between model and
% saturated (data) model
nn = sum(sum(n));
x(x<=0)=.5/nn; % replace expected zeros with small number
x(x>=1)=(nn-.5)/nn; % replace expected ones with number close to 1
a = [x 1-x]; 
a = a.*repmat(sum(n,2),1,size(a,2)); % expected counts
g2 = 2*n.*log(n./a); g2(n<=0) = 0; % g-squared
f = sum(sum(g2));
f(f<0)=0;
