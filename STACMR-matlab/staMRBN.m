function [x, G2, f] = staMRBN (data, E)
% function [x, G2, f] = staMRBN (data, E)
% fits partial order monotonic regression model to binomial data
% data = nsubj x nvar cell array of ncond x 2 matrices of counts (hits, misses) 
% where ncond = no. of conditions
% E is the partial order specified as either:
% (a) cell array, 
% (b) adjacency matrix,
% (c) inequality matrix
% returns:
% x = nsubj x nvar cell array of ncond model means (proportions)
% G2 = nsub-vector of g-squared fit of model
% f = nsub-vector of least squares fit of model
%
% *************************************************************************
% Revised 20 June 2015 following discussions with Oleg Sysoev
% revised 10 March 2017 to incorporate java commands
% *************************************************************************
%
tol = 1e-10;
if nargin < 2; E = []; end
if isempty(E), E = {}; end
if ~iscell(E) && isvector(E), E = {E}; end
if isempty(E)
    disp ('Warning: Partial order undefined.');
end

if ~iscell(data)
    y = BNgen2cell(data); y = binSTATS(y);
elseif isstruct(data{1})
    y = data;
else
    y = binSTATS (data);
end

x = cell(size(y)); f = zeros(size(y)); ml = f; 
for isub = 1:size(y,1)
    for ivar = 1:size(y,2)
        W = y{isub,ivar}.weights; if isvector(W), W = diag(W); end
        [x{isub,ivar}, f(isub,ivar)] = jMR (y{isub,ivar}.means, W, E);
        x{isub,ivar}(x{isub,ivar}<0) = 0;
        x{isub,ivar}(x{isub,ivar}>1) = 1;
        ml(isub,ivar) = mleBN(y{isub,ivar}.count,x{isub,ivar});
    end
end
f = sum(f,2); % sum over variables
G2 = sum(ml,2); G2(G2<=tol)=0;




