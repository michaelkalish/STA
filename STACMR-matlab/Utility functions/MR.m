function [x, fit] = MR (y, W, E, x0)
% function [x, fit] = MR (y, W, E, x0)
% fits basic monotonic regression model to vector y with weight matrix W
% according to partial order E
% x0 is optional starting point
% returns:
% x = best fitting MR values to y-means
% fit = least squares fit
%
% last updated: 19 March 2014
% uses quadprog instead of lsqlin to increase speed (a bit)
%
% Modified 29 July 2015 to be algorithm of choice (because of var-cov
% weight matrix)
%
if isrow(y) % convert data to column vector
    y = y';
end
switch nargin
    case 3
        x0 = repmat(mean(y),numel(y),1);  % starting point is a vector of means
    case 2
        E = {}; x0 = repmat(mean(y),numel(y),1); 
    case 1
        W = eye(numel(y)); E = {}; x0 = repmat(mean(y),numel(y),1);  % default identity matrix
end
if isempty(E)
    E = {};
end
if isempty (x0)
    x0 = repmat(mean(y),numel(y),1);
end
if isempty(W)
    W = eye(numel(y));
end
if isvector(W)
    W = diag(W); % convert vector to matrix if required
end

if iscell(E)
    adj = cell2adj(1:numel(y),E); % convert to adjacency matrix
    if sum(sum(adj ~= 0)) > 0
        A = adj2ineq (adj); % convert adjacency matrix to a set of inequalities
    else
        A = [];
    end
else
    adj = E;
    if sum(sum(adj == -1)) > 0 % check if already an inequality matrix
        A = adj;
    else
        A = adj2ineq (adj); % convert adjacency matrix to a set of inequalities
    end
end

% define b matrix
if ~isempty(A)
    b = zeros(size(A,1),1);
else
    b = [];
end

% do monotonic regression 
%options = optimset ('LargeScale','off', 'display','off'); % turn off display and use medium algorithm to avoid warning
%C = sqrtm(W);
%d = C*y;
%[x, fit] = lsqlin(C,d,A,b,[],[],[],[],x0,options);
if sum(sum(adj ~= 0)) == 0
    x = y;
else
    options = optimset ('LargeScale','off', 'display','off', 'algorithm','interior-point-convex'); % turn off display and use medium algorithm to avoid warning
    u = -y'*W;
    x = quadprog(W,u,A,b,[],[],[],[],x0,options);
end
fit = (x - y)'*W*(x - y);