function y = vec2STATS (x,s,n)
% function y = vec2STATS (x,s,n)
% converts data in vector format to pseudo-staSTATS format (i.e. structured
% cell array)
% x = ncond*nvar matrix of means
% s = ncond*nvar matrix of standard deviations
% n = ncond*nvar matrix of n's
%
ncond = size(x,1); nvar = size(x,2);
if nargin < 3, n = ones(ncond,nvar); end
if nargin < 2, s = ones(ncond,nvar); end
y = cell(1,nvar);
for ivar=1:nvar
    y{ivar}.means = x(:,ivar);
    y{ivar}.n = diag(n(:,nvar)); % put in matrix format
    y{ivar}.cov = diag(s(:,ivar).^2); % convert to variances
    y{ivar}.regcov = y{ivar}.cov; % assumes no shrinkage
    y{ivar}.shrinkage=0; % like I said
    y{ivar}.weights = y{ivar}.n.*y{ivar}.regcov^-1; % now we have weights
    y{ivar}.lm = y{ivar}.regcov; % not used but may as well park something there
    y{ivar}.shrink=0;
end
