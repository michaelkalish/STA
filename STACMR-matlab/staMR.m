function [xPrime, fit, shrinkage] = staMR (data, E, shrink, reverse)
% function [xPrime, fit, shrinkage] = staMR (data, E, shrink, reverse)
% fits monotonic regression model to data according to partial order
% data is cell array of data or structured output from staSTATS
% E is partial order
% shrink is parameter to control shrinkage of covariance matrix;
% 0 = no shrinkage; 1 = diagonal matrix; -1 = calculate optimum
% returns:
% x = best fitting MR values to y-means
% f = total fit statistic
% *************************************************************************
% Last updated: 4 Sept 2015
% modified 3 June 2016 to deal with cell array data format
% *************************************************************************
%
tol = 1e-10;
if nargin < 4, reverse = false; end
if nargin < 3, shrink = -1; end
if nargin < 2, E = {}; end

if isempty(E), E = {}; end
if ~iscell(E), E = {E}; end
if isempty(shrink), shrink = -1; end
if isempty(reverse), reverse = 0; end

if iscell(data)
    if isstruct(data{1})
        y = data; % if structured then already in stats form
    else
        y = staSTATS(data, shrink); % otherwise assume within-subjects data and get stats
    end
else
    if isstruct(data) 
        y = cell(1,1); y{1} = data; % data is structured output for one variable
    else
        celldata = gen2cell (data); % convert from old "general format"
        y = staSTATS (celldata, shrink);
    end
end

nvar = numel(y);
shrinkage = zeros(numel(y{1}.shrinkage),nvar); 
for ivar=1:nvar, shrinkage(:,ivar) = y{ivar}.shrinkage; end

xPrime = cell(1,nvar); 
fits = zeros(nvar,1); 
for ivar = 1:nvar
    [xPrime{ivar}, fits(ivar)] = jMR (y{ivar}.means, y{ivar}.weights, E, reverse);
end
fit = sum(fits);
if fit < tol, fit=0; end



