function [x, f, stats] = staCMR (data, E, shrink, tolerance, proc)
% function [x, f, stats] = staCMR (data, E, shrink, tolerance, proc)
% conducts conjoint monotonic regression for STA
% data is cell array of data or structured output from staSTATS 
% E is an optional partial order
% shrink is parameter to control shrinkage of covariance matrix:
% 0 = no shrinkage; 1 = diagonal matrix; -1 = calculate optimum
% tolerance (reserved for future implementation - set to zero)
% proc is the number processors available: default (-1) use all available
% returns:
% x = best fitting CMR values to data means
% f = fit statistic
% stats = structured variable of statistics from fit procedure
% this includes:
% shrinkage = returned covariance shrinkage values for each group if
% computed
% nitr = number of iterations
%
% ************************************************************
% last modified: 10 February 2017
% checked 14 August 2017
% ************************************************************
%
% enter defaults
if nargin < 5, proc = -1; end
if nargin < 4, tolerance=0; end
if nargin < 3, shrink=-1; end
if nargin < 2, E = {}; end
if isempty(E)
    E = {}; 
elseif ~iscell(E) && isvector(E)
    E = {E};
end
if isempty(shrink), shrink=-1; end
if isempty(tolerance), tolerance=0; end
if isempty(proc), proc=-1; end

% call function
[x, f, stats] = staCMRx (data, [], E, shrink, tolerance, proc);

function [x, f, stats] = staCMRx (data, model, E, shrink, tolerance, proc, approx)
% function [x, f, stats] = staCMRx (data, model, E, shrink, tolerance, proc, approx)
% Multidimensional CMR
% data is cell array of data or structured output from staSTATS 
% model is a nvar*k matrix specifying the linear model; default = ones(nvar,1)
% E is partial order
% shrink is parameter to control shrinkage of covariance matrix;
% 0 = no shrinkage; 1 = diagonal matrix; -1 = calculate optimum
% returns:
% x = best fitting CMR values to y-means
% f = fit statistic
% shrinkage = estimated shrinkage of covariance matrix
% *************************************************************************
% previously modified 22 January 2016
% modified 12 April 2016 to circumvent null rows in model
% modified 3 June 2016 to deal with cell array data format
% modified 25 Aug 2016 to rectify a few odd things
% modified 12 Nov 2016 to include tolerance (0 to 1)
% *************************************************************************

if nargin < 7, approx = 0; end
if nargin < 6, proc = -1; end
if nargin < 5, tolerance=0; end
if nargin < 4, shrink=-1; end
if nargin < 3, E = {}; end
if nargin < 2, model=[]; end

if isempty(E)
    E = {}; 
elseif ~iscell(E) && isvector(E)
    E = {E};
end
if isempty(shrink), shrink=-1; end
if isempty(tolerance), tolerance=0; end
if isempty(proc), proc=-1; end

if iscell(data)
    if isstruct(data{1})
        y = data; % if structured then already in stats form
    else
        y = staSTATS(data, shrink); % otherwise assume within-subjects data and get stats
    end
else
    celldata = gen2cell (data); % convert from  "general format"
    y = staSTATS (celldata, shrink);
end

nvar = numel(y); 
if isempty(model), model = ones(nvar,1); end % default model for STA

shrinkage = zeros(numel(y{1}.shrinkage),numel(y)); 
for ivar=1:nvar, shrinkage(:,ivar) = y{ivar}.shrinkage; end

% decompose model if there are null rows to short circuit result
x = cell(1,nvar);
s = any(model,2); i = find(s==1);

nitr=[]; adjStar=[]; fl=[]; fu=[];fuf=[]; remaining=[];
if rank(model(i,:))==numel(i)
    xx = cell(1,numel(i)); f = 0;
    for k=1:numel(i)
        xx{k} = y{i(k)}.means';
    end
else
    if approx==0
        [xx, f, nitr, adjStar, fl, fu, fuf, remaining] = jCMRx (y(i), model(i,:), E, proc, tolerance); % call CMRx for non-zero rows of model
    else
        [xx, f] = jCMRx1 (y(i), model(i,:), E);
    end
end

x(i) = xx;
j = find(s==0); % replace zero rows with weighted means
if ~isempty(j)
    for k=1:numel(j)
        u = j(k);
        w = y{u}.weights;
        if ~isvector(w)
            w = diag(w); if iscolumn(w), w=w'; end
        end
        m = y{u}.means; if isrow(m), m=m'; end
        z = w*m/sum(w);
        x{u} = z*ones(size(m));
        d = m - x{u};
        f = f + d'*diag(w)*d;
    end
end
if f < 0, f = 0; end

if approx==0
    stats.shrinkage = shrinkage;
    stats.nitr = nitr;
    stats.adj = adjStar;
    stats.fl = fl; 
    stats.fu = fu;
    stats.fuf = fuf;
    stats.remaining = remaining;
else
    stats=[];
end
