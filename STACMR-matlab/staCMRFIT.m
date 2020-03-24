function [p, datafit, fits, pars, maxbad, times] = staCMRFIT (data, varargin)
% [p, datafit, fits, pars, maxbad, times] = staCMRFIT (data, varargin)
% multidimensional version of CMRfits. Replaces CMRfitsx.
% nsample = no. of Monte Carlo samples (about 10000 is good)
% data = data structure (cell array, general, or structured)
% model is a nvar * k matrix specifying the linear model; default = ones(nvar,1)
% E = optional partial order model e.g. E={[1 2] [3 4 5]} indicates that
% condition 1 <= condition 2 and condition 3 <= condition 4 <= condition 5
% shrink is parameter to control shrinkage of covariance matrix (if input is not stats form);
% 0 = no shrinkage; 1 = diagonal matrix; -1 = calculate optimum
% returns:
% p = empirical p-value
% datafit = observed fit of monotonic (1D) model
% fits = nsample vector of fits of Monte Carlo samples (it is against this
% distribution that datafit is compared to calculate p)
% pars = nvar cell array of bootstrap means
% *************************************************************************
% Last modified: 24 February 2017
% additional changes 4 March 2020
% *************************************************************************
%
if iscell(data)
    if isstruct(data{1})
        y = data; 
    else
        y = data; % cell array format
    end
else
    y = gen2cell(data); % convert from general format
end
nvar = size(y,2);
model = ones(nvar,1);
nsample = 1;
E = {};
shrink = -1;
proc = -1;
cheapP = 0;
approximate = 0;
mrTol = 0;
ranseed = -1;
showStatus = 0;
tol = 1e-6;

for i = 1 : 2 : length(varargin)-1
    name = varargin{i}; 
    value = varargin{i+1};
    switch name
        case {'model', 'mod', 'm'}
            model = value; % model
        case {'nsamples', 'nsample', 'ns', 'n'} 
            nsample = value; % number of samples
        case {'partial', 'part', 'p', 'order', 'ord', 'o'}
            E = value; % partial order
        case {'shrink', 's'}
            shrink = value; % shrink value
        case {'workers', 'work', 'w'}
            proc = value; % no. of workers
        case {'cheapP', 'cheap', 'c'}
            cheapP = value; % cheapP option
        case {'approximate', 'approx', 'a'}
            approximate = value; % approximate solution option
        case {'mrTol', 'tol', 't'}
            mrTol = value; % approximate solution option
        case {'ranseed','ran','r'}
            ranseed = value; % random number seed
        case {'showstatus','show','sh'}
            showStatus = value; % showstatus option
    end
end

if ~iscell(E)
    E = adj2cell(E); % convert from adjacency matrix form
end

[p, datafit, fits, pararray, maxbad, times] = jCMRfitsx(nsample, y, model, E, shrink, proc, cheapP, approximate, mrTol, ranseed, showStatus);
% additional changes March 2020
datafit(datafit<tol) = 0;
fits(fits < tol) = 0;
pars = cell(1,nvar); % pars contains bootstrap means
for ivar=1:nvar
    pars{ivar} = squeeze(pararray(:,ivar,:));
end



