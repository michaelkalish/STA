function [p, datafit, fits, pars] = staCMRFITBN (data, varargin)
% [p, datafit, fits, pars] = staCMRFITBN (data, varargin)
% Binomial version of staCMRFIT
% nsample = no. of Monte Carlo samples (about 10000 is good)
% data = data structure (cell array, general, or structured)
% model is a nvar * k matrix specifying the linear model; default = ones(nvar,1)
% E = optional partial order model e.g. E={[1 2] [3 4 5]} indicates that
% condition 1 <= condition 2 and condition 3 <= condition 4 <= condition 5
%
% returns:
% p = empirical p-value
% datafit = observed fit of monotonic (1D) model
% fits = nsample vector of fits of Monte Carlo samples (it is against this
% distribution that datafit is compared to calculate p)
% % Note: These are g-squared values (not least squares)
% % pars = nvar cell array of bootstrap means nsample x ncond x nsub
% *************************************************************************
% Last modified: 24 April 2017
% pars output added 4 March 2020
% *************************************************************************
%
if ~iscell(data)
    y = BNgen2cell(data);
else
    y = data;
end

% defaults
nsample = 1;
E = {};
nvar = size(y,2);
model = ones(nvar,1);
approximate = 0;
showStatus = 0;
ranseed = -1;
tol = 1e-6;

% set arguments if specified
for i = 1 : 2 : length(varargin)-1
    name = varargin{i}; 
    value = varargin{i+1};
    switch name
        case {'nsamples', 'nsample', 'ns', 'n'} 
            nsample = value; % number of samples
        case {'partial', 'part', 'p'}
            E = value; % partial order
        case {'model','mod','m'}
            model = value; % model
        case {'approximate','approx','a'}
            approximate = value; % approximate solution if true
        case {'showstatus','show','s'}
            showStatus = value; % Luke's update window
        case {'ranseed','ran','r'}
            ranseed = value; % random number seed
    end
end

if ~iscell(E)
    E = adj2cell(E); % convert from adjacency matrix form
end

[p, datafit, fits, pararray] = jCMRxBNfits(nsample, y, E, model, approximate, showStatus, ranseed);
% additional changes March 2020
datafit(datafit<tol) = 0;
fits (fits < tol) = 0;
pars = cell(1,nvar); % pars contains bootstrap means
for ivar=1:nvar
	z = squeeze(pararray(:,:,ivar,:));
    pars{ivar} = permute(z, [1 3 2]);
end

