function [p, datafit, fits, pars] = staMRFITBN (data, varargin)
% [p, datafit, fits, pars] = staMRFITBN (data, varargin)
% Binomial version of staMRFIT
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
% distribution that datafit is compared to calculate p).
% Note: These are g-squared values (not least squares)
% pars = nsub x nvar cell array of bootstrap means nsample x ncond 
% *************************************************************************
% Last modified: 24 February 2017
% updated 10 March 2017
% modified to return pars 24 March 2020
% *************************************************************************
%
tol = 1e-10;
if ~iscell(data)
    y = BNgen2cell(data);
else
    y = data;
end
nsub = size(y,1);
nvar = size(y,2);
nsample = 1;
E = {};

for i = 1 : 2 : length(varargin)-1
    name = varargin{i}; 
    value = varargin{i+1};
    switch name
        case {'nsamples', 'nsample', 'ns', 'n'} 
            nsample = value; % number of samples
        case {'partial', 'part', 'p', 'order','ord', 'o'}
            E = value; % partial order
    end
end

if ~iscell(E)
    E = adj2cell(E); % convert from adjacency matrix form
end

if isempty(E)
    disp ('Warning: Partial order undefined.');
end

[p, datafit, fits, pararray] = jMRBNfits(nsample, y, E); % returns G2 values as datafit (and fits)

% clear up near zero values
datafit(datafit<tol)=0;
fits(fits<tol)=0;
% unpack bootstrap means
% pararray is nsample x nsub x nvar x ncond array
pars = cell(nsub,nvar); % pars contains bootstrap means
for ivar=1:nvar
    for isub = 1:nsub
        z = squeeze(pararray(:,isub,ivar,:));
        pars{isub,ivar} = z;
    end
end
