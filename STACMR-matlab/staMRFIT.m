function [p, datafit, fits, maxbad, times] = staMRFIT (data, varargin)
% function [p, datafit, fits] = staMRFIT (data, varargin)
% fits the MR model
% nsample = no. of Monte Carlo samples (about 10000 is good)
% data = data structure (cell array or general)
% E = optional partial order model e.g. E={[1 2] [3 4 5]} indicates that
% condition 1 <= condition 2 and condition 3 <= condition 4 <= condition 5
% shrink is parameter to control shrinkage of covariance matrix (if input is not stats form);
% 0 = no shrinkage; 1 = diagonal matrix; -1 = calculate optimum
% returns:
% p = empirical p-value
% datafit = observed fit of monotonic (1D) model
% fits = nsample vector of fits of Monte Carlo samples (it is against this
% distribution that datafit is compared to calculate p)
% *************************************************************************
% Last modified: 15 January 2017
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

nsample = 1;
E = {};
shrink = -1;
reverse = 0;
proc = -1;
mrTol = 0;
ranseed = -1;
showStatus = 0;

for i = 1 : 2 : length(varargin)-1
    name = varargin{i}; 
    value = varargin{i+1};
    switch name
        case {'nsamples', 'nsample', 'ns', 'n'} 
            nsample = value; % number of samples
        case {'partial', 'part', 'p', 'order', 'ord', 'o'}
            E = value; % partial order
        case {'shrink', 's'}
            shrink = value; % shrink value
        case {'reverse', 'rev'}
            reverse = value; % reverse value
        case {'workers', 'work', 'w', 'proc'}
            proc = value; % no. of workers
        case {'mrTol', 'tol', 't'}
            mrTol = value; % approximate solution option
        case {'ranseed','ran','seed'}
            ranseed = value; % random number seed
        case {'showstatus','show','sh'}
            showStatus = value; % showstatus option
    end
end
if isempty(E)
    disp('Warning: Partial order undefined.');
end
if ~iscell(E)
    E = adj2cell(E); % convert from adjacency matrix form
end
[p, datafit, fits, maxbad, times] = jMRfits(nsample, y, E, shrink, reverse, proc, mrTol, ranseed, showStatus);
