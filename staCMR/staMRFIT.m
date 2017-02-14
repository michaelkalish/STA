function [p, datafit, fits] = staMRFIT (nsample, data, E, shrink)
% function [p, datafit, fits] = staMRFIT (nsample, data, E, shrink)
% *************************************************************************
% Tests partial order or "trace" model specified by partial order, E,
% against unconstrained model
% Uses non-parametric bootstrap as suggested by Oleg Sysoev.
%
% computes empirical p-value for hypothesis that means across NCOND
% conditions are monotonically increasing
% nsample = no. of Monte Carlo samples (about 10000 is good)
% data = data - either within, general, or summary formats
% E = optional partial order model e.g. E={[1 2] [3 4 5]} indicates that
% condition 1 <= condition 2 and condition 3 <= condition 4 <= condition 5
% returns:
% p = empirical p-value
% datafit = observed fit of monotonic (1D) model
% fits = nsample vector of fits of Monte Carlo samples (it is against this
% distribution that datafit is compared to calculate p)
%
% *************************************************************************
% Last modified: 8 March 2016 to include summary data
%                3 Jan 2017 to reconcile with jMR
% *************************************************************************
%
if nargin < 4, shrink = -1; end

if iscell(data)
    if isstruct(data{1})
        type = -1; % data in summary format
        ys = data;
    else
        type = 0; % data in cell array format
        ys = staSTATS(data, shrink); % get stats
    end
else
    type = 1; % data in "general format" i.e. text file with one line per subject
    ys = staSTATS(data, shrink); % get stats
end
[~, datafit] = staMR (ys, E, shrink); % calculate observed data fit (difference between partial order and unconstrained models)

% initialise random number generator (to random start)
v = version; 
if str2double(v(1)) > 5
    rng ('default'); rng ('shuffle');
else
    rand('seed',sum(100*clock)); % legacy code for old versions of Matlab
end

% do non-parametric bootstrap resampling to estimate p-value
fits = zeros(nsample,1);
for isample=1:nsample
    datab1 = bootstrap (data, type); % bootstrap sample from data
    x = staMR (datab1, E, shrink); % fit partial order model
    datam = movedata(data, ys, x); % move data to best fitting partial order model means
    datab2 = bootstrap(datam, type); % bootstrap sample from partial order model data
    [~, fits(isample)] = staMR(datab2, E, shrink); % calculate and store observed data fit (difference between partial order and unconstrained models)  
end
p = sum(fits >= datafit)/nsample; % calculate p
 
function ym = movedata (data, ys, x1)
% function ym = movedata (data, ys, x1)
% moves data from means specified by ys to means specified by x1
if iscell(data)
    if isstruct(data{1})
        type = -1; % summary stats
    else
        type = 0; % within-subjects format
    end
else
    type = 1; % general format
end
ym = data; 
if type == 0 % cell array format
    for ivar = 1:size(data,2)
        i2 = 0;
        for igroup=1:size(data,1)
            a = data{igroup,ivar}; 
            i1 = i2 + 1; i2 = i1 + size(a,2) - 1;
            m = repmat(ys{ivar}.means(i1:i2),size(a,1),1);
            q = repmat(x1{ivar}(i1:i2)',size(a,1),1); 
            ym{igroup,ivar} = a - m + q;
        end
    end
elseif type == 1 % general format
    cond = data(:,2); ucond = unique(cond); var = data(:,3); uvar=unique(var);
    within = data(:,4:end);
    for ivar = 1:numel(uvar)
        j = 1;
        for icond = 1:numel(ucond)
            k = find(var==uvar(ivar) & cond==ucond(icond));
            x = within(k,:); n = size(x,2);
            m = repmat(ys{ivar}.means(j:j+n-1),numel(k),1);
            ym(k,4:end) = x - m + repmat(x1{ivar}(j:j+n-1)',numel(k),1);
            j = j + n;
        end
    end
else
    ym = data; % summary format
    for ivar=1:numel(ys)
        ym{ivar}.means = x1{ivar};
    end
end


