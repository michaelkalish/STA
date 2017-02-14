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
% ************************************************************
%
% enter defaults
if nargin < 5, proc = -1; end
if nargin < 4, tolerance=0; end
if nargin < 3, shrink=-1; end
if nargin < 2, E = {}; end

if isempty(E), E = {}; end
if isempty(shrink), shrink=-1; end
if isempty(tolerance), tolerance=0; end
if isempty(proc), proc=-1; end

if ~iscell(E) && isvector(E)
    E = {E};
end
% call function
[x, f, stats] = staCMRx (data, [], E, shrink, tolerance, proc);

