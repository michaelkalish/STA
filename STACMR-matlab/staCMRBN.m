function [x, G2, f, nitr] = staCMRBN(data, E, model, approximate, proc)
% function [x, G2, f, nitr] = staCMRBN (data, E)
% fits partial order model to binomial data
% data = nsubj x nvar cell array of ncond x 2 matrices of counts (hits, misses) where ncond = no. of conditions
% E is the partial order specified as either (a) a cell array, (b) adjacency matrix
% returns:
% x = nsubj x nvar cell array of ncond model means (i.e., proportions)
% G2 = nsubj vector of G-square maximimum likelihood fits of model
% f = nsubj vector of least squares fits of model
% nitr = number of iterations
%
% *************************************************************************
% Revised 20 June 2015 following discussions with Oleg Sysoev
% modified again on 12 November 2015
% updated 10 March 2017 to incorporate java function
% updated 24 April 2017
% *************************************************************************
%
if ~iscell(data)
    y = BNgen2cell(data);
else
    y = data;
end
nvar = size(y,2);
if nargin < 5, proc = -1; end
if nargin < 4, approximate = 0; end
if nargin < 3, model = ones(nvar,1); end
if nargin < 2, E = {}; end

if isempty(E)
    E = {};
end
if ~iscell(E) && isvector(E)
    E = {E};
end

[x, G2, f, nitr] = jstaCMRBN(y, E, model, approximate);





