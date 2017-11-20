function output = staSTATSBN(data)
% function output = staSTATSBN(data)
% returns statistics on NSUB x NVAR cell array of binomial data
% where each element is a NCOND x 2 matrix (hits, misses)
% output is a NSUB x NVAR cell array:
% output.means = observed means
% output.n = number of observations
% output.weights = weights for monotonic regression
%
% *************************************************************************
% Revised 20 June 2015 following discussions with Oleg Sysoev
% Effectively returning binSTATS to previous verions (binSTATS_old)
% minor update 25 April 2017
% *************************************************************************
%
y = data;
if ~iscell(data)
    y = BNgen2cell(data);
end
output = cell(size(y));
for i=1:size(y,1)
    for j=1:size(y,2)
        out.count = y{i,j};
        n = sum(y{i,j},2);
        out.means = y{i,j}(:,1)./n;
        out.weights = n;
        out.n = n;
        output{i,j} = out;
    end
end
