function y = gen2cell(data)
% function y = gen2cell(data)
% converts data in "general format" to cell array suitable for input to
% staSTATS
% general format is defined as:
% column 1 = subject number (nsub)
% column 2 = between-subjects condition (ngroup)
% column 3 = dependent variable (nvar)
% columns 4 to end = values for each within-subjects condition (ncond)
% output (ys) is ngroup x nvar cell array in which each element is an nsub
% x ncond matrix of values.
%
% *************************************************************************
% written 18 August 2016
% *************************************************************************
%
group = data(:,2); ugroup = unique(group); ngroup = numel(ugroup);
var = data(:,3); uvar=unique(var); nvar = numel(uvar);
within = data(:,4:end); 

y = cell(ngroup, nvar);
for igroup = 1:ngroup
    for ivar = 1:nvar
        k = group==ugroup(igroup) & var==uvar(ivar);
        y{igroup, ivar} = within(k,:);
    end
end

    
