function output = BNgen2cell (data)
% function output = BNgen2cell (data)
% converts a data in general format into a cell array binomial data
% suitable for staCMRBN.m
%
if iscell(data)
    disp('Error: Data must be a matrix.');
    output = [];
else
    groups = data(:,1); ugroup = unique(groups); ngroup = length(ugroup);
    vars = data(:,3); uvar = unique(vars); nvar = length(uvar);
    counts = data(:,4:5);
    output = cell(ngroup,nvar);

    for igroup = 1:ngroup
        for ivar = 1:nvar
            k = groups==ugroup(igroup) & vars==uvar(ivar);
            output{igroup,ivar} = counts(k,:);
        end
    end
end
