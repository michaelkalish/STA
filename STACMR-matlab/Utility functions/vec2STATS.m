function y = vec2STATS (x)
% function y = vec2STATS (x)
% converts data in vector format to pseudo-staSTATS format (i.e. structured
% cell array)
% x = ncond*nvar matrix
ncond = size(x,1);
nvar = size(x,2);
y = cell(1,nvar);
for ivar=1:nvar
    y{ivar}.means = x(:,ivar);
    y{ivar}.n = ones(ncond);
    y{ivar}.cov = eye(ncond);
    y{ivar}.regov = eye(ncond);
    y{ivar}.shrinkage=0;
    y{ivar}.weights = eye(ncond);
    y{ivar}.lm = eye(ncond);
    y{ivar}.shrink=0;
end
