function yb = bootstrap (y, type)
% function yb = bootstrap (y)
% draws bootstrap sample from data
%
% type=0 is cell array format
% type>0 is general format
% type<0 is parametric bootstrap (data only in summary form)
%        assumes y is output from staSTATS or outSTATS
%        assumes normal distribution
if nargin==1
    if iscell(y)
        if isstruct(y{1})
            type = -1; 
        else
            type = 0;
        end
    else
        type = 1;
    end
end
yb = cell(size(y));
if type == 0 % y = cell array of NSUBJ x NCOND matrices
    for ivar = 1:size(y,2)
        for igroup = 1:size(y,1)
            a = y{igroup,ivar};
            nsub = size(a,1);
            yy = zeros(size(a)); b = zeros(1,nsub);
            v=nsub;
            while v==nsub
                for isub = 1:nsub
                    u = randperm(nsub);
                    yy(isub,:) = a(u(1),:);
                    b(isub)=u(1);
                end
                v = sum(b==b(1)); % check if only sampled a single subject
            end
            yb{igroup,ivar}=yy;
        end
    end
elseif type == 1 % y = "general format"
    cond = unique (y(:,2));
    var = unique (y(:,3));
    yb = y;
    for j=1:numel(var)
        for i=1:numel(cond)
            k = find(y(:,2)==cond(i) & y(:,3)==var(j));
            a = y(k,:);
            r = floor(rand(1,numel(k))*numel(k))+1;
            yb(k,:) = a(r,:);
        end
    end
else % parametric bootstrap (no raw data)
    yb = cell(1,numel(y));
    for ivar = 1:numel(y)
        if size(y{ivar}.n)==1
            nsub = y{ivar}.n;
        else
            nsub = y{ivar}.n(1,1);
        end
        yb{ivar} = mvnrnd(y{ivar}.means, y{ivar}.cov, nsub);
    end
end
