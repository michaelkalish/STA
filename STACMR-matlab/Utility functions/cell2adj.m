function adj = cell2adj (nodes, E)
% function adj = cell2adj  (nodes, E)
% converts a partial order model in cell array form to an adjacency matrix suitable for MR
% nodes = set of elements (usually 1:n)
% E = partial order on nodes
% returns corresponding adjacency matrix
% ****************************************************************************************
% modified by Luke Finlay 1 March 2016
% ****************************************************************************************
%
if nargin==1
    E={};
end
if ~iscell(E)
    E={E};
end
n=numel(nodes);
adj=zeros(n,n);
if n > 1
    if ~isempty (E)
        for i=1:numel(E)
            if ~isempty(E{i})
                if numel(E{i}) >= 2
                    u = E{i};
                    for j=1:(length(u) - 1)
                        adj(u(j),u(j+1)) = 1;
                    end
                end
            end
        end
    end
end
