function c = graphtransclosure (g, i, k)
% function c = graphtransclosure (g, i, k)
% returns transitive closure of directed acyclic graph g
% g is adjacency matrix (weights zero or one)
% implements Warshall's algorithm
c = g;
if nargin > 1 % update specified node
    if c(i,k) > 0
        c(i,:) = c(i,:) + c(k,:);
    end
else % update all nodes
    n = size(c,1);
    for k = 1:n
        for i = 1:n
            if c(i,k) > 0
                c(i,:) = c(i,:) + c(k,:);
            end
        end
    end
end
c(c > 1) = 1;
