function A = adj2ineq (adj)
% function A = adj2ineq (adj, d2)
% converts an adjacency matrix to a inequality matrix (A) for MR
% Last updated 28 October 2014
%
A = []; 
[i, j] = find(adj); % find indices of non-zero elements of adj
if ~isempty(i)
    n = numel(i);
    A = zeros(n, size(adj,1)); % pre-allocate A
    t = 1:n; t = t';
    k = sub2ind(size(A),t,i); A(k) = 1; % set row index to 1 
    k = sub2ind(size(A),t,j); A(k) = -1; % set column index to -1
end


