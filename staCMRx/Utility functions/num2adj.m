function adj = num2adj(x)
% function adj = num2adj(x)
% converts a vector of numbers into an adjacency matrix
if isrow(x)
    x = x'; % make sure x is a column vector
end
u = repmat(x,1,numel(x));
adj = u - u'; % difference matrix
adj(adj > 0) = 0; adj(adj < 0) = 1; % convert to adjacency
