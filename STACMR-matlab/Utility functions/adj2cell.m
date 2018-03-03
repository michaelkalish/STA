function E = adj2cell(adj)
% converts adjacency matrix to cell array
[row, col] = find(adj==1);
E=cell(1,numel(row));
for i=1:numel(row)
    E{i} = [row(i) col(i)];
end
