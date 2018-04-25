function r = resampleBN (y)
% function r = resampleBN (y)
% y is output from binSTATS or binSTATS2
r = y;
for i=1:size(y,1)
    for j=1:size(y,2)
        if isstruct(y{i,j})
            b = binornd (y{i,j}.n, y{i,j}.means);
            r{i,j}.count = [b r{i,j}.n - b];
            r{i,j}.means = b./r{i,j}.n;
        else
            n = sum(y{i,j},2);
            b = binornd (n, y{i,j}(:,1)./n);
            r{i,j} = [b n-b];
        end
    end
end
