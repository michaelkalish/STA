function G2 = g2BN (pred, data)
% function G2 = g2BN (pred, data)
% returns vector of G-squared from binary data
% pred is predicted means (e.g., output from staCMRBN)
%
if ~iscell(data)
    y = BNgen2cell(data); y = binSTATS(y);
else
    if isstruct(data{1})
        y = data;
    else
        y = binSTATS(data);
    end
end
G2 = zeros(size(y));
for isub = 1:size(y,1)
    for ivar = 1:size(y,2)
        n = y{isub,ivar}.n; x = pred{isub,ivar}; c = y{isub,ivar}.count;
        k=find(x <= 0);
        if ~isempty(k)
            x(k)=.5/n(k); % replace expected zeros with small number
        end
        k=find(x >= 1);
        if ~isempty(k)
            x(k)=(n(k)-.5)/nn(k); % replace expected ones with number close to 1
        end
        a = [x 1-x]; 
        a = a.*repmat(n,1,size(a,2)); % expected counts
        g2 = 2*c.*log(c./a); g2(n<=0) = 0; % g-squared
        f = sum(sum(g2)); f(f<0)=0;
        G2(isub,ivar) = f;
    end
end

