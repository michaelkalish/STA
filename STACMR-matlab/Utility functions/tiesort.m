function ix = tiesort (xx, yy)
% sorts y values in increasing x-order
% within blocks of tied x-values, sorts y values in increasing y-order
bignumber=100;
x=round(xx*bignumber)/bignumber;
y=round(yy*bignumber)/bignumber;
t=1:numel(x);
z=[x y t'];
a=sortrows(z,[1 2]);
ix=a(:,3);
