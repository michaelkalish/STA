function errorbarx (x, y, x_low, x_high)
% function errorbarx (x, y, x_low, x_high)
% plots errorbars in x direction
if nargin==3, x_high=x_low; end
for i=1:numel(x)
    a(1) = x(i)-abs(x_low(i));
    a(2) = x(i)+abs(x_high(i));
    b(1) = y(i); 
    b(2) = y(i);
    h = plot (a, b, 'k');
    DeleteLegendEntry (h);
end
yy=ylim; ticklength=(yy(2)-yy(1))/25;
% plot ticks
for i=1:numel(x)
    a(1) = x(i)-abs(x_low(i));
    a(2) = x(i)+abs(x_high(i));
    b(1) = y(i); 
    b(2) = y(i);
    for j=1:2
        aa(1)=a(j); 
        aa(2)=a(j);
        bb(1)= b(j)-ticklength/2;
        bb(2)= b(j)+ticklength/2;
        h = plot (aa, bb, 'k');
        DeleteLegendEntry (h);
    end
end
