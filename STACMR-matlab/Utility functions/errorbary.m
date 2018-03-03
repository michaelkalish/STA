function errorbary (x, y, y_low, y_high)
% function errorbary (x, y, y_low, y_high)
% plots error bars in y direction
if nargin==3, y_high=y_low; end
for i=1:numel(x)
    a(1) = x(i);
    a(2) = x(i);
    b(1) = y(i)-abs(y_low(i)); 
    b(2) = y(i)+abs(y_high(i));
    h = plot (a, b, 'k');
    DeleteLegendEntry (h);
end
% plot ticks
xx=xlim; ticklength=(xx(2)-xx(1))/25;
for i=1:numel(x)
    a(1) = x(i);
    a(2) = x(i);
    b(1) = y(i)-abs(y_low(i)); 
    b(2) = y(i)+abs(y_high(i));
    for j=1:2
        bb(1)=b(j); 
        bb(2)=b(j);
        aa(1)= a(j)-ticklength/2;
        aa(2)= a(j)+ticklength/2;
        h = plot (aa, bb, 'k');
        DeleteLegendEntry (h);
    end
end