function staPLOT (data, varargin)
% function staPLOT (data, varargin)
% generates a state-trace plot
% generalization of stPLOT to n dependent variables
% data = data
% input arguments ('keyword',value):
%	predictions (prediction, pred, p) is a d-element cell array that defines a set of fitted values (output from staMR or staCMRx)
%	vars (v, dv) is a 2-vector that defines two dependent variables for plotting (default = [1, 2] )
%	groups (group, g) is a cell array that defines a structure on the conditions. For example if there are 2 between-subjects conditions and 3 within subjects conditions then one group structure might be {1:3, 4:6}. At present, only one factor can be identified in this way.
%	labels (label, lab, l) is a cell array that defines the labels of the levels defined by groups
%	axislabels (axislabel, axislab, axis) is a cell array that defines the labels of the x- and y-axes
%	axislimits (axislimit, axislim, limits, limit, lim) is a cell array that defines the upper and lower limits of the axes
%	axisticks (axistick, ticks, tick, t) is a cell array that defines tick marks for x- and y-axes (if specified then axislimits is redundant)
%	location (loc) is a Matlab value for the location of the legend (default = ‘northwest’)
%	color  (col, c)  is a 3-element RGB vector specifying the color of the plot of predicted values (default is black)
%	line is a flag to indicate if the predicted values are connected by a line. If set a line is drawn (default) otherwise no line is drawn.
%
% ***************************************************************************************************
% Modified 23 November 2016; 3 Jan 2017
% ***************************************************************************************************
%

pred= [];
vars = [1 2];
groups = [];
labels = [];
axislabels = [];
axislimits = [];
axisticks = [];
color = [];
lineflag = 1;
location = 'northwest';
shrink = -1;

for i = 1 : 2 : length(varargin)-1
    name = varargin{i}; 
    value = varargin{i+1};
    switch name
        case {'predictions', 'prediction', 'pred', 'p'}
            pred= value; % model predictions
        case {'vars', 'v', 'dv'} 
            vars = value; % 2-element vector identifying x and y DVs
        case {'groups', 'group', 'g'}
            groups = value; % cell array of main division into levels of one IV
        case {'labels', 'label', 'lab', 'l'}
            labels = value; % cell array corresponding labels for levels of IV
        case {'axislabels', 'axislabel', 'axislab', 'axis'}
            axislabels = value; % standard deviation of target distribution; if < 0 return estimate
        case {'axislimits', 'axislimit', 'axislim', 'limits', 'limit', 'lim'}
            axislimits = value; % cell array of xlim and ylim values
        case {'axisticks', 'axistick', 'ticks', 'tick', 't'}
            axisticks = value; % cell array of x and y ticks
        case {'location','loc'}
            location = value; % where to put the legend
        case {'color','col','c'}
            color = value; % colour of prediction plot as RGB vector
        case {'line'}
            lineflag = value; %0=plot points only, otherwise line (default)
        case {'shrink','s'}
            shrink = value; %shrinkage for staSTATS, default=-1 (estimated)
    end
end

if iscell(data)
    if isstruct(data{1})
        ys = data; % if structured then already in stats form
    else
        ys = staSTATS(data, shrink); % otherwise assume within-subjects data and get stats
    end
else
    celldata = gen2cell (data); % convert from old "general format"
    ys = staSTATS (celldata, shrink);
end

nvar = numel(ys);

if isempty(groups)
    groups = {1:numel(ys{vars(1)}.means)};
end

if isempty(labels)
    labels=cell(1,numel(groups));
    for i=1:numel(labels)
        labels{i} = ['Condition ' num2str(i)];
    end
end

if isempty(axislabels)
    axislabels = cell(1,nvar);
    for i=1:numel(axislabels)
        axislabels{i} = ['Outcome Variable ' num2str(i)];
    end
end

x = ys{vars(1)}.means; 
y = ys{vars(2)}.means; 

if numel(ys{vars(1)}.n) > 1
    cx = sqrt(diag(ys{vars(1)}.cov)./diag(ys{vars(1)}.n)); % between-subjects error bars
    cy = sqrt(diag(ys{vars(2)}.cov)./diag(ys{vars(2)}.n));
else
    cx = sqrt(diag(ys{vars(1)}.lm)/ys{vars(1)}.n); % within-subjects error bars
    cy = sqrt(diag(ys{vars(2)}.lm)/ys{vars(2)}.n);
end

% plot data
plotdata (x, y, groups, 0);

% plot error bars
errorbarx (x, y, cx, cy);
errorbary (x, y, cx, cy);

% plot model
if ~isempty(pred)
    xm = pred{vars(1)}; ym = pred{vars(2)};
    ix = tiesort(xm,ym);
    if isempty(color)
        if lineflag==0
            plot (xm(ix), ym(ix), 'kx');
        else
            plot (xm(ix), ym(ix), 'k--'); 
            
        end
    else
        if lineflag==0
            plot (xm(ix), ym(ix), 'x', 'Color', color);
        else
            plot (xm(ix), ym(ix)', 'Color',  color);
        end
    end
end

% re-plot data over error bars
plotdata (x, y, groups, 1);
hold off

% set ticks
if ~isempty(axisticks)
    set (gca, 'XTick', axisticks{vars(1)});
    set (gca, 'YTick', axisticks{vars(2)});
    if isempty(axislimits) % if limits unspecified, infer from ticks (if specified)
        axislimits{1} = [axisticks{vars(1)}(1), axisticks{vars(1)}(end)];
        axislimits{2} = [axisticks{vars(2)}(1), axisticks{vars(2)}(end)];
    end
end

% set axis limits
if ~isempty(axislimits)
    xlim(axislimits{vars(1)}); ylim(axislimits{vars(2)});
else
    xlim('auto'); ylim('auto');
end

% output axis labels
a = 'xlabel('; a = [a '''' axislabels{vars(1)} '''' ');']; eval(a);
a = 'ylabel('; a = [a '''' axislabels{vars(2)} '''' ');']; eval(a);

% output legend
a = 'legend(';
for i=1:numel(labels)
    a = [a, '''' labels{i} '''' ','];
end
if ~isempty(pred)
    a = [a, '''Predictions'','];
end
a = [a, '''location''', ',',  '''' location '''', ');']; 
eval(a);

function plotdata(x,y,groups,flag)
msize=8; colour = [.7 .7 .7];
for igroup = 1:numel(groups)
    a = x(groups{igroup});
    b = y(groups{igroup});
    if igroup==1
        h=plot (a, b, 'ko', 'markerfacecolor', colour, 'markersize', msize); hold on
    elseif igroup==2
        h=plot (a, b, 'ko', 'markerfacecolor', 'w', 'markersize', msize);
    elseif igroup==3
        h=plot (a, b, 'k^', 'markerfacecolor', colour, 'markersize', msize);
    elseif igroup==4
        h=plot (a, b, 'k^', 'markerfacecolor', 'w', 'markersize', msize);
    elseif igroup==5
        h=plot (a, b, 'ks', 'markerfacecolor', colour, 'markersize', msize);
    else 
        h=plot (a, b, 'ks', 'markerfacecolor', 'w', 'markersize', msize);
    end
    if flag
        DeleteLegendEntry (h);
    end 
end

function errorbarx (x, y, cx, cy)
% plots errorbars in x direction
yrange=max(y)-min(y)+3*max(cy); 
ticklength=yrange/40;
for i=1:numel(x)
    a(1) = x(i)-cx(i);
    a(2) = x(i)+cx(i);
    b(1) = y(i); 
    b(2) = y(i);
    h = plot (a, b, 'k');
    DeleteLegendEntry (h);
% plot ticks
    for j=1:2
        aa(1)=a(j); 
        aa(2)=a(j);
        bb(1)= b(j)-ticklength/2;
        bb(2)= b(j)+ticklength/2;
        h = plot (aa, bb, 'k');
        DeleteLegendEntry (h);
    end
end

function errorbary (x, y, cx, cy)
% plots error bars in y direction
xrange=max(x)-min(x)+3*max(cx); 
ticklength=xrange/40;
for i=1:numel(x)
    a(1) = x(i);
    a(2) = x(i);
    b(1) = y(i)-cy(i); 
    b(2) = y(i)+cy(i);
    h = plot (a, b, 'k');
    DeleteLegendEntry (h);
% plot ticks
    for j=1:2
        bb(1)=b(j); 
        bb(2)=b(j);
        aa(1)= a(j)-ticklength/2;
        aa(2)= a(j)+ticklength/2;
        h = plot (aa, bb, 'k');
        DeleteLegendEntry (h);
    end
end

function DeleteLegendEntry (h)
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

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

