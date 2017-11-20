function staPLOTBN (data, varargin)
% function staPLOTBN (data, varargin)
% generates a state-trace plot for binomial data
% generalization of stPLOT to n dependent variables
% data = data
% pred = predicted model estimates (cell array output from staCMR)
% vars = 2-element vector identifying x and y DVs (if empty then
% first two DVs are used)
% groups = cell array of main division into levels of one IV
% labels = cell array corresponding labels for levels of IV
% axislabels = cell array of axis labels
% axislimits = cell array of xlim and ylim values
%
% For example, suppose we have one IV with two levels (labelled 'High' and 'Low') and another
% IV with three levels (which we could label '1','2','3', but this will not feature on the plot.
% There are thus 6 conditions organised so that the first 3 are 'High' and the next 3 are 'Low'.
% Suppose too that we have previously run staCMR with the call, [x f] = staCMR (data, E), where 
% E is an optional partial order.
% Finally, suppose the 2 DVs can be labelled, 'Mean RT' for DV 1 and 'Mean Error' for DV 2 and that
% plausible ranges for DVs 1 and 2 are 300-1000 and 0.5-1.0, respectively.
% Then the appropriate call to staPLOT is:
% staPLOT (data, x, {[1 2 3] [4 5 6]}, {'High' 'Low'}, {'Mean RT' 'Mean Error'}, {[300 1000] [.5 1]};
% ***************************************************************************************************
% Modified 10 November 2015 (and again in April 2017)
% ***************************************************************************************************
%
if ~iscell(data)
    ys = BNgen2cell(data); ys = binSTATS (ys); 
elseif isstruct(data{1})
    ys = data; % if structured then already in stats form
else
    ys = binSTATS(data);
end

model = [];
vars = [1 2];
subj = 1;
groups = [];
labels = [];
axislabels = [];
axislimits = {};

for i = 1 : 2 : length(varargin)
    name = varargin{i};
    value = varargin{i+1};
    switch name
        case {'predicted','pred', 'p'}
            model = value; % model estimates
        case {'vars', 'v', 'dv'} 
            vars = value; % 2-element vector identifying x and y DVs
        case {'subj', 's'} 
            subj = value; % subject number
        case {'groups', 'g'}
            groups = value; % cell array of main division into levels of one IV
        case {'labels', 'label', 'lab', 'l'}
            labels = value; % cell array corresponding labels for levels of IV
        case {'axislabels', 'axis', 'alab'}
            axislabels = value; % standard deviation of target distribution; if < 0 return estimate
        case {'axislimits', 'limits', 'lim'}
            axislimits = value; % cell array of xlim and ylim values
    end
end

if isempty(groups)
    groups = {1:numel(ys{vars(1)}.means)};
    nogroups = 1;
else
    nogroups = 0;
end

if isempty(labels)
    labels=cell(1,numel(groups));
    for i=1:numel(labels)
        labels{i} = ['Condition ' num2str(i)];
    end
end

if isempty(axislabels)
    axislabels = cell(1,2);
    for i=1:numel(axislabels)
        axislabels{i} = ['Outcome Variable ' num2str(i)];
    end
end

x = ys{subj,vars(1)}.means; 
y = ys{subj,vars(2)}.means; 

cx = sqrt(x.*(1-x)./ys{subj,vars(1)}.n);
cy = sqrt(y.*(1-y)./ys{subj,vars(2)}.n);

% plot data
plotdata (x, y, groups, 0);

% plot error bars
errorbarx (x, y, cx, cy);
errorbary (x, y, cx, cy);

% plot model
if ~isempty(model)
    xm = model{subj,vars(1)}; ym = model{subj,vars(2)};
    ix = tiesort(xm,ym);
    plot (xm(ix), ym(ix), 'r--'); 
end

% re-plot data over error bars
plotdata (x, y, groups, 1);
hold off

% set axis limits
if ~isempty(axislimits)
    xlim(axislimits{1}); ylim(axislimits{2});
else
    xlim('auto'); ylim('auto');
end

% output axis labels
a = 'xlabel('; a = [a '''' axislabels{1} '''' ');']; eval(a);
a = 'ylabel('; a = [a '''' axislabels{2} '''' ');']; eval(a);

% output legend
if nogroups==0
    a = 'legend(';
    for i=1:numel(labels)
        a=[a '''' labels{i} '''' ','];
    end
    a=[a '''location'', ''northwest'');']; eval(a);
end

function plotdata(x, y, groups, flag)
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

