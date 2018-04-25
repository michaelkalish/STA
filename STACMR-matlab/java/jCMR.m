function [xStar, fStar, nitr, FL, FU] = jCMR(y, E)
% function [xStar, fStar, nitr] = jCMR(y, E)
% CMR that calls Java for calculations
%
import au.edu.adelaide.fxmr.model.CMRSolver;
import au.edu.adelaide.fxmr.model.CMRProblemMaker;

p = CMRProblemMaker();
solver = CMRSolver();

if nargin==1
    E={};
end
if iscell(E)
    for i=1:numel(E)
        p.addRangeSet(E{i});
    end
else
    %
    disp('Currently only cell structure allowed')
    return
end

nvar = numel(y);
for ivar = 1:nvar
    p.addMeanArray(y{ivar}.means)
end

for ivar = 1:nvar
    w = y{ivar}.weights;
    if isempty(w)
        w = eye(numel(y{1}.means)); % extract weights
    end
    p.addWeightArray(w);
end

xStar = cell(1,numel(y)); fStar = 0;
if numel(y{1}.means)==1 % deal with trivial case of single data point
    for i=1:numel(y)
        xStar{i}=y{i}.means;
    end
else
    problem = p.getProblem();
    solution = solver.solve(problem);
    fStar = solution.getFStar();
    xStar = solution.getXStar();
    nitr = solution.getIters().size();
end