function [p, datafit, fits, maxbad, times] = jCMRfitsx(nsample,y, model, E, shrink, proc, cheapP, mrTol)
% function [xStar, fStar, nitr] = jCMR(y, E)
% CMR that calls Java for calculations
%
import au.edu.adelaide.fxmr.model.*;

if nargin<4
    E={};
end
if nargin<5
    shrink = -1;
end
if nargin < 6
    proc=0;
end
if nargin < 7
    mrTol=0;
end

if iscell(y)
    if isstruct(y{1})
        problem = CMRxFitsProblemMaker();
        nvar = numel(y);
        ncond = numel(y{1});
    else
        problem = CMRxFitsGMProblemMaker();
        problem.setShrink(shrink);
        ngroup = size(y,1); nvar = size(y,2);
        ncond = size(y{1,1},2)*ngroup;
    end
else
    problem = CMRxFitsGMProblemMaker();
    problem.setShrink(shrink);
    ngroup = numel(unique(y(:,2)));
    nvar = numel(unique(y(:,3)));
    ncond = ngroup * (size(y,2) - 3);
end

if nargin<3
    model=ones(ncond,1);
end


% convert E to simpleConstraint objects
if iscell(E)
    if ~isempty(E)
        if iscell(E{1})
            %cells of cells
            for i=1:nvar
                Ecur = E{i};
                index = problem.initAdj();
                for j=1:numel(Ecur)
                    problem.addAdj(ncond, index, Ecur{j});
                end
            end
        else
            %Assume all adjs are same, from a cell
            index = problem.initAdj();
            for j=1:numel(E)
                problem.addAdj(ncond, index, E{j});
            end
            problem.dupeAdj(nvar);
        end
    end
else
    if ndims(E)==1
        %2D matrix (repeat nvar times)
        problem.addAdj(E);
        problem.dupeAdj(nvar);
    else
        %3D Matrix
        for i=1:size(E,3)
            problem.addAdj(E(:,:,i));
        end
    end
end

problem.setModel(model);

if iscell(y)
    if isstruct(y{1})
        ns = zeros(1,nvar);
        for ivar = 1:nvar
            problem.addMeanArray(y{ivar}.means);
            problem.addCov(y{ivar}.cov);
            w = y{ivar}.weights;
            if isempty(w)
                w = eye(numel(y{ivar}.means));
            end
            problem.addWeightArray(w);
            if size(y{ivar}.n)==1
                ns(ivar) = y{ivar}.n;
            else
                ns(ivar) = y{ivar}.n(1,1);
            end
        end
        problem.setN(ns);
        
        sol = CMRxFits(nsample,problem.getProblem(),shrink,proc,cheapP,mrTol,mrTol*1000);
    else
        for group=1:ngroup
            for var=1:nvar
                problem.addCell(group,var,y{group,var})
            end
        end
        sol = problem.solve(nsample,proc,cheapP,mrTol,mrTol*1000);
    end
else
    problem.setGM(y);
    sol = problem.solve(nsample,proc,cheapP,mrTol,mrTol*1000);
end

p = sol.getP();
fits = sol.getFits();
datafit = sol.getDataFit();
maxbad = zeros(1,nsample);
bads = sol.getBadnesses();
for i=1:nsample
    maxbad(i) = bads(i).ordinal();
end
times=sol.getTimes();