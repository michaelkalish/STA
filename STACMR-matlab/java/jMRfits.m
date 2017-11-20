function [p, datafit, fits, maxbad, times] = jMRfits(nsample, y, E, shrink, reverse, proc, mrTol, seed, showStatus)
% function [p, datafit, fits, maxbad, times] = jMRfits(nsample, y, E, shrink, proc)
% jMRfits that calls the same java code as jCMRfits but with a parameter to
% use the non-coupled code (model is ignored and never set).

import au.edu.adelaide.fxmr.model.*;

if nargin<3, E={}; end
if nargin<4, shrink = -1; end
if nargin<5, reverse=false; end
if nargin<6, proc=0; end
if nargin<7, mrTol=0; end
if nargin<8, seed=-1; end
if nargin<9, showStatus=false; end

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


if iscell(y)
    if isstruct(y{1})
        ns = zeros(1,nvar);
        for ivar = 1:nvar
            problem.addMeanArray(y{ivar}.means);
            problem.addCov(y{ivar}.regcov);
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
        sol = CMRxFits(nsample,problem.getProblem(),shrink,proc,false,true,mrTol,mrTol*1000,false,reverse, seed, showStatus);
    else
        for group=1:ngroup
            for var=1:nvar
                problem.addCell(group,var,y{group,var})
            end
        end
        sol = problem.solve(nsample,proc,false,true,mrTol,mrTol*1000,false,reverse,seed,showStatus);
    end
else
    problem.setGM(y);
    sol = problem.solve(nsample,proc,false,true,mrTol,mrTol*1000,false,reverse,seed,showStatus);
end

p = sol.getP();
fits = sol.getFits();
datafit = sol.getDataFit();
maxbad = zeros(1,nsample);
bads = sol.getBadnesses();
for i=1:bads.length
    maxbad(i) = bads(i).ordinal();
end
times=sol.getTimes();