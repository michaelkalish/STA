function [xStar, fStar, adjStar] = jCMRx1(y, model, E, mrTol)
import au.edu.adelaide.fxmr.model.*;

p = CMRxProblemMaker();

solver = CMRxSolver();
% Setting this will cause the solver to only return the first feasible solution that is found.
solver.setOnlyFeas(true);

if nargin > 3
    solver.setMrTolerance1(mrTol);
    solver.setMrTolerance2(mrTol * 1000);
end
if nargin < 3, E = {}; end
if nargin < 2, model=[]; end


nvar = numel(y); if isempty(model), model=ones(nvar,1); end % default model
ncond = numel(y{1});
% convert E to simpleConstraint objects
if iscell(E)
    if ~isempty(E)
        if iscell(E{1})
            %cells of cells
            for i=1:nvar
                Ecur = E{i};
                index = p.initAdj();
                for j=1:numel(Ecur)
                    p.addAdj(ncond, index, Ecur{j});
                end
            end
        else
            %Assume all adjs are same, from a cell
            index = p.initAdj();
            for j=1:numel(E)
                p.addAdj(ncond, index, E{j});
            end
            p.dupeAdj(nvar);
        end
    end
else
    if ndims(E)==1
        %2D matrix (repeat nvar times)
        p.addAdj(E);
        p.dupeAdj(nvar);
    else
        %3D Matrix
        for i=1:size(E,3)
            p.addAdj(E(:,:,i));
        end
    end
end

for ivar = 1:nvar
    p.addMeanArray(y{ivar}.means)
    w = y{ivar}.weights;
    if isempty(w)
        w = eye(numel(y{1}.means));
    end
    p.addWeightArray(w);
end

p.setModel(model);

xStar = cell(1,nvar);
adjStar = cell(1,nvar);
fStar = 0;
if numel(y{1}.means)==1 % deal with trivial case of single data point
    for i=1:numel(y)
        xStar{i}=y{i}.means;
    end
else
    problem = p.getProblem();
    solution = solver.solve(problem, false);
    fStar = solution.getFStar();
    for i=1:numel(y)
        xStar{i} = solution.getXStar(i-1); %i-1 converts to zero indexing
        adjStar{i} = solution.getAdjMatrix(i-1); %i-1 converts to zero indexing
    end
end