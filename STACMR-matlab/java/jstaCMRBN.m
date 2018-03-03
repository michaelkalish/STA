function [x, G2, f, nitr] = jstaCMRBN(data, E, model, approximate)
% function [x, G2, f, nitr] = staCMRBN (data, E)
% data is an NSUB x NVAR cell array of binomial data
% where each element is an NCOND x 2 matrix (hits, misses)
import au.edu.adelaide.fxmr.model.bin.BinCMRxSolver;
import au.edu.adelaide.fxmr.model.bin.BinCMRProblemMaker;

if ~iscell(data)
    disp('Currently only cell structure allowed for data')
    return
end

nSubj = size(data,1);
nVar = size(data,2);

p = BinCMRProblemMaker(nSubj, nVar);
if nargin<2
    E={};
end
if nargin<3 || size(model,1) ~= nVar
    model=ones(nVar,1);
end
if nargin<4
    approximate=false;
end

if iscell(E)
    for i=1:numel(E)
        p.addRangeSet(E{i});
    end
else
    %
    disp('Currently only cell structure allowed for E')
    return
end

p.setModel(model);

for s=1:nSubj
    for v=1:nVar
        %Minus 1 to zero index
        if isstruct(data{s,v})
            %Stats version
            p.setElement(s-1,v-1,data{s,v}.count)
        else
            p.setElement(s-1,v-1,data{s,v})
        end
    end
end

solver = BinCMRxSolver();
solver.setOnlyFeas(approximate);

solutions = solver.solve(p.getProblems());

x = cell(nSubj,nVar); f = zeros(nSubj,1); G2 = f; nitr=f;
for isub = 1:nSubj
    f(isub) = solutions(isub).getFStar();
    G2(isub) = solutions(isub).getG2Star();
    nitr(isub) = solutions(isub).getIter();
    for ivar = 1:nVar
        x{isub,ivar} = solutions(isub).getXStar(ivar-1);
    end
end
