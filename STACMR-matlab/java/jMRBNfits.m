function [p, datafit, fits, pars] = jMRBNfits (nsample, data, E)
import au.edu.adelaide.fxmr.model.bin.BinCMRFits;
import au.edu.adelaide.fxmr.model.bin.BinCMRProblemMaker;

if ~iscell(data)
    disp('Currently only cell structure allowed for data')
    return
end

nSubj = size(data,1);
nVar = size(data,2);
proc = -1;

pm = BinCMRProblemMaker(nSubj, nVar);
if nargin==2
    E={};
end
if iscell(E)
    for i=1:numel(E)
        pm.addRangeSet(E{i});
    end
else
    %
    disp('Currently only cell structure allowed for E')
    return
end

for s=1:nSubj
    for v=1:nVar
        if isstruct(data{s,v})
            %Stats version
            pm.setElement(s-1,v-1,data{s,v}.count)
        else
            pm.setElement(s-1,v-1,data{s,v})
        end
    end
end

problem = pm.getBaseProblem();

fObj = BinCMRFits(nsample, problem, proc, true);

p = fObj.getP();

datafit = fObj.getBaseFitDiff();
fits = fObj.getFits()';
pars = fObj.getXStars();
