function [p, datafit, fits, pars] = jCMRBNfits (nsample, data, E, proc, seed)
import au.edu.adelaide.fxmr.model.bin.BinCMRFits;
import au.edu.adelaide.fxmr.model.bin.BinCMRProblemMaker;

if ~iscell(data)
    disp('Currently only cell structure allowed for data')
    return
end

nSubj = size(data,1);
nVar = size(data,2);

pm = BinCMRProblemMaker(nSubj, nVar);
if nargin<=2
    E={};
end
if nargin<=3
    proc=-1;
end
if nargin<=4   
    seed=-1;
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

fObj = BinCMRFits(nsample, problem, proc, false, seed);

p = fObj.getP();

datafit = fObj.getBaseFitDiff();
fits = fObj.getFits()';
pars = fObj.getXStars();
