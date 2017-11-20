function [p, datafit, fits] = jCMRxBNfits (nsample, data, E, model, approximate, showStatus)
import au.edu.adelaide.fxmr.model.bin.BinCMRxFits;
import au.edu.adelaide.fxmr.model.bin.BinCMRProblemMaker;

if ~iscell(data)
    disp('Currently only cell structure allowed for data')
    return
end

nSubj = size(data,1);
nVar = size(data,2);

pm = BinCMRProblemMaker(nSubj, nVar);
if nargin<3
    E={};
end
if nargin<4 || size(model,1) ~= nVar
    model=ones(nVar,1);
end
if nargin<5
    approximate=false;
end
if nargin<6
    showStatus=false;
end
if iscell(E)
    for i=1:numel(E)
        pm.addRangeSet(E{i});
    end
else
    disp('Currently only cell structure allowed for E')
    return
end

pm.setModel(model);

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
proc = 0;

fObj = BinCMRxFits(nsample, problem, proc, approximate, showStatus);

p = fObj.getP();

datafit = fObj.getBaseFitDiff();
fits = fObj.getFits()';