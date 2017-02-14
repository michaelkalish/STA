function [x, fit] = jMR(y, W, E, mrTol)
% function [x, fit] = jMR(y, W, E)
import au.edu.adelaide.fxmr.model.mr.MRSolverAJOptimiser;
import au.edu.adelaide.fxmr.model.mr.MRProblemMaker;

p = MRProblemMaker(y);
solver = MRSolverAJOptimiser();

nvar = numel(y);
if nargin==1
    W = eye(nvar);
    E={};
elseif nargin==2
    E={};
elseif nargin==4
    %Do not use tol < 1e-11
    solver.setTolerance(mrTol, mrTol * 1000);
end

p.setWeightArray(W);


% convert E to adjacency matrix (null if not specified) - this code will
% need some work in the future - the 3d matrix isn't used explicitly
if iscell(E)
    %This could be improved...
    adj = cell2adj2(1:nvar, E);
else
    adj = E;
end

%Add constraints one by one (this should be done for the jCMR and jCMRx)
[pos,neg]=find(adj);
n = numel(pos);
for i=1:n
   p.addConstraint(pos(i),neg(i));
end

problem = p.getProblem();
solution = solver.solve(problem);
fit = solution.getfVal();
x = solution.getxVector();
