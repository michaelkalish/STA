function output = staSTATS (data, shrink, warning)
% function output = staSTATS (data, shrink, warning)
% returns summary statistics on data
% data is a ngroup * nvar cell array of matrices:
% each row of data corresponds to one of ngroup groups of participants
% each column of data corresponds to one of nvar dependent variables (DVs)
% each cell of data is a nsub_ij x nwithin_ij matrix of observations
% where, nsub_ij = no. of participants in group i for variable j,
% nwithin_ij = no. of within-subjects conditions in group i for variable j 
% This means that the total number of conditions (ncond) is sum_i=1:N (nwithin_i1)
% Note: The same participants may or may not be measured under each DV.
% This is currently ignored. Each DV is assumed to be independent of all
% other DVs. Also assumed that each DV is measured under all conditions
% (otherwise STA/SDA is meaningless)
% shrink is option for covariance matrix shrinkage (0=none, 1=diag,
% -1=estimated)
% if warning is set then a warning message is printed if nans are detected
%
% staSTATS returns a 1 * nvar cell array means, cov, nsub and weights
% output{ivar}.means = 1 x ncond vector of observed means for DV ivar
% output{ivar}.n = ncond x ncond matrix of number of subjects
% output{ivar}.cov = ncond x ncond observed covariance matrix
% output{ivar}.regcov = regularized covariance matrix
% output{ivar}.shrinkage = shrinkage parameter for regularization
% output{ivar}.weights = 1 x ncond weight matrix for monotonic regression
% output{ivar}.lm = ncond x ncond matrix of loftus-masson within-subject variance
% output{ivar}.shrink = shrinkage factor (between 0 and 1)
% if shrink==1 then diagonalise weight matrix, if shrink==0 then do nothing
%
% *************************************************************************
% last modified 13 October 2016 (to treat nans as missing)
% *************************************************************************

if nargin < 3, warning = 0; end
if nargin < 2, shrink = -1; end %  estimate shrinkage by default

if ~iscell(data)
    y = gen2cell(data); % convert from general format
else
    y = data; % cell array format
end

ngroup = size(y,1); nvar = size(y,2);

% calculate statistics
output = cell(1,nvar); 
for ivar = 1:nvar
    output{ivar}.means = [];
    output{ivar}.cov = [];
    output{ivar}.regcov = [];
    output{ivar}.n = [];
    output{ivar}.weights = [];
    output{ivar}.lm = [];
    output{ivar}.shrinkage = [];
    output{ivar}.nanflag = [];
    output{ivar}.bad = [];
    for igroup = 1:ngroup
        yy = y{igroup, ivar};
%        u = sum(isnan(yy),2); yy=yy(u==0,:); % delete nans
        out.nanflag = sum(sum(isnan(yy))); % check for nans
        out.means = mean(yy,1,'omitnan'); % means excluding nans
        out.means(isnan(out.means)) = 0; % Means can not be NAN (means whole column is NAN & bad will be 2 later)
        a=yy;a(isnan(yy))=0;a(~isnan(yy))=1; out.n=a'*a; % number of observations
        out.cov = cov(yy,1,'partialrows'); out.cov(isnan(out.cov))=0;  % covariance
        eigCov = eig(out.cov);
        if cond(out.cov) < 1e6 && min(eigCov) > 0
            [out.regcov, out.shrinkage] = shrinkDiag(yy, shrink); % shrink if required
            out.bad = 0;
        else
            [out.regcov, out.shrinkage] = shrinkDiag(yy, 1); % diagonalise bad covariance matrix
            out.bad = 1;
        end
        eigRegCov = eig(out.regcov);
        if cond(out.regcov) > 1e6 || min(eigRegCov) <= 0
            out.regcov = trace(out.regcov)*eye(size(out.regcov,1))/size(out.regcov,1); % make sure all diagonal elements are nonzero
            out.bad = 2;
        end
        outNmin1 = max(out.n,eye(size(out.regcov,1)));
        out.weights =outNmin1.*out.regcov^-1; % weights
        out.lm = loftusmasson(yy); % loftusmasson within subjects variance
% add to output
        output{ivar}.means = [output{ivar}.means, out.means];
        output{ivar}.cov = blkdiag(output{ivar}.cov, out.cov);
        output{ivar}.regcov = blkdiag(output{ivar}.regcov, out.regcov);
        output{ivar}.n = blkdiag(output{ivar}.n, out.n);
        output{ivar}.weights = blkdiag(output{ivar}.weights, out.weights);
        output{ivar}.lm = blkdiag(output{ivar}.lm, out.lm);
        output{ivar}.shrinkage = [output{ivar}.shrinkage, out.shrinkage];
        output{ivar}.nanflag = [output{ivar}.nanflag, out.nanflag];
        output{ivar}.bad = [output{ivar}.bad, out.bad];
    end
    if sum(output{ivar}.nanflag) > 0 && warning ~= 0
        fprintf ('Warning: %d NANs detected for variable %d', sum(output{ivar}.nanflag), ivar);
        fprintf ('. These are treated as missing.\n');
    end
    if output{ivar}.bad > 0 & warning ~= 0
        fprintf ('Bad covariance matrix detected for variable %d. Type = %d\n', ivar, output{ivar}.bad);
    end
end

function [sigma,shrinkage]=shrinkDiag(x,shrink)
% function [sigma,shrinkage]=shrinkDiag(x,shrink)
% x (t*n): t iid observations on n random variables
% sigma (n*n): invertible covariance matrix estimator
%
% Shrinks towards diagonal matrix
% if shrink is specified, then this constant is used for shrinkage

% de-mean returns
t=size(x,1);
meanx=mean(x,1,'omitnan');

% compute sample covariance matrix
sample = cov(x,1,'partialrows'); sample(isnan(sample))=0;

% compute prior
prior=diag(diag(sample)); 

if (nargin < 2 || shrink == -1) % compute shrinkage parameters
  
  % what we call p 
  y=(x-meanx(ones(t,1),:)).^2; y(isnan(y))=0;
  phiMat = y'*y/t-sample.^2; 
  phi=sum(sum(phiMat)); 
  
  % what we call r
  rho=sum(diag(phiMat)); 
  
  % what we call c
  gamma=norm(sample-prior,'fro')^2; 

  % compute shrinkage constant
  kappa=(phi-rho)/gamma;
  shrinkage=max(0,min(1,kappa/t));
    
else % use specified constant
  shrinkage=shrink;
end

% compute shrinkage estimator
sigma=shrinkage*prior+(1-shrinkage)*sample;

function r = loftusmasson(yy)
% returns Loftus & Masson within subjects error based on Mean Square residual
y_cond = repmat(mean(yy,1,'omitnan'), size(yy,1), 1); 
y_subj = repmat(mean(yy,2,'omitnan'), 1, size(yy,2)); 
m = mean(reshape(yy,1,numel(yy)),'omitnan'); 
y_mean = repmat(m, size(yy,1), size(yy,2)); 
ya = yy - y_cond - y_subj + y_mean; ya(isnan(ya))=0; 
ss = sum(sum(ya.*ya));
df = (size(yy,1)-1)*(size(yy,2)-1);
ms_resid = ss/df;
r = diag(repmat(ms_resid,size(yy,2),1));
