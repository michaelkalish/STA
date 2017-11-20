function [p_agg, datafit_agg, fits_agg] = aggBN (nsamples, datafit, fits)
% function [p_agg, datafit_agg, fits_agg] = aggBN (nsamples, datafit, fits)
% aggregates across subject-level BN fits
% fits is nsub x ns matrix of fits (from CMRBNfits)
datafit_agg = sum(datafit);
fits_agg = zeros(1,nsamples);
nsub = size(fits,1); ns=size(fits,2);
for isample=1:nsamples
    j = floor(rand(nsub,1)*ns)+1;
    for i=1:nsub
        fits_agg(isample) = fits_agg(isample) + fits(i,j(i));
    end
end
k = find(datafit_agg <= fits_agg);
p_agg=numel(k)/nsamples;
