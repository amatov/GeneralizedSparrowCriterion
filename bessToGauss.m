function [res,jac]=bessToGauss(sigma_ga,hB,height)

% bessToGauss prepares the parameters for the 
% non-linear fitting
%
% SYNOPSIS [res,jac]=bessToGauss(sigma_ga,hB,height)
%
% INPUT    sigma_ga  :  sigma of the GK
%          hB        :  the PSF
%          height    :  height of the GK
%
% OUTPUT   res       : retrun residual
%          jac       : jacobian
%
% Alexandre Matov, January 7th, 2003

sZ=size(hB,1); %size of the PSF
cT=((sZ-1)/2)+1; %center

% Radius of the GK
R=(sZ-1)/2;

[jac,hG]=gauss2Djac(sigma_ga,R,height);

% return residual 
res=hB-hG;

res=res(:);
