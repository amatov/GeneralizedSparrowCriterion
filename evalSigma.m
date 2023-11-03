function sigma_ga=evalSigma(sigma_ga,height)

% evalSigma fits Bessel PSF with a GK
% (sigma_ga,height) = (1,1) default
% the height is fixed, for sigma_ga we give the initial value
%
% SYNOPSIS      sigma_ga=evalSigma(sigma_ga,height)
%
% 
% INPUT         sigma_ga    : sigma of the GK
%               height      : height of the GK
%
% OUTPUT        sigma_ga    : sigma after fitting
%
% DEPENDENCES   evalSigma uses { psf2D1, bessToGauss,
%                              lsqnonlin, gauss2Djac }
%               evalSigma is used by {  }
%
% Alexandre Matov, January 7th, 2003

options = optimset('Jacobian','on','Display','on');

if nargin == 0
    sigma_ga=1;
    height=1;
end

lb=0;
ub=10;

% SYNOPSIS PSF = psf2D(pixSize,NA,lambda)
hB=psf2D1(0.067,1.4,0.595);

sZ=size(hB,1);
cT=((sZ-1)/2)+1;
R=(sZ-1)/2;% Radius for the GK

sigma_ga=lsqnonlin(@bessToGauss,sigma_ga,lb,ub,options,hB,height);

[jac,hG]=gauss2Djac(sigma_ga,R,height);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEBUG

MaxValueCenterDataPSF=max(hB(:))
MaxValueCenterModelGK=max(hG(:))
HEIGHT_FIXED=height
SIGMA_PARAMETER=sigma_ga

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures

figure(1)
surf(hB);
title('PSF (pixSize,NA,lambda) = (0.067,1.4,0.595)');
axis([0 25 0 25 0 1.5]);

figure(2)
surf(hG);
title('fitting GK -> gauss2D11(sigma,height)');
axis([0 25 0 25 0 1.5]);

figure(3);
surf(hB-hG);
title('substruction PSF-GK; after optimization h - fixed; sigma - found using Matlab function LSQNONLIN');
axis([0 25 0 25 0 1.5]);
