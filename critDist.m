function [d,A1]=critDist(px)

% critDist displays the distance at which
% two GK merge as a function of its
% relative intensity
%
% SYNOPSIS   [d,A1]=critDist(px)
%
% INPUT      px       : initial distance
%
% OUTPUT     A1       : relative intensity
%            d        : critical disntance
% 
%
% DEPENDENCES   critDist uses { }
%               critDist is used by { }
%
% Alexandre Matov, January 7th, 2003


s=1.77; %d from 3.54
% 3.54 is the value for which two EQUAL kernels fuse
% for sigma 1.77(1.33PSF*1.33GaussFiltering)

if nargin == 0
    px = 5.2; % or 6
end

d=3.54:.01:px; % begin from the distance for ratio 1:1

P1=1/2*d-1/2*sqrt((d.^2-4*(s.^2)));
P2=1/2*d+1/2*sqrt((d.^2-4*(s.^2)));

A1=(P1.*exp(1/2*(s.^2+P2.*d-d.^2)./(s.^2)-1/2*(s.^2-P2.*d)./(s.^2)))./P2;
A2=(P2.*exp(1/2*(s.^2+P1.*d-d.^2)./(s.^2)-1/2*(s.^2-P1.*d)./(s.^2)))./P2;
    
figure,plot(A1,d,'g--') % plot in reverse order, so eventhough we can solve only A=f(d), we still can visualize d=f(A)
ylabel('Critical Distance')
xlabel('Relative Intencity')
    

