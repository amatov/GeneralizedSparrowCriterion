function [jac,exy]=gauss2Djac(sigma_ga,R,height)

% gauss2Djac
% applies a two dimensional Gaussian filter
% calculates the Jacobian
%
% SYNOPSIS   [jac,exy]=gauss2Djac(sigma_ga,R,height)
%
% INPUT      R          :   radius of the GK
%            sigma_ga   :   sigma of the GK
%            height     :   height of the GK
%
% OUTPUT     IG         :   filtered image
%
% REMARKS       difference with function Gauss2D - the Gaussian mask
%               is not mormalized (so the sum is not equal to one);
%               gauss2Djac also calculates the Jacobian
%
%
% DEPENDENCES   gauss2Djac uses {  }
%               gauss2Djac is used by { evalSigma }
%
% Alexandre Matov, January 7th, 2003

Debug=0;

exy=zeros(2*R+1);
x=([-R:R])./sigma_ga;
y=([-R:R])./sigma_ga;

px=1/2*(x.^2);
py=1/2*(y.^2);

[p1,p2]=meshgrid(px,py);
p=p1+p2;

%ex=exp(-px);
%ey=exp(-py);

exy=height*exp(-p);
%exy=height*ex'*ey;

%jacox=exy.*2./sigma_ga.*p;
%jacox=jacox(:);
%jacoy=exy.*2./sigma_ga.*p;
%jacoy=jacoy(:);

jac=zeros(289,1);
%jac=[jacox;jacoy];
jac=-exy.*2./sigma_ga.*p;
jac=jac(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEBUG
if Debug == 1
    surf(exy);
    axis([0 2*R+1 0 2*R+1 0 height]);
    
    SIZE=size(exy);
    MAX=max(exy(:));
    t=1:size(jac,1);
    %figure,plot(t,jac);
    %figure,plot(jac);
    sigma_ga;
    %axis([0 2*R+1 0 2*R+1 0 height]);
end