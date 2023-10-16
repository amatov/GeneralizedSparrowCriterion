function [jac,exy]=Gauss2D11center(sigma_ga,R,height,centr)

Debug=0;

if nargin == 0
    sigma_ga = 1.33;
    R = 15;
    height = 1;
    centr = [0 0];
end

%exy=zeros(2*R+1);
%NE ZNAM DALI DA BUHAM ZEROS PARVO

x=([-R:R]-centr(1))./sigma_ga;
y=([-R:R]-centr(2))./sigma_ga;

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
    figure,surf(exy);
    axis([0 2*R+1 0 2*R+1 0 height]);
    
    SIZE=size(exy);
    MAX=max(exy(:));
    t=1:size(jac,1);
    %figure,plot(t,jac);
    %figure,plot(jac);
    sigma_ga;
    %axis([0 2*R+1 0 2*R+1 0 height]);
end