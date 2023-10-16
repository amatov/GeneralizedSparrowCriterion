function [DER1,x]=twoGK(A,sigma,d,x)

% twoGK plots two overlapped Gaussian Kernels
% and the derivative of the function of the their sum
% 
%
% SYNOPSIS      [DER1,x]=twoGK(A,sigma,d,x)
%
% 
% INPUT         A     : the ratio between the two kernels
%               sigma : sigma
%               d     : distance between the centers
%               x     : x argument (ex.'x=-7:0.1:11;')
%
% OUTPUT        DER1  : the first derivative of their sum
%               x     : x argument 
%
% DEPENDENCES   twoGK uses {  }
%               twoGK is used by {  }
%
% Alexandre Matov, January 7th, 2003

if nargin==0
    A=1.5;
    sigma=.5;
    d=5.3;
    x=-7:0.1:11;
end

GK1=A*exp(-x.^2/2*sigma^2)/sqrt(2*pi*sigma^2);
GK2=exp(-(x-d).^2/2*sigma^2)/sqrt(2*pi*sigma^2);
GK3=exp(-(x-d+.6).^2/2*sigma^2)/sqrt(2*pi*sigma^2);

I1=GK1+GK2;
I2=GK1+GK3;
plot(x,I1)
hold on;
plot(x,GK1,'--');
plot(x,GK2,'--');
hold off;
figure,plot(x,I2)
hold on;
plot(x,GK1,'--');
plot(x,GK3,'--');
hold off;

DER1=-A*exp(-x.^2/(2*sigma^2)).*x/(sigma^3*sqrt(2*pi))-exp(-(x-d).^2/(2*sigma^2)).*(x-d)/(sigma^3*sqrt(2*pi));
DER1=DER1*2;

figure,plot(x,DER1);