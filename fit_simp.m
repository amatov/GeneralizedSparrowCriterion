function diff = fit_simp(x,X,Y)
% This function is called by lsqnonlin.
% x is a vector which contains the coefficients of the
% equation.  X and Y are the option data sets that were
% passed to lsqnonlin.

A=x(1);
B=x(2);
%C=x(3);
% D=x(4);
% E=x(5);

%diff = A + B.*exp(C.*-X)-Y; 
diff = A + B./X-Y; 



