function x=fit_ex(s,r)

% Define the data sets that you are trying to fit the function to

X=s;%15:-1:3;
%Y=2.0.*exp(5.0.*X)+3.0.*exp(2.5.*X)+1.5.*rand(size(X));

Y=r;%[3.9500    3.9500    3.8500    3.9000    4.0000    4.0500    4.1000    4.1000    4.1000    4.3000    4.2000    4.0500    4.7000];
%Y = [3.6 7.7 9.3 4.1 8.6 2.8 1.3 7.9 10.0 5.4];

% Initialize the coefficients of the function.
X0=[1 1]';

% Set an options file for LSQNONLIN to use the medium-scale algorithm 

% options = optimset('Largescale','off');

% Calculate the new coefficients using LSQNONLIN.
x=lsqnonlin('fit_simp',X0,[],[],[],X,Y);
%x=lsqcurvefit('fit_simp',X0,[],[],[],X,Y);

% Plot the original and experimental data.

%Y_new = x(1) + x(2).*exp(x(3).*-X); 
Y_new = x(1) + x(2)./X; 

figure,h=plot(X,Y,'r',X,Y_new,'b--');
xlabel('SNR','FontSize',30);
ylabel('distance','FontSize',30);
title('critical distance as function of signal-to-noise ratio','FontSize',30);
set(h,'LineWidth',3)
% AXIS(h,'FontSize',30)
% a=get(gca)
% set(a,'FontSize',30)
legend('simulation results','fitting curve')
axis([3 15 3.7 5]);
%axis([2 50 0 1.2]);