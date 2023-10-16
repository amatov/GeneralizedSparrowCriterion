function fitDistCurve(d,A1)


cf=polyfit(A1,d,6);


ffd=cf(1).*A1 .^6+cf(2).*A1.^5+cf(3).*A1.^4+cf(4).*A1.^3+cf(5).*A1.^2+cf(6).*A1+cf(7);

figure,plot(A1,ffd,'r-')

pv = polyval(ffd,A1); 

hold on

plot(A1,pv,'g-')


% cf=polyfit(A1,d,2);
% 
% ffd2=cf(1).*A1.^2+cf(2).*A1+cf(3);
% 
% plot(A1,ffd2)