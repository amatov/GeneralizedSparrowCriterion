function compFits


RelInt=1:12;

r1=5.2-3.5.*exp(.5-RelInt); 
r2=(4.1814.*RelInt)./(2.2159+RelInt)+2.5613; 

plot(RelInt,r1)
hold on
plot(RelInt,r2,'r-')

