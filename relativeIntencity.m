function relativeIntencity


r_plot=   [3.74 4.56 5.1 5.6 5.6 5.65 5.65 5.75  5.9  6.05 6.25 6.5 6.5 6.55 6.55 6.55 6.55 6.55 6.55 6.55 6.55 6.55 6.55 6.55 6.55 6.55 6.55 6.55];
ReIn_plot=[1    2    3   4   5   6    7    8     9    10   11   12  13  14   15   16   17   18   19   20   21   22   23   24   25   26   27   28];

coef_neg=fit_simp_neg(ReIn_plot,r_plot)

% coef=fit_ex(ReIn_plot,r_plot)