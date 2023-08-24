function [g,geq]=nonlcon(x)
[stress,Q_reduced] = sol_TenBarTruss(x(1), x(2))
yield_stress=250e+6;
disp_limit=0.02;
g(1)=max(abs(stress)-yield_stress)./1e+6;
g(2)=sqrt((Q_reduced(3))^2 + (Q_reduced(4))^2)-disp_limit;
geq=[];%等式拘束條件，在本問題中沒有等式拘束條件，故為空集合