clear all;clc;
x0=[0.1;0.05]; %起始點
A=[]; %線性不等式拘束條件的係數矩陣
b=[]; %線性不等式拘束條件的係數向量 AX <= b
Aeq=[]; %線性不等式拘束條件的係數向量
beq=[]; %線性等式拘束條件的係數向量 AeqX = beq
ub=[0.5;0.5]; %設計空間的 upper bounds
lb=[0.001;0.001]; %設計空間的 lower bounds
options = optimset ('display','off','Algorithm','sqp');%演算法的參數設定
[x,fval,exitflag]=fmincon(@(x)obj(x) ,x0,A,b,Aeq,beq,lb,ub,@(x)nonlcon(x),options);
%執行fmincon，輸出最佳解,r; 最佳目標函數值,fval; 收斂情形,exitflag
%obj為目標函數，nonlcon 為 (非線性) 拘束條件