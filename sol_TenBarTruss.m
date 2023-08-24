function [stress,Q_reduced] = sol_TenBarTruss(x1, x2)

    %定義各參數數值    
    nod_coor=[18.28 9.14;18.28 0;9.14 9.14;9.14 0;0 9.14;0 0];    %nodes coordinates
    num_nod=size(nod_coor,1);    %number of nodes
    ele_con=[3 5;1 3;4 6;2 4;3 4;1 2;4 5;3 6;2 3;1 4];    %elements' connection
    num_ele=size(ele_con,1);    %number of elements
    ele_dof=[5 6 9 10;1 2 5 6;7 8 11 12;3 4 7 8;5 6 7 8;1 2 3 4;7 8 9 10;5 6 11 12;3 4 5 6;1 2 7 8];    %elements' dof
    %A, E
    A(1)=x1^(2)*pi;
    A(2)=x1^(2)*pi;
    A(3)=x1^(2)*pi;
    A(4)=x1^(2)*pi;
    A(5)=x1^(2)*pi;
    A(6)=x1^(2)*pi;
    A(7)=x2^(2)*pi;
    A(8)=x2^(2)*pi;
    A(9)=x2^(2)*pi;
    A(10)=x2^(2)*pi;
    E=200*10^9;
        
    % 開一個空白的剛性矩陣 (stiffness matrix)
    % 計算 stiffness matrix
    K=zeros(2*num_nod);
    for e=1:num_ele
        L(e)=sqrt((nod_coor(ele_con(e,2),1)-nod_coor(ele_con(e,1),1))^2+(nod_coor(ele_con(e,2),2)-nod_coor(ele_con(e,1),2))^2);
        C=(nod_coor(ele_con(e,2),1)-nod_coor(ele_con(e,1),1))/L(e);
        S=(nod_coor(ele_con(e,2),2)-nod_coor(ele_con(e,1),2))/L(e);
        k=(A(e)*E/L(e)*[C*C C*S -C*C -C*S;C*S S*S -C*S -S*S;-C*C -C*S C*C C*S; -C*S -S*S C*S S*S]);
        
    ele_dof_vec=ele_dof(e,:);
    for i=1:4
        for j=1:4
        K(ele_dof_vec(1,i),ele_dof_vec(1,j))=K(ele_dof_vec(1,i),ele_dof_vec(1,j))+k(i,j);
        end
    end
    end

    % 建立力矩陣
    F=zeros(2*num_nod,1);
    F(4)=-10^7;
    F(8)=-10^7;
   
    % 建立空白位移矩陣
    Q=zeros(2*num_nod,1);
    % 計算位移量 Q_reduced 
    F_reduced=F(1:8);
    K_reduced=K(1:8,1:8);
    Q_reduced=Q(1:8);
    Q_reduced=K_reduced^(-1)*F_reduced;
    Q_reduced(12,1)=0;

    % 計算應力 (stress) 
    for e=1:num_ele
        L(e)=sqrt((nod_coor(ele_con(e,2),1)-nod_coor(ele_con(e,1),1))^2+(nod_coor(ele_con(e,2),2)-nod_coor(ele_con(e,1),2))^2);
        C=(nod_coor(ele_con(e,2),1)-nod_coor(ele_con(e,1),1))/L(e);
        S=(nod_coor(ele_con(e,2),2)-nod_coor(ele_con(e,1),2))/L(e);
        stress(e)=(E/L(e))*[-C -S C S]*Q_reduced((ele_dof(e,:))');
    end

    % 計算反作用力 R
    K_reaction=K(9:12,1:12);
    R=K_reaction*Q_reduced;
end

