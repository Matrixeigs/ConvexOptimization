clc;clear;

Epsilon = 1e-3;
k=1;
p=0.01;
na=0;
u=na/p;

%%Pn  [Pbess Pg Pex] 4*3
Pn = sdpvar(4,3);
P_hold = [0 0 0 0];
bess_n = [0.3 0.32 0.25 0.35];
g_n    = [0.3269    0.0478    0.8488
    0.2744    0.0541    0.7098
    0.3000    0.0503    0.8208
    0.2930    0.0496    0.8222];
ex_n   = [0.02 0.04 0.02 0.06];
obj_Solved = zeros(1,1000);
err = zeros(1,1000);
err_s = zeros(4,1000);

Pl = [500 640 600 490];
Ps = [220 300 540 240];
Pw = [280 320 0 400];
Pgmm = [100 300;100 300;100 400;50 250];
%%Pgmm [Pgmin Pgmax]            4x2
Pbess = [200 230 230 200];
%%Pbess [Pbessmax]            4x1
Pn_Solved =  double(Pn);


while 1
    for n=1:4
        obj = bess_n(n) * Pn(n,1)^2 ;
        obj = obj + g_n(n,1) * Pn(n,2)^2 + g_n(n,2) * Pn(n,2) + g_n(n,3) ;
        obj = obj + ex_n(n) * Pn(n,3)^2 ;
        
        obj = obj + p/2*(norm( Pn(n,3) + sum(P_hold) - P_hold(n) + u,2)^2);
        
        cons = [];
        cons = cons + [Pn(1,1)-Pn(1,2)+Pn(1,3)-Pw(1)-Ps(1)+Pl(1) ==0];
        cons = cons + [Pn(2,1)-Pn(2,2)+Pn(2,3)-Pw(2)-Ps(2)+Pl(2) ==0];
        cons = cons + [Pn(3,1)-Pn(3,2)+Pn(3,3)-Pw(3)-Ps(3)+Pl(3) ==0];
        cons = cons + [Pn(4,1)-Pn(4,2)+Pn(4,3)-Pw(4)-Ps(4)+Pl(4) ==0];
        cons = cons + [Pgmm(:,1) <= Pn(:,2) <= Pgmm(:,2)];
        cons = cons + [-Pbess' <= Pn(:,1) <= Pbess'];
        
        
        ops = sdpsettings('solver','Cplex');
        sol = optimize(cons,obj,ops);
        
        err_s(k,n) = norm(p*(double(Pn(n,3))-P_hold(n)),2);
        P_hold(n) = double(Pn(n,3));
        Pn_Solved(n,:) =  double(Pn(n,:));
    end
    
    for n=1:4
        obj_Solved(k) = obj_Solved(k) + bess_n(n) * Pn_Solved(n,1)^2 + g_n(n,1) * Pn_Solved(n,2)^2 + g_n(n,2) * Pn_Solved(n,2) + g_n(n,3) + ex_n(n) * Pn_Solved(n,3)^2;
    end
    err(k) = norm(sum(P_hold));
    k = k + 1;
    u = u + sum(P_hold);
    if (err(k-1) <= Epsilon) && (err_s(k-1,1) <= Epsilon) && (err_s(k-1,2) <= Epsilon) && (err_s(k-1,3) <= Epsilon) && (err_s(k-1,4) <= Epsilon)
        break;
    end
end
figure
obj_Solved(k-1)                 %%目标函数
plot(obj_Solved(1:k-1))
figure
plot(err(1:k-1))                %%原始误差
figure
plot(err_s(1:k-1,:))            %%对偶误差
