mpc.gend = [
    1  150 455 671 10.1  0.000299 80  120 400;
    2  150 455 574 10.2  0.000183 80  120 300;
    3  20  130 374 8.8   0.001126 130 130 105;
    4  20  130 374 8.8   0.001126 130 130 100;
    5  150 470 461 10.4  0.000205 80  120 90;
    6  135 460 630 10.1  0.000301 80  120 400;
    7  135 465 548 9.8   0.000364 80  120 350;
    8  60  300 227 11.2  0.000338 65  100 95;
    9  25  162 173 11.2  0.000807 60  100 105;
    10 25  160 175 10.7  0.001203 60  100 110;
    11 20  80  186 10.2  0.003586 80  80  60;
    12 20  80  230 9.9   0.005513 80  80  40;
    13 25  85  225 13.1  0.000371 80  80  30;
    14 15  55  309 12.1  0.001929 55  55  20;
    15 15  55  323 12.4  0.004447 55  55  20;
    ];

%% 火电机组的禁止区间
%%unit   Pi1L  Pi1U  Pi2L  Pi2U  Pi3L  Pi3U
mpc.prohz=[
    2    185  225  305  335  420  450;
    5    180  200  305  335  390  420;
    6    225  230  365  395  430  455;
    12   30   40   55   65   0    0  ;
    
    ];
%% B loss coefficients matrix网损b系数矩阵
mpc.Bij = [
    14  12   7   -1    -3   -1   -1   -1    -3     5    -3    -2    4    3   -1;
    12  15  13    0    -5   -2    0    1    -2    -4    -4     0    4   10   -2;
    7  13  76   -1   -13   -9   -1    0    -8   -12   -17     0  -26  111  -28;
    -1   0  -1   34    -7   -4   11   50    29    32   -11     0    1    1  -26;
    -3  -5 -13   -7    90   14   -3  -12   -10   -13     7    -2   -2  -24   -3;
    -1 -20  -9   -4    14   16    0   -6    -5    -8    11    -1   -2  -24    3;
    -1   0  -1   11    -3    0   15   17    15     9    -5     7    0   -2   -8;
    -1   1   0   50   -12   -6   17  168    82    79   -23   -36    1    5  -78;
    -3  -2  -8   29   -10   -5   15   82   129   116   -21   -25    7  -12  -72;
    -5  -4 -12   32   -13   -8    9   79   116   200   -27   -34    9  -11  -88;
    -30 -40 -17  -11     7   11   -5  -23   -21   -27   140     1    4  -38  168;
    -2   0   0    0    -2   -1    7  -36    25   -34     1    54   -1   -4   28;
    4   4 -26    1    -2   -2    0    1     7     9     4    -1  103 -101   28;
    3  10  111   1   -24  -17   -2    5   -12   -11   -38    -4 -101  578  -94;
    -1  -2  -28 -26    -3    3   -8  -78   -72   -88   168    28   28  -94  128;
    
    ];
mpc.B0j = [
    -1  -2  28  -1  -1  -3  -2  -2  6  39  -17  0  -32  67  -64;
    
    ];
%%发电机组启停状态
mpc.Pgk = [
    1  1  1  1
    1  1  1  1
    1  1  1  1
    1  1  1  1
    1  1  1  1
    1  1  1  1
    1  1  1  1
    1  1  1  1
    1  1  1  1
    1  1  1  1
    1  1  1  1
    1  1  1  1
    1  1  0  0
    1  1  0  0
    1  1  1  1
    ];
mpc.B00 = 50;
mpc.Bij = mpc.Bij * 1e-4;
mpc.B0j = mpc.B0j * 1e-4;
mpc.B00 = mpc.B00 * 1e-4;
mpc_case = mpc;


% [Unit,Pgmin,Pgmax,ai,bi,ci,URi,DRi,Pgo] = idx_gend;
Unit = mpc_case.gend(:,1);
Pgmin = mpc_case.gend(:,2);
Pgmax = mpc_case.gend(:,3);
ai = mpc_case.gend(:,4);
bi = mpc_case.gend(:,5);
ci = mpc_case.gend(:,6);
URi = mpc_case.gend(:,7);
DRi = mpc_case.gend(:,8);
Pgo = mpc_case.gend(:,9);
Bij = mpc_case.Bij;
B0j = mpc_case.B0j;
B00 = 0.005;
Bij = Bij * 0.0001;
B0j = B0j * 0.0001;
%unit   Pg1L  Pg1U  Pg2L  Pg2U  Pg3L  Pg3U
Unit2 = mpc_case.prohz(:,1);
Pg1L = mpc_case.prohz(:,2);
Pg1U = mpc_case.prohz(:,3);
Pg2L = mpc_case.prohz(:,4);
Pg2U = mpc_case.prohz(:,5);
Pg3L = mpc_case.prohz(:,6);
Pg3U = mpc_case.prohz(:,7);
Pgk = mpc_case.Pgk;

ng = size(mpc_case.gend,1);
n2 = size(mpc_case.prohz,1);
n3 = size(mpc_case.Bij,1);
n4 = size(mpc_case.Bij,2);

T = 4;
Ig = binvar(n2,4);
objsum = 0;
objex = 1:T;
Pg1 = zeros(15,T);


bigM = max(Pgmax);
for t=1:T
    Pg = sdpvar(ng,1);
    con = [];
    obj = 0;
    for a=1:ng
        if(Pgk(a,t)==1)
            obj = obj + ci(a) * Pg(a)^2 + bi(a) * Pg(a) + ai(a);
        end
    end
    %爬坡约束条件
    if(t>=2)
        for d = 1:ng
            if(Pgk(d,t)==1)
                con = con + [ Pg(d) - Pg1(d,t-1) <= URi(d)];
                con = con + [ Pg1(d,t-1) - Pg(d) <= DRi(d)];
            end
        end
    end
    %爬坡约束条件2（可与上式写在一起）
    if (t>=2)
        for e = 1:ng
            if(Pgk(e,t)==1)
                maxPg = max(Pgmin(e),Pg1(e,t-1) - DRi(e));
                minPg = min(Pgmax(e),Pg1(e,t-1) + URi(e));
                con = con + [ Pg(e) <= minPg];
                con = con + [ Pg(e) >= maxPg];
            end
        end
    end
    
    %最大最小值限制
    con = con + [Pg(:,1)<=Pgmax];
    con = con + [Pg(:,1)>=Pgmin];
    
    %     for b=1:ng
    %         con = con + [ Pg(b,t) >= Pgmin(b)];
    %         con = con + [ Pg(b,t) <= Pgmax(b)];
    %     end
    
    %     Ploss1 = 0;
    %     Ploss2 = 0;
    %     for i=1:n3
    %         for j=1:n4
    %             Ploss1 = Ploss1 + Pg(i,t)*Bij(i,j)*Pg(j,t);
    %         end
    %         Ploss2 = Ploss2 + B0j(i)*Pg(i,T);
    %     end
    
    %     Ploss = Ploss1 + Ploss2 + B00;
    
    %损耗约束
    con = con +[sum(Pg(:,1).*Pgk(:,t)) == 2630 + (Pg(:,1).*Pgk(:,t))'*Bij*(Pg(:,1).*Pgk(:,t)) + B0j*(Pg(:,1).*Pgk(:,t))+B00];
    
    %禁止区间
    for c = 1:n2
        %gend_id = find(Unit(:,1)==Unit2(c));
        %         con = con + [ Pg(Unit2(c),t) <=Pg1L(c)|Pg(Unit2(c),t) >=Pg1U(c)];
        %         con = con + [ Pg(Unit2(c),t) <=Pg2L(c)|Pg(Unit2(c),t) >=Pg2U(c)];
        %         con = con + [ Pg(Unit2(c),t) <=Pg3L(c)|Pg(Unit2(c),t) >=Pg3U(c)];
        if(Pgk(c,t)==1)
            if Pg1L(c)>Pgmin(Unit2(c))
                con = con +[ Pg(Unit2(c)) <= Ig(c,1) * Pg1L(c) + (1-Ig(c,1))*bigM];
                con = con +[ Pg(Unit2(c)) >= Ig(c,1) * Pgmin(Unit2(c))- (1-Ig(c,1))*bigM];
            end
            if Pg2L(c)>Pg1U(c)
                con = con +[ Pg(Unit2(c)) <= Ig(c,2) * Pg2L(c)+ (1-Ig(c,2))*bigM];
                con = con +[ Pg(Unit2(c)) >= Ig(c,2) * Pg1U(c)- (1-Ig(c,2))*bigM];
            end
            if Pg3L(c)>Pg2U(c)
                con = con +[ Pg(Unit2(c)) <= Ig(c,3) * Pg3L(c)+ (1-Ig(c,3))*bigM];
                con = con +[ Pg(Unit2(c)) >= Ig(c,3) * Pg2U(c)- (1-Ig(c,3))*bigM];
            end
            if Pgmax(Unit2(c))>Pg3U(c)
                con = con +[ Pg(Unit2(c)) <= Ig(c,4) * Pgmax(Unit2(c))+ (1-Ig(c,4))*bigM];
                con = con +[ Pg(Unit2(c)) >= Ig(c,4) * Pg3U(c)- (1-Ig(c,4))*bigM];
            end
            con = con +[sum(Ig(c,:))==1];
        end
    end
    
    
    
    % options = sdpsettings('solver','bnb');
    % sol = solvesdp(con,obj,options);
    sol = solvesdp(con,obj);
    obj = double(obj);
    objsum = objsum + obj; %对每次求解obj进行累加求和
    objex(t) = double(obj);
    Pg1(:,t) = double(Pg(:,1)).*Pgk(:,t);
end
%Pg = double(Pg).*Pgk;

