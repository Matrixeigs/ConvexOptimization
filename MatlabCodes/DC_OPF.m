%% Test DCOPF


%% Parameters
% Data structure
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

mpc = loadcase('case30');

% % Define the variables
nb = size(mpc.bus,1);
nl = size(mpc.branch,1);
ng = size(mpc.gen,1);

% 1) Define the upper and lower boundary
Theta_lb = -180*ones(nb,1);
Theta_ub = 180*ones(nb,1);
Pg_lb = mpc.gen(:,PMIN);
Pg_ub = mpc.gen(:,PMAX);

x_lb=[Theta_lb;Pg_lb];
x_ub=[Theta_ub;Pg_ub];
%% 2) Define the objective functions
c = [zeros(nb,1);mpc.gencost(:,5)];
q = [zeros(nb,1);mpc.gencost(:,4)];
Q = 2*diag(q);
%% 3) Constraint set
[B, Bf, Pbusinj, Pfinj] = makeBdc(mpc.baseMVA,mpc.bus,mpc.branch);
Cg = sparse(mpc.gen(:, GEN_BUS), (1:ng)', mpc.gen(:, GEN_STATUS) > 0, nb, ng);
Aeq = [B,Cg];
beq = mpc.bus(:,PD);

A =[Bf, zeros(nl,ng);];
b = mpc.branch(:,RATE_A);

A=[A;-Bf, zeros(nl,ng)];
b=[b;mpc.branch(:,RATE_A)];

%% Solve the problem
 [x,fval,exitflag] = quadprog(Q,c,A,b,Aeq,beq,x_lb,x_ub);

 Theta = x(1:nb);
 Pg = x(nb+1:end);
 Pf = Bf * Theta;

