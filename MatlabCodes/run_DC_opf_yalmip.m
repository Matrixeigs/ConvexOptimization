function result = run_DC_opf_yalmip(mpc_case)
%% DC optimal power flow modelling based on YALMIP
if nargin<1
    mpc_case = loadcase('case30');
end
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
%%
Pd = mpc_case.bus(:,PD);
%%% Step 1: Define the number of variables %%%%
nb = size(mpc_case.bus,1);
nl = size(mpc_case.branch,1);
ng = size(mpc_case.gen,1);
% Define the variables
Pg = sdpvar(ng,1); % Output of generators
Pij = sdpvar(nl,1); % Power flow from the start bus to the end bus
Theta = sdpvar(nb,1); % Angle of each bus

%%% Step 2: Define the objective function %%%%
c0 = mpc_case.gencost(:,6);
c = mpc_case.gencost(:,5);
q = mpc_case.gencost(:,4);
obj = 0;
for g=1:ng
    obj = obj + q(g) * Pg(g)^2 + c(g) * Pg(g) + c0(g);
end

%%% Step 3: Constraint set %%%%
% 3.1)nodal power balance
cons = [];
for i = 1:nb
    % Find the index of generators at bus i
    gen_id = find(mpc_case.gen(:,GEN_BUS)==i);
    branch_id_f = find(mpc_case.branch(:,F_BUS)==i);
    branch_id_t = find(mpc_case.branch(:,T_BUS)==i);
    % KCL
    cons = cons + [sum(Pg(gen_id,1)) - Pd(i) == sum(Pij(branch_id_f,1)) - sum(Pij(branch_id_t,1)) ];
end
% 3.2) power flow constraints
cons = cons + [ Pij<= mpc_case.branch(:,RATE_A)];
cons = cons + [ Pij>= -mpc_case.branch(:,RATE_A)];

% 3.3) voltage angle constraints
cons = cons + [ Theta<= 180*ones(nb,1)];
cons = cons + [ Theta>= -180*ones(nb,1)];
% if you need to set one bus as the reference bus
cons = cons + [Theta(mpc_case.bus(:,BUS_TYPE)==REF)==0];
% 3.4) relation between power flow and angle
for i = 1:nl
   cons = cons + [Pij(i) == (Theta(mpc_case.branch(i,F_BUS))-Theta(mpc_case.branch(i,T_BUS)))/mpc_case.branch(i,BR_X)]; 
end


%%% Solving the problem

sol = solvesdp(cons,obj);

obj = double(obj);
Pg = double(Pg);
Pij = double(Pij);
Theta = double(Theta);


end

%% References:
% [ ]
% [ ]