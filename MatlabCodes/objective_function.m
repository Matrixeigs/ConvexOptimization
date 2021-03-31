function obj = objective_function(x)
mpc = loadcase('case30');
nb = size(mpc.bus,1);
nl = size(mpc.branch,1);
ng = size(mpc.gen,1);


c = [zeros(nb,1);mpc.gencost(:,5)];
q = [zeros(nb,1);mpc.gencost(:,4)];
Q = 2*diag(q);

obj = c'*x + x'*Q*x/2;

end