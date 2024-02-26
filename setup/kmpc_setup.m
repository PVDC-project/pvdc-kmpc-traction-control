%% Koopman MPC controller setup
function [] = kmpc_setup(N,Ts,R,kappa_ref,compile_for_simulink)
%% system dynamics
load kmpc_data.mat Alift Blift Clift;

% add the integral state (e_int_dot = s - kappa_ref*R*w) using Euler method
Alift = [Alift(1:2,1:2), zeros(2,1), Alift(1:2,3:end);
     Ts -Ts*kappa_ref*R 1 zeros(1,size(Alift,1)-2);
     Alift(3:end,1:2) zeros(size(Alift,1)-2,1) Alift(3:end,3:end)];

Blift = [Blift(1:2,:)
         zeros(1,size(Blift,2));
         Blift(3:end,:)];

Clift = [eye(3), zeros(3,size(Clift,2)-2)];

ny = size(Clift,1);         % number of outputs
[nz, nu] = size(Blift);     % number of states and inputs in the Koopman model

% get dense form matrices (Nc=Np)
[F,Phi,~] = dense_prediction_matrices(Alift,Blift,Clift,N,N);

%% MPC parameters
% cost matrices
w_x1 = 1e3 / 1^2;       % wheel slip velocity tracking error weight
w_x3 = 1e5 / 1^2;       % integral state weight
w_u = 1e-2 / (250^2);   % torque reduction weight
Q = diag([w_x1 0 w_x3]);    % state cost
R = w_u;                    % input cost

%% define problem using YALMIP
u = sdpvar(nu,N);
z0 = sdpvar(nz,1);
T_ref = sdpvar;  % for online torque limiting
Y = F*z0 + Phi*u(:);  % dense-form prediction

constraints = [];
objective = 0;

for k = 1:N
    y = Y((k-1)*ny+1:k*ny);  % output from k-th step
    es = y(1) - kappa_ref*R*y(2);  % s tracking error
    ey = [es;y(2);y(3)];  % state error (drive y(2) and y(3) to zero)
    objective = objective + ey'*Q*ey;  % tracking cost
    objective = objective + u(:,k)'*R*u(:,k);  % input cost
    constraints = [constraints, 0 <= u(:,k) <= T_ref];  % input constraint
end

%% solver settings
algorithm = 'PDIP';  % https://forces.embotech.com/Documentation/solver_options/index.html#solve-methods
codeoptions = ForcesGetDefaultOptions('kmpc',algorithm,'double');
codeoptions.printlevel = 1;  % summary line after each solve
codeoptions.overwrite = 1;  % always overwrite the solver
% codeoptions.accuracy.ineq = 1e-4;  % infinity norm of residual for inequalities
% codeoptions.accuracy.eq = 1e-4;    % infinity norm of residual for equalities
% codeoptions.accuracy.mu = 1e-4;    % absolute duality gap
% options.condense = 1; % enable state-elimination

% Simulink block options
if ~compile_for_simulink
    codeoptions.BuildSimulinkBlock = 0;
end
codeoptions.showinfo = 1;  % https://forces.embotech.com/Documentation/solver_options/index.html#solver-info-in-simulink-block

params = {z0,T_ref};
outputs = {u,Y};
parameter_names = {'z0','T_ref'};
output_names = {'uopt','yopt'};

%% generate code
solver_dir = './codegen';
if ~exist(solver_dir,'dir')
    mkdir(solver_dir)
end
addpath(solver_dir);
cd(solver_dir);

controller = optimizer(constraints, objective, [], params, outputs);
optimizerFORCES(constraints, objective, codeoptions, params, outputs, parameter_names, output_names);

cd('../');

%% test call
problem.z0 = rand(nz,1);
problem.T_ref = 250;

[~, exitflag, ~] = kmpc(problem);
assert(exitflag == 1, 'Some issue with FORCESPRO solver');
end