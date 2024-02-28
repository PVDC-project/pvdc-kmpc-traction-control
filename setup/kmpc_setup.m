%% Koopman MPC controller setup
function [] = kmpc_setup(N,Ts,Rw,kappa_ref,compile_for_simulink)
%% system dynamics
load kmpc_data.mat Alift Blift Clift PX;

% add the integral state (e_int_dot = s - kappa_ref*R*w) using Euler method
Alift = [Alift(1:2,1:2), zeros(2,1), Alift(1:2,3:end);
     Ts -Ts*kappa_ref*Rw 1 zeros(1,size(Alift,1)-2);
     Alift(3:end,1:2) zeros(size(Alift,1)-2,1) Alift(3:end,3:end)];

Blift = [Blift(1:2,:)
         zeros(1,size(Blift,2));
         Blift(3:end,:)];

Clift = [eye(3), zeros(3,size(Clift,2)-2)];

ny = size(Clift,1);         % number of outputs
[nz, nu] = size(Blift);     % number of states and inputs in the Koopman model

% get dense form matrices (Nc=Np=N)
[F,Phi] = dense_prediction_matrices(Alift,Blift,Clift,N);

%% MPC parameters
% cost matrices
w_x1 = 1e-2;           % wheel slip velocity tracking error weight (unscaled)
w_x3 = 1e4;      % integral state weight (unscaled)
w_u = 1e-2;         % torque reduction weight (scaled)
w_x1 = 1e-6; w_x3 = 1e-6; w_u = 1e6;  % test MPC functionality
Q = diag([w_x1 0 w_x3]);    % state cost, no penalty on wheel speed
R = w_u;                    % input cost

%% define problem using YALMIP
u = sdpvar(nu,N);
z0 = sdpvar(nz,1);
T_ref = sdpvar;  % for online torque limiting
Y = F*z0 + Phi*u(:);  % dense-form prediction

constraints = [];
objective = 0;

for k = 1:N
    y = Y((k-1)*ny+1:k*ny);  % output from k-th prediction step
    s = (y(1)-PX.xoffset(1)) / PX.gain(1);  % unscale the first state
    w = (y(2)-PX.xoffset(2)) / PX.gain(2);  % unscale the second state
    es = s - kappa_ref*Rw*w;  % s tracking error (unscaled)
    ey = [es;0;y(3)];  % state error (ignore y(2), drive y(3) to zero)
    objective = objective + ey'*Q*ey;  % tracking cost
    eu = T_ref-u(k);  % input torque reduction
    objective = objective + eu'*R*eu;  % input torque reduction cost
    constraints = [constraints, 0 <= u(k) <= T_ref];  % input constraint
end

%% solver settings
algorithm = 'PDIP';  % https://forces.embotech.com/Documentation/solver_options/index.html#solve-methods
codeoptions = ForcesGetDefaultOptions('kmpc',algorithm,'double');
codeoptions.printlevel = 0;  % summary line after each solve
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

% controller = optimizer(constraints, objective, [], params, outputs);
optimizerFORCES(constraints, objective, codeoptions, params, outputs, parameter_names, output_names);

cd('../');

%% test call
problem{1} = lifting_function([0;10;0]);  % z0
problem{2} = 1;  % T_ref (scaled)

[~, exitflag, ~] = kmpc(problem);
assert(exitflag == 1, 'Test call of FORCESPRO solver failed');
end